# -*- coding: utf-8 -*-

"""Running a recipe: its components as ordinary calculations, then the combination.

Each component is an ordinary single point at the recipe's geometry, run through
its own model chemistry (``ORCA:DFT@...``) by seamm_exec's :class:`Evaluator`:
one evaluator per component over all the structures, so the structures of a
component run in parallel, as tasks or over MDI as the Evaluator chooses. Each
component's result is the code's ordinary result with its ordinary calculation
record and settings fingerprint, so the fragment store finds and keeps components,
never composites, and a component shared by two recipes is computed once.

Energies are combined in kJ/mol, as the Evaluator returns them. Every function of
the ``combine`` language is linear in the energies, so the units do not matter
except for the empirical constants, which recipes give in E_h and which are
converted here. A recipe of one component therefore gives exactly that
component's energy, with no unit conversion of it.
"""

from dataclasses import dataclass, field
import logging
from pathlib import Path

from .recipes import GEOMETRY, PART_FUNCTIONS, RecipeError, _names_used

logger = logging.getLogger(__name__)

#: kJ/mol per E_h (CODATA 2018, as seamm_util's pint registry).
HARTREE_TO_KJ = 2625.4996394798254


@dataclass
class CompositeResult:
    """One structure's composite energy.

    Attributes
    ----------
    key : str
    ok : bool
    E0 : float or None
        kJ/mol.
    results : dict
        Every value ``combine`` defines (E0 among them), kJ/mol.
    components : dict
        Per component: ``energy`` (kJ/mol), ``how`` ("computed", or what the
        consumer's lookup said), ``record`` (the calculation record, if
        computed) and ``model_chemistry``.
    reason : str or None
        Why it failed.
    """

    key: str
    ok: bool
    E0: float = None
    results: dict = field(default_factory=dict)
    components: dict = field(default_factory=dict)
    reason: str = None


def energy_sources(recipe):
    """The names ``combine`` takes energies from: components, and ``geometry``
    when the geometry level's energy is used. Components named only inside
    hf_part/corr_part count too."""
    used = set()
    for _, node in recipe.statements:
        plain, parts = _names_used(node)
        used |= plain | parts
    defined = {name for name, _ in recipe.statements}
    constants = set(recipe.data.get("empirical") or {})
    return (used - defined - constants) & (set(recipe.components) | {GEOMETRY})


def single_source(recipe):
    """The one energy E0 depends on, or None when it needs several (or the HF
    part of one, which a plain energy does not give)."""
    sources = energy_sources(recipe)
    if len(sources) != 1:
        return None
    for _, node in recipe.statements:
        for sub in _walk_calls(node):
            if sub in PART_FUNCTIONS:
                return None
    return next(iter(sources))


def _walk_calls(node):
    import ast

    for sub in ast.walk(node):
        if isinstance(sub, ast.Call) and isinstance(sub.func, ast.Name):
            yield sub.func.id


def component_levels(recipe):
    """{source: level spec} for every energy the recipe combines, the geometry
    level included when ``combine`` uses it. Phase 1 needs a complete ``level``
    for each; method + basis on the recipe's chosen code comes with the code
    choice (phase 2)."""
    levels = {}
    for source in sorted(energy_sources(recipe)):
        settings = recipe.geometry if source == GEOMETRY else recipe.components[source]
        level = settings.get("level")
        if level is None:
            raise RecipeError(
                f"The recipe {recipe.name}: '{source}' gives a method and basis but "
                "no 'level'; choosing the code for it is not implemented yet."
            )
        levels[source] = level
    return levels


def resolve_level(level, available=None):
    """The ``_model_chemistry`` wrapper of a level spec, as the Model Chemistry
    step resolves it (any basis for an offered owner/type/method)."""
    from model_chemistry_step.model_chemistry import (
        availability_problem,
        discover_model_chemistries,
        match_model_chemistry,
    )

    if available is None:
        available = discover_model_chemistries(periodic_only=False)
    mc = available.get(level) or match_model_chemistry(level, available)
    if mc is None:
        raise RecipeError(availability_problem(level, available, periodic=False))
    return mc


def evaluate_recipe(
    node,
    recipe,
    structures,
    *,
    directory=None,
    lookup=None,
    constants=None,
    available=None,
    evaluator_factory=None,
    name="Composite",
):
    """The composite energy of each structure.

    Parameters
    ----------
    node : seamm.Node
        The step running it (its flowchart, executor, directory).
    recipe : composite_step.recipes.Recipe
    structures : {str: configuration}
        The structures, by key (filesystem-safe keys: they name the tasks).
    directory : Path, optional
        Where the components' tasks go, ``<directory>/<component>``.
    lookup : callable, optional
        ``lookup(source, model_chemistry, key, configuration)`` -> ``(energy
        kJ/mol, how)`` or None: a value the consumer already has (from the
        fragment store or an earlier run), which is then not computed.
    constants : {str: float}, optional
        The recipe's empirical values, E_h; converted to kJ/mol here.
    available : dict, optional
        The discovered model chemistries (to avoid discovering them again).
    evaluator_factory : callable, optional
        Makes the evaluators, ``factory(node, mc, properties=, directory=,
        name=)``; seamm_exec.Evaluator by default (tests replace it).
    name : str

    Returns
    -------
    {str: CompositeResult}
        In the order of ``structures``.
    """
    if evaluator_factory is None:
        import seamm_exec

        evaluator_factory = seamm_exec.Evaluator

    levels = component_levels(recipe)
    mcs = {source: resolve_level(level, available) for source, level in levels.items()}

    known = set(recipe.data.get("empirical") or {})
    constants = dict(constants or {})
    missing = known - set(constants)
    if missing:
        raise RecipeError(
            f"The recipe {recipe.name} needs the empirical values "
            + ", ".join(sorted(missing))
        )
    constants_kj = {k: v * HARTREE_TO_KJ for k, v in constants.items()}

    out = {key: CompositeResult(key=key, ok=True) for key in structures}
    for source, mc in mcs.items():
        to_run = []
        for key, configuration in structures.items():
            if not out[key].ok:
                continue
            got = lookup(source, mc, key, configuration) if lookup else None
            if got is not None:
                energy, how = got
                out[key].components[source] = {
                    "energy": energy,
                    "how": how,
                    "record": None,
                    "model_chemistry": mc["level"],
                }
            else:
                to_run.append(key)
        if not to_run:
            continue
        where = Path(directory) / source if directory is not None else None
        with evaluator_factory(
            node,
            mc,
            properties=("energy",),
            directory=where,
            name=f"{name}_{source}",
        ) as evaluator:
            for key in to_run:
                evaluator.submit(structures[key], key=key)
            for result in evaluator.results():
                composite = out[result.key]
                if not result.ok or result.energy is None:
                    composite.ok = False
                    composite.reason = (
                        f"the {source} component ({mc['level']}) failed: "
                        f"{result.reason}"
                    )
                    continue
                record = (result.data or {}).get("calculation")
                composite.components[source] = {
                    "energy": result.energy,
                    "how": "computed",
                    "record": record if isinstance(record, dict) else None,
                    "model_chemistry": mc["level"],
                }
        for key in to_run:
            if out[key].ok and source not in out[key].components:
                out[key].ok = False
                out[key].reason = f"the {source} component returned no result"

    for composite in out.values():
        if not composite.ok:
            continue
        energies = {s: c["energy"] for s, c in composite.components.items()}
        try:
            composite.results = recipe.combine(energies, constants_kj)
        except (RecipeError, ArithmeticError) as e:
            composite.ok = False
            composite.reason = f"the combination failed: {e}"
            continue
        composite.E0 = composite.results["E0"]
    return out


# -- The single-component task route for any consumer of a model chemistry ---------


def _recipe_of(model_chemistry, data_path=None):
    from .recipes import RecipeLibrary

    options = model_chemistry.get("options") or {}
    name = options.get("recipe") or model_chemistry.get("method")
    version = options.get("recipe_version")
    extra = [Path(p) / "recipes" for p in (data_path or ())]
    return RecipeLibrary(extra_directories=extra).get(name, version)


def delegate(model_chemistry):
    """(recipe, source, the source's ``_model_chemistry``, its provider) for a
    recipe whose E0 depends on one energy; RecipeError otherwise."""
    recipe = _recipe_of(model_chemistry)
    source = single_source(recipe)
    if source is None:
        raise RecipeError(
            f"The composite {model_chemistry['level']} combines several "
            "calculations; it runs through composite_step.evaluate_recipe, not as "
            "one task per structure."
        )
    mc = resolve_level(component_levels(recipe)[source])
    return recipe, source, mc, provider_of(mc["step"])


_PROVIDERS = {}


def provider_of(step):
    """The plug-in class (``get_task``/``analyze_task``) of a step, by its
    stevedore name in ``org.molssi.seamm``, as model-chemistry discovery finds
    it."""
    if step not in _PROVIDERS:
        import stevedore

        manager = stevedore.ExtensionManager(
            namespace="org.molssi.seamm", invoke_on_load=False
        )
        for extension in manager:
            _PROVIDERS[extension.name] = extension.plugin
    if step not in _PROVIDERS:
        raise RecipeError(f"No step '{step}' is installed")
    return _PROVIDERS[step]


def delegated_get_task(configuration, model_chemistry, **kwargs):
    """``get_task`` of a one-energy recipe: the component's own task."""
    _, _, mc, provider = delegate(model_chemistry)
    return provider.get_task(configuration, mc, **kwargs)


def delegated_analyze_task(result, model_chemistry, configuration, **kwargs):
    """``analyze_task`` of a one-energy recipe: the component's own results, with
    the energy replaced by the recipe's E0 (kJ/mol) and the composite noted."""
    recipe, source, mc, provider = delegate(model_chemistry)
    if "gradients" in (kwargs.get("properties") or ()):
        from seamm_exec import AnalysisError

        raise AnalysisError(
            f"The composite {model_chemistry['level']} has no gradients yet."
        )
    data = dict(provider.analyze_task(result, mc, configuration, **kwargs))
    energy = data.get("energy")
    if energy is not None:
        results = recipe.combine({source: energy})
        data["component energy"] = energy
        data["energy"] = results["E0"]
        data["composite"] = {
            "recipe": recipe.name,
            "version": recipe.version,
            "fingerprint": recipe.fingerprint(),
            "component": source,
            "component model chemistry": mc["level"],
            "results": results,
        }
    return data
