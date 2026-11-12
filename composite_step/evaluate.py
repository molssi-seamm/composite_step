# -*- coding: utf-8 -*-

"""Running a recipe: its components as ordinary calculations, then the combination.

Each component is an ordinary single point at the recipe's geometry, run through
its own model chemistry (``ORCA:DFT@...``) by seamm_exec's :class:`Evaluator`:
one evaluator per component over all the structures, so the structures of a
component run in parallel, as tasks or over MDI as the Evaluator chooses. Each
component's result is the code's ordinary result with its ordinary calculation
record and settings fingerprint, so the fragment store finds and keeps components,
never composites, and a component shared by two recipes is computed once.

Energies are combined in kJ/mol, as the Evaluator returns them. Every result of a
``combine`` block scales like an energy (homogeneous of degree 1 in the energies,
checked when the recipe is read: see :data:`composite_step.recipes.DEGREES`), so
the result is the same in any unit; only the empirical constants, which recipes
give in E_h, are converted here. A recipe of one component therefore gives exactly that
component's energy, with no unit conversion of it.
"""

from dataclasses import dataclass, field
import json
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
    options=None,
):
    """The composite energy of each structure.

    The components run one after another, cheapest first (:func:`cost_key`),
    each over all the structures at once. A structure that fails a component is
    not sent to the later ones.

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
    options : {str: dict}, optional
        Per structure (by key), the options it is submitted with to every
        component (``charge``, ``multiplicity``, ``atom_indices``, ...).

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
    mcs = dict(sorted(mcs.items(), key=lambda item: cost_key(item[1], item[0])))
    options = options or {}
    if directory is None:
        directory = getattr(node, "directory", None)

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
                evaluator.submit(structures[key], key=key, options=options.get(key))
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


#: Rough relative cost of a model chemistry's type and method, for ordering the
#: components: a structure that fails a cheap component then never pays for an
#: expensive one.
_TYPE_COST = {"FF": 0, "SQM": 1, "MLFF": 1, "DFT": 2, "CMP": 2, "QC": 3}
_METHOD_COST = (
    ("CCSDT", 9),
    ("CCSD(T)", 7),
    ("QCISD(T)", 7),
    ("CCSD", 6),
    ("MP4", 5),
    ("MP3", 4),
    ("MP2", 3),
    ("HF", 1),
)


def cost_key(model_chemistry, name=""):
    """A sort key that puts cheap components first: the type, the method
    (DLPNO cheaper, F12 dearer), then the basis set's cardinal number."""
    level = str(model_chemistry.get("level") or "")
    owner_type, _, rest = level.partition("@")
    mtype = owner_type.split(":")[-1].upper()
    method, _, basis = rest.partition("/")
    method = method.upper()
    cost = 0.0
    for word, value in _METHOD_COST:
        if word in method:
            cost = value
            break
    if "DLPNO" in method:
        cost -= 2.5
    if "F12" in method:
        cost += 1.5
    cardinal = 0
    upper = basis.upper()
    for letter, n in (("D", 2), ("T", 3), ("Q", 4), ("5", 5), ("6", 6)):
        if (
            f"V{letter}Z" in upper
            or f"-{letter}ZVP" in upper
            or f"{letter}ZVP" in upper
        ):
            cardinal = n
    return (_TYPE_COST.get(mtype, 2), cost, cardinal, name)


def composite_record(recipe, model_chemistry, components):
    """The calculation record of a composite energy, in the form of a code's
    (:func:`seamm_util.fingerprint.calculation_record`): its fingerprint hashes the
    recipe's fingerprint and every component's settings fingerprint, so equal
    fingerprints mean the same recipe of the same calculations. The components'
    own records are kept under ``components``.

    Parameters
    ----------
    components : {str: dict}
        Per component, its ``record`` (the code's calculation record, or None
        when the energy came from a store or an earlier run without one) and
        ``model_chemistry``.
    """
    from seamm_util.fingerprint import calculation_record

    fingerprints = {
        source: (c.get("record") or {}).get("fingerprint")
        for source, c in components.items()
    }
    record = calculation_record(
        "SEAMM composite",
        str(recipe.version),
        {
            "recipe": recipe.fingerprint(),
            "components": json.dumps(fingerprints, sort_keys=True),
        },
        hashed_version=str(recipe.version),
    )
    geometry = {
        (c.get("record") or {}).get("geometry fingerprint") for c in components.values()
    } - {None}
    record["geometry fingerprint"] = geometry.pop() if len(geometry) == 1 else None
    record["recipe"] = recipe.name
    record["model chemistry"] = model_chemistry.get("level")
    record["components"] = {
        source: {
            "model chemistry": c.get("model_chemistry"),
            "fingerprint": fingerprints[source],
            "record": c.get("record"),
        }
        for source, c in components.items()
    }
    return record


def composite_data(recipe, model_chemistry, components, results):
    """``data`` as consumers get it on either route: ``energy`` (E0, kJ/mol),
    ``composite`` (the recipe, the components' energies and how they were
    obtained, every value ``combine`` defines) and ``calculation`` (the composite
    record, :func:`composite_record`)."""
    return {
        "energy": results["E0"],
        "composite": {
            "recipe": recipe.name,
            "version": recipe.version,
            "fingerprint": recipe.fingerprint(),
            "model chemistry": model_chemistry.get("level"),
            "components": {
                source: {
                    "energy": c["energy"],
                    "how": c.get("how", "computed"),
                    "model chemistry": c.get("model_chemistry"),
                }
                for source, c in components.items()
            },
            "results": dict(results),
        },
        "calculation": composite_record(recipe, model_chemistry, components),
    }


# -- The single-component task route for any consumer of a model chemistry ---------


def _recipe_of(model_chemistry, data_path=None):
    """The recipe of a ``SEAMM:CMP@`` model chemistry: from the library (shipped,
    and the ``recipes`` directories of ``data_path``), else from the file that
    model-chemistry discovery found it in (``recipe_path`` in its options), so a
    user's or a job's recipe is found wherever the model chemistry is run."""
    from .recipes import RecipeLibrary, load_recipe

    options = model_chemistry.get("options") or {}
    name = options.get("recipe") or model_chemistry.get("method")
    version = options.get("recipe_version")
    extra = [Path(p) / "recipes" for p in (data_path or ())]
    try:
        return RecipeLibrary(extra_directories=extra).get(name, version)
    except KeyError:
        path = options.get("recipe_path")
        if not path or not Path(path).is_file():
            raise
        recipe = load_recipe(path, source=options.get("recipe_source", "user"))
        if recipe.name != name or (version is not None and recipe.version != version):
            raise
        return recipe


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
        record = data.get("calculation")
        components = {
            source: {
                "energy": energy,
                "how": "computed",
                "record": record if isinstance(record, dict) else None,
                "model_chemistry": mc["level"],
            }
        }
        data.update(composite_data(recipe, model_chemistry, components, results))
    return data


# -- The evaluator of a multi-component recipe, for every consumer ------------------


class CompositeEvaluator:
    """Evaluates a composite model chemistry at many structures, with the
    interface of :class:`seamm_exec.Evaluator`, which returns it from
    ``Evaluator(...)`` through :meth:`CompositeStep.make_evaluator` when the recipe
    combines several calculations.

    Each component runs as an ordinary Evaluator of its own program, over all the
    submitted structures (:func:`evaluate_recipe`), with that program's own
    resources: a consumer's ``resources`` are not passed on, since an F12
    component needs other memory than a DFT one. The job's target and the task
    set's options are passed on.

    Energies only: a composite has no gradients or stress yet.

    The arguments ``node``, ``model_chemistry``, ``properties``, ``path``,
    ``directory``, ``target``, ``task_set_options``, ``resources`` and ``name`` are
    as for :class:`seamm_exec.Evaluator`. ``path`` (for tests only) is passed to
    every component's evaluator; otherwise each chooses its own.

    Parameters
    ----------
    recipe : composite_step.recipes.Recipe, optional
        Default: the model chemistry's.
    component_lookup : callable, optional
        ``lookup(source, model_chemistry, key, configuration)`` -> ``(energy
        kJ/mol, how)`` or None: a component energy the consumer already has
        (a fragment store, an earlier run), which is then not computed.
    available : dict, optional
        The discovered model chemistries.
    evaluator_factory : callable, optional
        For tests: makes the components' evaluators.
    """

    path = "batch"
    mdi_capable = False

    def __init__(
        self,
        node,
        model_chemistry=None,
        recipe=None,
        *,
        properties=("energy", "gradients"),
        path=None,
        directory=None,
        target=None,
        task_set_options=None,
        resources=None,
        name="Composite",
        component_lookup=None,
        available=None,
        evaluator_factory=None,
    ):
        if model_chemistry is None:
            model_chemistry = node.get_variable("_model_chemistry")
        self.node = node
        self.model_chemistry = model_chemistry
        self.options = model_chemistry.get("options") or {}
        self.recipe = recipe if recipe is not None else _recipe_of(model_chemistry)
        self.properties = tuple(properties)
        extra = sorted(set(self.properties) - {"energy"})
        if extra:
            raise ValueError(
                f"The composite method {model_chemistry.get('level')} gives energies "
                f"only; it cannot give {' or '.join(extra)}."
            )
        self._path = path
        self.directory = directory
        self.target = target
        self.task_set_options = dict(task_set_options or {})
        self.resources = resources  # not passed on, see above
        self.name = name
        self.component_lookup = component_lookup
        self.available = available
        self._factory = evaluator_factory
        self.constants = _constants(self.recipe)
        self._submitted = {}
        self._done = set()
        self._cancelled = False
        self._running = []

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()
        return False

    def close(self):
        pass

    def submit(self, configuration, key=None, *, options=None):
        """Add a structure; returns its key (unique, filesystem-safe)."""
        if key is None:
            key = f"s{len(self._submitted) + 1:06d}"
        key = str(key)
        if key in self._submitted:
            raise ValueError(f"Duplicate key '{key}'")
        self._submitted[key] = (configuration, dict(options or {}))
        return key

    def cancel(self):
        """Stop: the component evaluator running now is cancelled, and no other
        component starts."""
        self._cancelled = True
        for evaluator in self._running:
            cancel = getattr(evaluator, "cancel", None)
            if cancel is not None:
                cancel()

    @property
    def cancelled(self):
        return self._cancelled

    def expected_fingerprints(self, configuration, **kwargs):
        """None: a composite is never stored as such, only its components (whose
        fingerprints the components' programs predict)."""
        return None

    def _make_component_evaluator(self, node, mc, *, properties, directory, name):
        if self._cancelled:
            raise RuntimeError("cancelled")
        if self._factory is not None:
            evaluator = self._factory(
                node, mc, properties=properties, directory=directory, name=name
            )
        else:
            import seamm_exec

            evaluator = seamm_exec.Evaluator(
                node,
                mc,
                properties=properties,
                path=self._path,
                directory=directory,
                target=self.target,
                task_set_options=self.task_set_options,
                name=name,
            )
        self._running = [evaluator]
        return evaluator

    def results(self):
        """Compute what has been submitted; yield an
        :class:`seamm_exec.evaluator.EvaluatorResult` for each structure, once all
        its components are done."""
        from seamm_exec.evaluator import EvaluatorResult

        pending = {k: v for k, v in self._submitted.items() if k not in self._done}
        if not pending:
            return
        try:
            out = self._evaluate(pending)
        except RuntimeError:
            if not self._cancelled:
                raise
            out = {}
        for key in pending:
            if key not in out:
                self._done.add(key)
                yield EvaluatorResult(
                    key=key, ok=False, reason="cancelled", path="batch"
                )
        for key, composite in out.items():
            self._done.add(key)
            if not composite.ok:
                yield EvaluatorResult(
                    key=key, ok=False, reason=composite.reason, path="batch"
                )
                continue
            data = composite_data(
                self.recipe,
                self.model_chemistry,
                composite.components,
                composite.results,
            )
            yield EvaluatorResult(
                key=key, ok=True, energy=composite.E0, data=data, path="batch"
            )

    def _evaluate(self, pending):
        out = evaluate_recipe(
            self.node,
            self.recipe,
            {key: configuration for key, (configuration, _) in pending.items()},
            directory=self.directory,
            lookup=self.component_lookup,
            constants=self.constants,
            available=self.available,
            evaluator_factory=self._make_component_evaluator,
            name=self.name,
            options={key: options for key, (_, options) in pending.items()},
        )
        self._running = []
        return out


def _constants(recipe):
    """The recipe's empirical values (E_h) that are plain numbers. Terms that
    depend on the structure (higher-level corrections, spin-orbit tables) are not
    implemented yet, and a recipe that needs them is refused."""
    constants = {}
    for name, value in (recipe.data.get("empirical") or {}).items():
        if isinstance(value, (int, float)) and not isinstance(value, bool):
            constants[name] = float(value)
        elif isinstance(value, dict) and isinstance(value.get("value"), (int, float)):
            constants[name] = float(value["value"])
        else:
            raise RecipeError(
                f"The recipe {recipe.name}: the empirical term '{name}' depends on "
                "the structure, which is not implemented yet."
            )
    return constants
