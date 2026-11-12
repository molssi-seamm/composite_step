# -*- coding: utf-8 -*-

"""Composite-method recipes: the file format, its checks, the combination of the
component energies, the capability check and the thermal terms.

A recipe is a declarative YAML file that says exactly how an energy is made from
ordinary calculations: the geometry and frequency level, the single-point
components at that geometry, and an expression (``combine``) that combines them.
The result is itself a model chemistry, ``SEAMM:CMP@<recipe>``.

This module has no GUI and no SEAMM flowchart dependencies, so it can be used and
tested on its own. See ``NOTES_composite_methods.md`` (the design) for the reasons.

Units: a recipe's empirical constants are in hartree (E_h). Every result of
``combine`` must scale like an energy (degree 1 in the energies; checked when the
recipe is read), so the components may be combined in any one unit, E_h or kJ/mol,
with the constants in the same unit. The thermal terms likewise.
"""

import ast
import copy
from dataclasses import dataclass, field
import hashlib
import inspect
import json
import logging
import math
from pathlib import Path

import yaml

logger = logging.getLogger(__name__)

#: Characters a recipe name may not contain: they are the separators of the
#: model-chemistry grammar ``[owner:]type@method/basis|...``.
FORBIDDEN_IN_NAME = ":@/|"

#: The owner and type of SEAMM's recipes as model chemistries.
OWNER = "SEAMM"
TYPE = "CMP"

#: The name of the geometry level's own energy inside ``combine``.
GEOMETRY = "geometry"

#: The top-level fields of a recipe file.
FIELDS = {
    "recipe",
    "version",
    "description",
    "reference",
    "doi",
    "requires",
    "geometry",
    "frequencies",
    "frozen_core",
    "frozen_core_definitions",
    "components",
    "combine",
    "empirical",
    "atoms",
    "allow_mixed_reference",
    "uncertainty",
    "validation",
}
REQUIRED = {"recipe", "version", "geometry", "combine"}

#: Fields that describe a recipe without changing its energy; they are left out
#: of the fingerprint.
DESCRIPTIVE = {"description", "reference", "doi", "uncertainty", "validation"}

#: The fields of a component (one ordinary single point at the recipe's geometry).
COMPONENT_FIELDS = {
    "level",  # a complete model chemistry, e.g. ORCA:DFT@PBE0/def2-TZVP
    "method",  # or method + basis, run on the recipe's code
    "basis",
    "code",  # an explicit per-component override of the recipe's code
    "reference",  # RHF, UHF, ROHF
    "frozen_core",
    "density_fitting",
    "relativistic",
    "nucleus",  # point | finite, required with a relativistic Hamiltonian
    "f12",
    "sp_shells",
    "convergence",
    "options",  # anything else the code must be told, passed through as given
    "omit",  # parts deliberately left out, e.g. ["CABS singles"]
}


class RecipeError(ValueError):
    """A recipe file that is not valid, or a combination that cannot be made."""


# -- The combination: a small expression language ---------------------------------


@dataclass
class ComponentEnergy:
    """The energy of one component, with the parts a combination may need.

    Parameters
    ----------
    total : float
        The total energy.
    hf : float, optional
        The HF (reference) part of a correlated energy, in the same unit. Taken
        from ``parts["HF"]`` when not given.
    parts : {str: float}, optional
        The named parts of the energy (the program's ``energy parts``), in the same
        unit as ``total``.
    factors : {str: float}, optional
        Dimensionless factors, e.g. ``"(T) scale f"``.
    reference : str, optional
        The reference the correlation is built on: "RHF", "UHF" or "QRO".
    """

    total: float
    hf: float = None
    parts: dict = field(default_factory=dict)
    factors: dict = field(default_factory=dict)
    reference: str = None

    def __post_init__(self):
        if self.hf is None and "HF" in self.parts:
            self.hf = self.parts["HF"]

    @property
    def correlation(self):
        """The correlation energy alone: E - HF - CABS singles (an F12 run's CABS
        singles belong to neither the HF part nor the correlation)."""
        if self.hf is None:
            return None
        return self.total - self.hf - self.parts.get("CABS singles", 0.0)

    def part(self, name):
        """A named part: stored, or one of :data:`DERIVED_PARTS`."""
        if name in self.parts:
            return self.parts[name]
        if name == "HF" and self.hf is not None:
            return self.hf
        if name in DERIVED_PARTS:
            return DERIVED_PARTS[name](self)
        raise RecipeError(
            f"no part '{name}' (the parts are "
            f"{', '.join(sorted(self.parts)) or 'none'})"
        )


def _scaled_triples(triples):
    def scaled(c):
        if "(T) scale f" not in c.factors:
            raise RecipeError("no (T) scale factor f")
        return c.factors["(T) scale f"] * c.part(triples)

    return scaled


#: Parts computed at use from the stored ones, never stored themselves.
DERIVED_PARTS = {
    "CCSD-F12 correlation": lambda c: c.part("CCSD correlation")
    + c.part("F12 correction"),
    "DLPNO-CCSD-F12 correlation": lambda c: c.part("DLPNO-CCSD correlation")
    + c.part("F12 correction"),
    "HF+CABS": lambda c: c.part("HF") + c.part("CABS singles"),
    "(T*)": _scaled_triples("(T)"),
    "DLPNO-(T0*)": _scaled_triples("DLPNO-(T0)"),
}
#: The parts that are the reference's energy, whose references must agree
HF_TYPE_PARTS = {"HF", "HF+CABS"}


def cbs_power(e_x, e_y, x, y, alpha=3.0):
    """Two-point inverse-power extrapolation (Helgaker 1997, Halkier 1998):
    E(n) = E_CBS + A n^-alpha, solved from cardinal numbers x and y."""
    if x == y:
        raise RecipeError("cbs_power needs two different cardinal numbers")
    xa = x**alpha
    ya = y**alpha
    return (xa * e_x - ya * e_y) / (xa - ya)


def hf_extrapolate_exp(e_x, e_y, n, alpha):
    """Two-point extrapolation of the HF energy in the form
    E(n) = E_CBS + B (n + 1) exp(-alpha sqrt(n)), from cardinal numbers n and n + 1.

    G4 uses alpha = 1.63 (Curtiss 2007); Karton and Martin (2006) use alpha = 9.
    """
    f_x = (n + 1) * math.exp(-alpha * math.sqrt(n))
    f_y = (n + 2) * math.exp(-alpha * math.sqrt(n + 1))
    if f_x == f_y:
        raise RecipeError("hf_extrapolate_exp: the two terms are equal")
    return (f_x * e_y - f_y * e_x) / (f_x - f_y)


def schwenke(e_x, e_y, factor):
    """Schwenke's (2005) two-point extrapolation: E_CBS = e_x + F (e_y - e_x)."""
    return e_x + factor * (e_y - e_x)


def scale(e, factor):
    """An energy times a factor."""
    return e * factor


#: The functions ``combine`` may call, with the number of positional arguments
#: they take (min, max).
FUNCTIONS = {
    "cbs_power": (cbs_power, 4, 5),
    "hf_extrapolate_exp": (hf_extrapolate_exp, 4, 4),
    "schwenke": (schwenke, 3, 3),
    "scale": (scale, 2, 2),
}
#: The degree of each function's arguments and result in the energies: 1 for an
#: energy, 0 for a pure number (a cardinal number, an exponent, a factor). Every
#: function is homogeneous of degree 1 in its energy arguments, so a combination
#: gives the same result in any energy unit. A new function must be too: a fitted
#: exponent, for example, enters only through a ratio of energy differences
#: (degree 0), as in Feller's three-point form.
DEGREES = {
    "cbs_power": ((1, 1, 0, 0, 0), 1),
    "hf_extrapolate_exp": ((1, 1, 0, 0), 1),
    "schwenke": ((1, 1, 0), 1),
    "scale": ((1, 0), 1),
    "hf_part": ((1,), 1),
    "corr_part": ((1,), 1),
    "part": ((1, 0), 1),
}
#: Functions of a component itself, not of its energy.
PART_FUNCTIONS = {"hf_part", "corr_part", "part"}

_BINARY = {
    ast.Add: lambda a, b: a + b,
    ast.Sub: lambda a, b: a - b,
    ast.Mult: lambda a, b: a * b,
    ast.Div: lambda a, b: a / b,
}
_UNARY = {ast.USub: lambda a: -a, ast.UAdd: lambda a: a}


def parse_combine(text):
    """Parse a ``combine`` block into a list of (name, expression AST).

    Each line is ``name = expression``. An expression may use numbers, the names
    of components and of earlier results, + - * / and parentheses, and the
    functions in :data:`FUNCTIONS` and :data:`PART_FUNCTIONS`. Nothing else is
    accepted: the text is parsed with :mod:`ast` and never given to ``eval``.
    """
    if not isinstance(text, str) or not text.strip():
        raise RecipeError("'combine' must be a non-empty text")
    try:
        tree = ast.parse(text, mode="exec")
    except SyntaxError as e:
        raise RecipeError(f"'combine' is not valid: {e.msg} (line {e.lineno})")
    statements = []
    for node in tree.body:
        if (
            not isinstance(node, ast.Assign)
            or len(node.targets) != 1
            or not isinstance(node.targets[0], ast.Name)
        ):
            raise RecipeError(
                f"'combine' line {node.lineno}: each line must be 'name = expression'"
            )
        _check_expression(node.value, node.lineno)
        statements.append((node.targets[0].id, node.value))
    return statements


def _check_expression(node, line):
    """Refuse anything in an expression that is not in the whitelist."""
    if isinstance(node, ast.Constant):
        if isinstance(node.value, bool) or not isinstance(node.value, (int, float)):
            raise RecipeError(f"'combine' line {line}: only numbers are allowed")
    elif isinstance(node, ast.Name):
        pass
    elif isinstance(node, ast.BinOp):
        if type(node.op) not in _BINARY:
            raise RecipeError(f"'combine' line {line}: only + - * / are allowed")
        _check_expression(node.left, line)
        _check_expression(node.right, line)
    elif isinstance(node, ast.UnaryOp):
        if type(node.op) not in _UNARY:
            raise RecipeError(f"'combine' line {line}: only unary + and - are allowed")
        _check_expression(node.operand, line)
    elif isinstance(node, ast.Call):
        if not isinstance(node.func, ast.Name):
            raise RecipeError(f"'combine' line {line}: only plain function calls")
        name = node.func.id
        if name == "part":
            if (
                len(node.args) != 2
                or node.keywords
                or not isinstance(node.args[0], ast.Name)
                or not isinstance(node.args[1], ast.Constant)
                or not isinstance(node.args[1].value, str)
            ):
                raise RecipeError(
                    f"'combine' line {line}: part() takes a component and the name "
                    'of a part, e.g. part(ccsdt, "(T)")'
                )
            return
        if name in PART_FUNCTIONS:
            if (
                len(node.args) != 1
                or node.keywords
                or not isinstance(node.args[0], ast.Name)
            ):
                raise RecipeError(
                    f"'combine' line {line}: {name}() takes the name of one component"
                )
            return
        if name not in FUNCTIONS:
            raise RecipeError(
                f"'combine' line {line}: unknown function '{name}'; the functions "
                f"are {', '.join(sorted([*FUNCTIONS, *PART_FUNCTIONS]))}"
            )
        _, lo, hi = FUNCTIONS[name]
        n = len(node.args) + len(node.keywords)
        if not lo <= n <= hi:
            raise RecipeError(
                f"'combine' line {line}: {name}() takes {lo}"
                + (f" to {hi}" if hi != lo else "")
                + f" arguments, not {n}"
            )
        for arg in node.args:
            _check_expression(arg, line)
        for kw in node.keywords:
            if kw.arg is None:
                raise RecipeError(f"'combine' line {line}: no ** arguments")
            _check_expression(kw.value, line)
    else:
        raise RecipeError(
            f"'combine' line {line}: '{type(node).__name__}' is not allowed"
        )


def _names_used(node):
    """The names an expression refers to, with those used as a component itself
    (the argument of hf_part/corr_part) separately. Function names are not names."""
    plain, parts, skip = set(), set(), set()
    for sub in ast.walk(node):
        if isinstance(sub, ast.Call):
            skip.add(id(sub.func))
            if sub.func.id in PART_FUNCTIONS:
                parts.add(sub.args[0].id)
                skip.add(id(sub.args[0]))
    for sub in ast.walk(node):
        if isinstance(sub, ast.Name) and id(sub) not in skip:
            plain.add(sub.id)
    return plain, parts


def parts_used(statements):
    """{component: set of what is used of it}: "hf_part", "corr_part" and the
    names given to part()."""
    used = {}
    for _, node in statements:
        for sub in ast.walk(node):
            if (
                isinstance(sub, ast.Call)
                and isinstance(sub.func, ast.Name)
                and sub.func.id in PART_FUNCTIONS
            ):
                which = sub.args[0].id
                what = sub.func.id if sub.func.id != "part" else sub.args[1].value
                used.setdefault(which, set()).add(what)
    return used


def evaluate(statements, components, constants=None, allow_mixed_reference=False):
    """Evaluate parsed ``combine`` statements.

    Parameters
    ----------
    statements : [(str, ast.AST)]
        From :func:`parse_combine`.
    components : {str: ComponentEnergy or float}
        The energies of the components (and of the geometry level), E_h.
    constants : {str: float}, optional
        Other named values: empirical terms, spin-orbit, E_h.
    allow_mixed_reference : bool
        Allow HF-type parts of components with different references (RHF, UHF,
        QRO) in one combination; refused otherwise.

    Returns
    -------
    {str: float}
        Every result the statements define, in order; ``E0`` among them.
    """
    comps = {
        k: v if isinstance(v, ComponentEnergy) else ComponentEnergy(float(v))
        for k, v in components.items()
    }
    values = {k: c.total for k, c in comps.items()}
    values.update(constants or {})
    results = {}
    references = {}

    def ev(node):
        if isinstance(node, ast.Constant):
            return float(node.value)
        if isinstance(node, ast.Name):
            if node.id not in values:
                raise RecipeError(f"'combine': no value for '{node.id}'")
            return values[node.id]
        if isinstance(node, ast.BinOp):
            return _BINARY[type(node.op)](ev(node.left), ev(node.right))
        if isinstance(node, ast.UnaryOp):
            return _UNARY[type(node.op)](ev(node.operand))
        if isinstance(node, ast.Call):
            name = node.func.id
            if name in PART_FUNCTIONS:
                which = node.args[0].id
                if which not in comps:
                    raise RecipeError(f"'combine': {name}() of unknown '{which}'")
                component = comps[which]
                if name == "part":
                    wanted = node.args[1].value
                    try:
                        value = component.part(wanted)
                    except RecipeError as e:
                        raise RecipeError(f"'combine': the component '{which}': {e}")
                    if wanted in HF_TYPE_PARTS:
                        references[which] = component.reference
                    return value
                value = component.hf if name == "hf_part" else component.correlation
                if value is None:
                    raise RecipeError(
                        f"'combine': the component '{which}' has no HF part"
                    )
                if name == "hf_part":
                    references[which] = component.reference
                return value
            function = FUNCTIONS[name][0]
            args = [ev(a) for a in node.args]
            kwargs = {kw.arg: ev(kw.value) for kw in node.keywords}
            return function(*args, **kwargs)
        raise RecipeError(f"'combine': cannot evaluate {type(node).__name__}")

    for name, node in statements:
        value = ev(node)
        results[name] = value
        values[name] = value
    if "E0" not in results:
        raise RecipeError("'combine' must define E0")
    kinds = {r for r in references.values() if r is not None}
    if len(kinds) > 1 and not allow_mixed_reference:
        raise RecipeError(
            "'combine' uses the HF parts of components with different references ("
            + ", ".join(f"{k} {references[k]}" for k in sorted(references))
            + "); say 'allow_mixed_reference: true' in the recipe if that is meant"
        )
    return results


# -- The recipe ---------------------------------------------------------------------


@dataclass
class Recipe:
    """A composite method read from a recipe file.

    Use :func:`load_recipe` or :meth:`from_dict`, which check it.
    """

    data: dict
    path: Path = None
    source: str = "shipped"  # shipped, user or job
    statements: list = field(default_factory=list, repr=False)

    @classmethod
    def from_dict(cls, data, path=None, source="shipped"):
        data = copy.deepcopy(data)
        statements = validate(data)
        return cls(data=data, path=path, source=source, statements=statements)

    # The fields
    @property
    def name(self):
        return self.data["recipe"]

    @property
    def version(self):
        return self.data["version"]

    @property
    def description(self):
        return self.data.get("description", "")

    @property
    def model_chemistry(self):
        """The level spec of this recipe, ``SEAMM:CMP@<name>``."""
        return f"{OWNER}:{TYPE}@{self.name}"

    @property
    def geometry(self):
        return self.data["geometry"]

    @property
    def components(self):
        return self.data.get("components") or {}

    @property
    def zpe_scale(self):
        """(factor, source) for the ZPE at the frequency level."""
        freq = self.data.get("frequencies") or {}
        return float(freq.get("scale_zpe", 1.0)), freq.get("scale_source", "")

    @property
    def requires(self):
        return list(self.data.get("requires") or [])

    def combine(self, components, constants=None):
        """E0 and the other results from the components' energies (E_h).

        ``components`` must include the geometry level's energy as ``geometry``
        when the combination uses it.
        """
        return evaluate(
            self.statements,
            components,
            constants,
            allow_mixed_reference=bool(self.data.get("allow_mixed_reference")),
        )

    @property
    def uses_default_frozen_core(self):
        """Whether a component leaves the frozen core to the program's default."""
        return any(c.get("frozen_core") == "default" for c in self.components.values())

    @property
    def needs_parts(self):
        """The components whose parts ``combine`` uses (hf_part, corr_part,
        part)."""
        return set(parts_used(self.statements))

    def thermal(self, E0, geometry):
        """H(0 K) and H(T) from E0 and the frequency level's values.

        Parameters
        ----------
        E0 : float
            The recipe's electronic energy.
        geometry : dict
            ``E`` (the frequency level's electronic energy), ``ZPE`` (unscaled)
            and optionally ``H`` (its enthalpy at T: E + ZPE + thermal), in the
            same units as E0.

        Returns
        -------
        dict
            ``E0``, ``ZPE`` (unscaled), ``ZPE scaled``, ``H0`` and, with ``H``,
            ``H``. H0 = E0 + s ZPE and H = E0 + s ZPE + [H - E - ZPE] at the
            frequency level: only the ZPE is scaled.
        """
        s, _ = self.zpe_scale
        ZPE = geometry["ZPE"]
        result = {"E0": E0, "ZPE": ZPE, "ZPE scaled": s * ZPE}
        # The order of the additions follows the thermochemistry step's, so that
        # the results agree with it to the last bit (the scaling is applied as a
        # shift of the unscaled sums there)
        shift = (s - 1.0) * ZPE
        result["H0"] = E0 + ZPE
        if s != 1.0:
            result["H0"] += shift
        if geometry.get("H") is not None:
            result["H"] = E0 + (geometry["H"] - geometry["E"])
            if s != 1.0:
                result["H"] += shift
        return result

    def fingerprint(self, component_fingerprints=None):
        """A hash of the normalized recipe (its descriptive fields left out),
        and of the components' own fingerprints when given."""
        normal = {k: v for k, v in self.data.items() if k not in DESCRIPTIVE}
        payload = {"recipe": normal}
        if component_fingerprints:
            payload["components"] = dict(sorted(component_fingerprints.items()))
        text = json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str)
        return "sha256:" + hashlib.sha256(text.encode()).hexdigest()


def validate(data):
    """Check a recipe dict; return its parsed ``combine`` statements."""
    if not isinstance(data, dict):
        raise RecipeError("A recipe must be a mapping")
    unknown = set(data) - FIELDS
    if unknown:
        raise RecipeError(f"Unknown recipe fields: {', '.join(sorted(unknown))}")
    missing = REQUIRED - set(data)
    if missing:
        raise RecipeError(f"Missing recipe fields: {', '.join(sorted(missing))}")

    name = data["recipe"]
    if not isinstance(name, str) or not name.strip():
        raise RecipeError("'recipe' must be a non-empty name")
    bad = [c for c in FORBIDDEN_IN_NAME if c in name]
    if bad or name != name.strip():
        raise RecipeError(
            f"The recipe name '{name}' may not contain {' '.join(FORBIDDEN_IN_NAME)} "
            "or leading or trailing spaces"
        )
    version = data["version"]
    if isinstance(version, bool) or not isinstance(version, int) or version < 1:
        raise RecipeError("'version' must be a positive integer")

    _check_level("geometry", data["geometry"], allow_scale=False)

    freq = data.get("frequencies") or {}
    if not isinstance(freq, dict):
        raise RecipeError("'frequencies' must be a mapping")
    unknown = set(freq) - {"level", "method", "basis", "scale_zpe", "scale_source"}
    if unknown:
        raise RecipeError(f"Unknown 'frequencies' fields: {', '.join(sorted(unknown))}")
    if "scale_zpe" in freq:
        s = freq["scale_zpe"]
        if isinstance(s, bool) or not isinstance(s, (int, float)) or not 0 < s <= 1.5:
            raise RecipeError("'frequencies: scale_zpe' must be a number in (0, 1.5]")
        if not freq.get("scale_source"):
            raise RecipeError("'frequencies: scale_zpe' needs its 'scale_source'")

    components = data.get("components") or {}
    if not isinstance(components, dict):
        raise RecipeError("'components' must be a mapping of name: settings")
    for cname, comp in components.items():
        if cname == GEOMETRY:
            raise RecipeError(f"'{GEOMETRY}' is reserved for the geometry level")
        if not str(cname).isidentifier():
            raise RecipeError(f"The component name '{cname}' must be an identifier")
        _check_level(f"component '{cname}'", comp)

    empirical = data.get("empirical") or {}
    if not isinstance(empirical, dict):
        raise RecipeError("'empirical' must be a mapping")

    statements = parse_combine(data["combine"])
    known = set(components) | {GEOMETRY} | set(empirical)
    defined = set()
    # Energies (components, the geometry level, empirical terms) have degree 1
    degrees = {name: 1 for name in known}
    for target, node in statements:
        if target in known:
            raise RecipeError(f"'combine' redefines '{target}'")
        plain, parts = _names_used(node)
        undefined = (plain - known - defined) | (parts - set(components))
        if undefined:
            raise RecipeError(
                f"'combine' line for '{target}' uses unknown names: "
                + ", ".join(sorted(undefined))
            )
        degrees[target] = _degree(node, degrees, target)
        defined.add(target)
    if "E0" not in defined:
        raise RecipeError("'combine' must define E0")
    for target, _ in statements:
        if degrees[target] != 1:
            raise RecipeError(
                f"'combine': '{target}' is not an energy (it scales as energy^"
                f"{degrees[target]}); every result must scale like an energy, so "
                "that it is the same in any unit"
            )

    _check_core_and_relativity(data, components)
    _check_subtracted_pairs(statements, components)

    # An F12 component split into hf_part + corr_part must account for its CABS
    # singles, which belong to neither
    for which, used in parts_used(statements).items():
        comp = components.get(which)
        if comp is None or not {"hf_part", "corr_part"} <= used:
            continue
        level = str(comp.get("level") or comp.get("method") or "")
        if "F12" not in level.upper():
            continue
        omitted = comp.get("omit") or []
        if "CABS singles" not in used and "CABS singles" not in omitted:
            raise RecipeError(
                f"The component '{which}' is an F12 calculation split with hf_part "
                "and corr_part, which leave out its CABS singles: add "
                f'part({which}, "CABS singles"), or say omit: ["CABS singles"] '
                "for the component"
            )

    requires = data.get("requires") or []
    if not isinstance(requires, list) or not all(isinstance(r, str) for r in requires):
        raise RecipeError("'requires' must be a list of capabilities")
    return statements


RELATIVISTIC = ("none", "DKH2", "X2C")
NUCLEUS = ("point", "finite")
FROZEN_CORE_SPECIAL = ("none", "default")


def _check_core_and_relativity(data, components):
    """The frozen-core sets (orbitals per element) and each component's frozen
    core, relativistic Hamiltonian and nucleus model."""
    sets = data.get("frozen_core_definitions") or {}
    if not isinstance(sets, dict):
        raise RecipeError("'frozen_core_definitions' must be a mapping of name: set")
    for name, counts in sets.items():
        if name in FROZEN_CORE_SPECIAL:
            raise RecipeError(f"'{name}' is reserved and cannot name a frozen-core set")
        if not isinstance(counts, dict) or not all(
            isinstance(n, int) and not isinstance(n, bool) and n >= 0
            for n in counts.values()
        ):
            raise RecipeError(
                f"The frozen-core set '{name}' must give a number of orbitals (an "
                "integer, 0 or more) for each element"
            )
    for cname, comp in components.items():
        fc = comp.get("frozen_core")
        if fc is not None and fc not in FROZEN_CORE_SPECIAL and fc not in sets:
            raise RecipeError(
                f"The component '{cname}': no frozen-core set '{fc}' (the sets are "
                f"{', '.join(sorted(sets)) or 'none'}; or 'none', 'default')"
            )
        rel = comp.get("relativistic", "none")
        if rel not in RELATIVISTIC:
            raise RecipeError(
                f"The component '{cname}': relativistic must be one of "
                f"{', '.join(RELATIVISTIC)}, not {rel!r}"
            )
        nucleus = comp.get("nucleus")
        if rel != "none" and nucleus not in NUCLEUS:
            raise RecipeError(
                f"The component '{cname}' uses {rel}: give its nucleus model "
                f"explicitly, nucleus: {' or '.join(NUCLEUS)}"
            )
        if rel == "none" and nucleus is not None:
            raise RecipeError(
                f"The component '{cname}' gives a nucleus model without a "
                "relativistic Hamiltonian"
            )


def _level_parts(comp):
    """(method, basis) of a component: from its level, or its method and basis."""
    level = comp.get("level")
    if not level:
        return str(comp.get("method") or ""), str(comp.get("basis") or "")
    rest = str(level).split("@", 1)[-1]
    method, _, basis = rest.partition("/")
    return method, basis


def _differences(a, b):
    """The kinds of difference between two components: method, basis, frozen
    core, relativistic (Hamiltonian and nucleus), other settings."""
    kinds = set()
    (ma, ba), (mb, bb) = _level_parts(a), _level_parts(b)
    if ma != mb:
        kinds.add("method")
    if ba != bb:
        kinds.add("basis")
    if a.get("frozen_core") != b.get("frozen_core"):
        kinds.add("frozen core")
    if (a.get("relativistic", "none"), a.get("nucleus")) != (
        b.get("relativistic", "none"),
        b.get("nucleus"),
    ):
        kinds.add("relativistic")
    owner_a = str(a.get("level") or a.get("code") or "").split("@", 1)[0]
    owner_b = str(b.get("level") or b.get("code") or "").split("@", 1)[0]
    if owner_a != owner_b:
        kinds.add("program")
    skip = {
        "level",
        "code",
        "frozen_core",
        "relativistic",
        "nucleus",
        "method",
        "basis",
    }
    for key in (set(a) | set(b)) - skip:
        if a.get(key) != b.get(key):
            kinds.add(key)
    return kinds


#: The ways the two components of a difference may differ: one thing at a time
#: (a relativistic correction changes the Hamiltonian and, with it, to a
#: recontracted basis)
ALLOWED_PAIR_DIFFERENCES = (
    {"method"},
    {"basis"},
    {"frozen core"},
    {"relativistic"},
    {"relativistic", "basis"},
)


def _check_subtracted_pairs(statements, components):
    """A difference of two components (``a - b``) must change one thing only, so
    that it is the correction it is meant to be."""
    for _, node in statements:
        for sub in ast.walk(node):
            if not (
                isinstance(sub, ast.BinOp)
                and isinstance(sub.op, ast.Sub)
                and isinstance(sub.left, ast.Name)
                and isinstance(sub.right, ast.Name)
                and sub.left.id in components
                and sub.right.id in components
            ):
                continue
            a, b = sub.left.id, sub.right.id
            kinds = _differences(components[a], components[b])
            if not kinds:
                raise RecipeError(
                    f"'{a} - {b}': the two components are the same calculation"
                )
            if kinds not in ALLOWED_PAIR_DIFFERENCES:
                raise RecipeError(
                    f"'{a} - {b}': the two components differ in "
                    f"{', '.join(sorted(kinds))}; a difference must change one thing "
                    "(method, basis, frozen core, or the relativistic Hamiltonian "
                    "with its basis), so that it is the correction it is meant to be"
                )


def component_settings(recipe, comp):
    """The settings a component's program is given beyond its level --
    ``frozen_core`` (orbitals per element, "none" or "default"), ``relativistic``
    and ``nucleus`` -- and a label for the level that tells them apart, so two
    components of one level (all electrons and frozen core) are never mixed up."""
    settings = {}
    label = []
    fc = comp.get("frozen_core")
    if fc is not None:
        sets = recipe.data.get("frozen_core_definitions") or {}
        settings["frozen_core"] = sets.get(fc, fc)
        label.append(f"fc={frozen_core_key(settings['frozen_core'])}")
    rel = comp.get("relativistic", "none")
    if rel != "none":
        settings["relativistic"] = rel
        settings["nucleus"] = comp["nucleus"]
        label.append(f"{rel}({comp['nucleus']})")
    # A fixed order: the frozen core, then the Hamiltonian
    return settings, "".join(f"|{x}" for x in label)


def frozen_core_key(frozen_core):
    """A canonical key for a frozen core: "none", "default", or a short hash of
    the per-element counts, so the same counts give the same key whatever a recipe
    calls the set, and different counts never share one."""
    if frozen_core in FROZEN_CORE_SPECIAL:
        return frozen_core
    text = json.dumps(
        {str(k): int(v) for k, v in sorted(frozen_core.items())},
        sort_keys=True,
        separators=(",", ":"),
    )
    return hashlib.sha256(text.encode()).hexdigest()[:8]


def _check_level(what, level, allow_scale=False):
    if not isinstance(level, dict):
        raise RecipeError(f"The {what} must be a mapping")
    unknown = set(level) - COMPONENT_FIELDS - {"grid"}
    if unknown:
        raise RecipeError(f"Unknown fields in the {what}: {', '.join(sorted(unknown))}")
    if "level" in level:
        if "method" in level or "basis" in level:
            raise RecipeError(f"The {what} gives 'level' and 'method'/'basis'")
        if not isinstance(level["level"], str) or "@" not in level["level"]:
            raise RecipeError(
                f"The {what}: 'level' must be a model chemistry, e.g. "
                "'ORCA:DFT@PBE0/def2-TZVP'"
            )
    elif "method" not in level:
        raise RecipeError(f"The {what} needs 'level', or 'method' (and 'basis')")


def _degree(node, degrees, target):
    """The degree of an expression in the energies (energy^d), or a RecipeError
    when it has none: a sum of an energy and a number, a function given an
    argument of the wrong kind."""

    def fail(why):
        raise RecipeError(f"'combine' line for '{target}': {why}")

    def deg(n):
        if isinstance(n, ast.Constant):
            return 0
        if isinstance(n, ast.Name):
            return degrees[n.id]
        if isinstance(n, ast.UnaryOp):
            return deg(n.operand)
        if isinstance(n, ast.BinOp):
            left, right = deg(n.left), deg(n.right)
            if isinstance(n.op, (ast.Add, ast.Sub)):
                if left != right:
                    fail(
                        "adds or subtracts quantities of different kinds (an energy "
                        "and a number, or an energy and a product of energies)"
                    )
                return left
            if isinstance(n.op, ast.Mult):
                return left + right
            return left - right  # Div
        if isinstance(n, ast.Call):
            name = n.func.id
            wanted, result = DEGREES[name]
            given = [deg(a) for a in n.args]
            keywords = (
                list(inspect.signature(FUNCTIONS[name][0]).parameters)
                if name in FUNCTIONS
                else []
            )
            for kw in n.keywords:
                given.append(
                    (
                        keywords.index(kw.arg) if kw.arg in keywords else None,
                        deg(kw.value),
                    )
                )
            positional = [g for g in given if not isinstance(g, tuple)]
            for i, d in enumerate(positional):
                if d != wanted[i]:
                    fail(
                        f"argument {i + 1} of {name}() must be "
                        + ("an energy" if wanted[i] == 1 else "a number")
                    )
            for g in given:
                if isinstance(g, tuple):
                    i, d = g
                    if i is None or d != wanted[i]:
                        fail(f"a keyword argument of {name}() is of the wrong kind")
            return result
        fail(f"cannot find the degree of {type(n).__name__}")

    return deg(node)


def load_recipe(path, source="shipped"):
    """Read and check one recipe file."""
    path = Path(path)
    try:
        data = yaml.safe_load(path.read_text(encoding="utf-8"))
    except yaml.YAMLError as e:
        raise RecipeError(f"{path}: not valid YAML: {e}")
    try:
        return Recipe.from_dict(data, path=path, source=source)
    except RecipeError as e:
        raise RecipeError(f"{path}: {e}") from None


# -- The library of recipes ------------------------------------------------------


def shipped_directory():
    """The directory of the recipes released with composite_step."""
    import importlib.resources

    return Path(str(importlib.resources.files("composite_step") / "data" / "recipes"))


class RecipeLibrary:
    """The recipes available: the shipped set, the user directory and any job's.

    A local (user or job) recipe with the name of a shipped one is refused unless
    it carries another version, so a shipped name and version always mean the
    shipped recipe. Problems are collected in :attr:`problems`, not raised, so one
    bad local file does not hide the others.
    """

    def __init__(self, user_directory=None, extra_directories=(), shipped=True):
        self.recipes = {}  # (name, version) -> Recipe
        self.problems = []
        if shipped:
            self._add_directory(shipped_directory(), "shipped")
        if user_directory is not None:
            self._add_directory(Path(user_directory), "user")
        for directory in extra_directories:
            self._add_directory(Path(directory), "job")

    def _add_directory(self, directory, source):
        if not directory.is_dir():
            return
        for path in sorted(directory.glob("*.yaml")) + sorted(directory.glob("*.yml")):
            try:
                recipe = load_recipe(path, source=source)
            except RecipeError as e:
                self.problems.append(str(e))
                continue
            if source == "shipped" and recipe.uses_default_frozen_core:
                self.problems.append(
                    f"{path}: a shipped recipe must give its frozen cores explicitly, "
                    "not 'default'; ignored"
                )
                continue
            key = (recipe.name, recipe.version)
            if key in self.recipes:
                other = self.recipes[key]
                self.problems.append(
                    f"{path}: recipe {recipe.name} version {recipe.version} is already "
                    f"defined by {other.path} ({other.source}); ignored"
                )
                continue
            self.recipes[key] = recipe

    def names(self):
        return sorted({name for name, _ in self.recipes})

    def versions(self, name):
        return sorted(v for n, v in self.recipes if n == name)

    def get(self, name, version=None):
        """The recipe of that name: the given version, or the newest."""
        versions = self.versions(name)
        if not versions:
            raise KeyError(f"No recipe '{name}'")
        if version is None:
            version = versions[-1]
        try:
            return self.recipes[(name, version)]
        except KeyError:
            raise KeyError(
                f"No version {version} of the recipe '{name}' (have {versions})"
            ) from None

    def __iter__(self):
        return iter(self.recipes.values())

    def __len__(self):
        return len(self.recipes)


# -- Capabilities ----------------------------------------------------------------


def check_capabilities(recipe, available_levels=None, capabilities=None):
    """What a code lacks to run a recipe; an empty list if it can.

    Parameters
    ----------
    recipe : Recipe
    available_levels : collection of str, or callable, optional
        The model chemistries available (e.g. the keys of
        ``discover_model_chemistries()``), compared case-insensitively; or a
        function of a level spec returning whether it is available (e.g. one
        using ``match_model_chemistry``, which lets the basis be any).
    capabilities : collection of str, optional
        What the code offers for the recipe's ``requires`` list, as the same
        strings (``CCSD(T)``, ``basis:cc-pVTZ-F12``, ``F12/RI``, ...).
    """
    missing = []
    if available_levels is not None:
        if callable(available_levels):
            is_available = available_levels
        else:
            have = {str(x).lower() for x in available_levels}

            def is_available(spec):
                return spec.lower() in have

        levels = [("geometry", recipe.geometry)] + [
            (f"component '{k}'", v) for k, v in recipe.components.items()
        ]
        for what, level in levels:
            spec = level.get("level")
            if spec is not None and not is_available(spec):
                missing.append(f"{what}: the model chemistry {spec}")
    if capabilities is not None:
        have = {str(x).lower() for x in capabilities}
        for item in recipe.requires:
            if item.lower() not in have:
                missing.append(item)
    return missing
