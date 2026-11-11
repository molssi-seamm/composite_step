"""Tests of the recipe library: format, combination, capabilities, thermal terms."""

import math
import textwrap

import pytest

from composite_step import recipes as R


def write(directory, name, text):
    path = directory / name
    path.write_text(textwrap.dedent(text))
    return path


MINIMAL = {
    "recipe": "test",
    "version": 1,
    "geometry": {"level": "ORCA:DFT@PBE0/def2-SVP"},
    "combine": "E0 = geometry\n",
}


def recipe(**changes):
    data = dict(MINIMAL)
    data.update(changes)
    return R.Recipe.from_dict(data)


# -- the shipped recipes ------------------------------------------------------------


def test_shipped_recipes_load():
    library = R.RecipeLibrary()
    assert library.problems == []
    assert {"r2SCAN-3c-HO", "revDSD-r2SCAN-3c-HO"} <= set(library.names())
    default = library.get("revDSD-r2SCAN-3c-HO")
    assert default.model_chemistry == "SEAMM:CMP@revDSD-r2SCAN-3c-HO"
    assert default.geometry == {"level": "ORCA:CMP@R2SCAN-3C"}
    assert default.components["revdsd"]["level"] == (
        "ORCA:DFT@REVDSD-PBEP86-D4_2021/def2-TZVPP"
    )


def test_shipped_zpe_factor_is_the_tables():
    """The tiers pin the ZPE factor; it must be the one the thermochemistry step
    takes from seamm_thermochemistry's table today."""
    thermo = pytest.importorskip("seamm_thermochemistry")
    table = thermo.zpe_scale_factor("ORCA:CMP@R2SCAN-3C")
    for name in ("r2SCAN-3c-HO", "revDSD-r2SCAN-3c-HO"):
        s, source = R.RecipeLibrary().get(name).zpe_scale
        assert s == table.factor
        assert table.source in source


def _step_arithmetic(E_high, low, s):
    """The thermochemistry step's own arithmetic (cbh_approach.py, 2026.11.11.4):
    combine the energy level with the geometry level, then shift by the scaled
    ZPE. Kept literally, to check the recipe reproduces it to the last bit."""
    combined = {"E": E_high}
    combined["H0"] = E_high + low["ZPE"]
    combined["H298"] = E_high + (low["H298"] - low["E"])
    if s != 1.0:
        shift = (s - 1.0) * low["ZPE"]
        combined["H0"] += shift
        combined["H298"] += shift
    return combined


@pytest.mark.parametrize("name", ["r2SCAN-3c-HO", "revDSD-r2SCAN-3c-HO"])
def test_tiers_reproduce_the_step_bitwise(name):
    r = R.RecipeLibrary().get(name)
    s, _ = r.zpe_scale
    # kJ/mol, ethanol-like magnitudes
    low = {"E": -405_312.123456789, "ZPE": 205.4321987, "H298": -405_094.98765}
    E_revdsd = -405_398.55555321
    components = {"geometry": low["E"], "revdsd": E_revdsd}
    E0 = r.combine(components)["E0"]
    got = r.thermal(E0, {"E": low["E"], "ZPE": low["ZPE"], "H": low["H298"]})
    want = _step_arithmetic(
        E_revdsd if name == "revDSD-r2SCAN-3c-HO" else low["E"], low, s
    )
    assert got["H0"] == want["H0"]
    assert got["H"] == want["H298"]
    assert got["ZPE scaled"] == s * low["ZPE"]


# -- validation ---------------------------------------------------------------------


@pytest.mark.parametrize(
    "changes, message",
    [
        ({"recipe": "G4/MP2"}, "may not contain"),
        ({"recipe": "a:b"}, "may not contain"),
        ({"version": 0}, "positive integer"),
        ({"version": "1"}, "positive integer"),
        ({"bogus": 1}, "Unknown recipe fields"),
        ({"combine": "x = geometry\n"}, "must define E0"),
        ({"combine": "E0 = nothing\n"}, "unknown names: nothing"),
        ({"combine": "geometry = 1\nE0 = geometry\n"}, "redefines"),
        ({"geometry": {"basis": "def2-SVP"}}, "needs 'level', or 'method'"),
        ({"geometry": {"level": "PBE0"}}, "must be a model chemistry"),
        ({"frequencies": {"scale_zpe": 0.98}}, "needs its 'scale_source'"),
        ({"components": {"geometry": {"level": "A:B@C"}}}, "reserved"),
        ({"components": {"x y": {"level": "A:B@C"}}}, "identifier"),
    ],
)
def test_invalid_recipes(changes, message):
    with pytest.raises(R.RecipeError, match=message):
        recipe(**changes)


def test_missing_fields():
    with pytest.raises(R.RecipeError, match="Missing recipe fields: combine"):
        R.Recipe.from_dict({k: v for k, v in MINIMAL.items() if k != "combine"})


@pytest.mark.parametrize(
    "line",
    [
        "E0 = __import__('os')",
        "E0 = geometry.real",
        "E0 = [geometry][0]",
        "E0 = (lambda: 1)()",
        "E0 = eval('1')",
        "E0 = geometry ** 2",
        "E0 = 'text'",
        "E0 = True",
        "E0 = scale(geometry)",
        "E0 = hf_part(geometry + 1)",
        "import os\nE0 = geometry",
        "E0, x = geometry, 1",
    ],
)
def test_combine_whitelist(line):
    with pytest.raises(R.RecipeError):
        recipe(combine=line + "\n")


# -- the combination ----------------------------------------------------------------


def test_cbs_power_is_exact_for_its_form():
    E, A = -76.3, 0.25
    e3, e4 = E + A * 3**-3, E + A * 4**-3
    assert R.cbs_power(e3, e4, 3, 4) == pytest.approx(E, abs=1e-12)
    e3, e4 = E + A * 3**-5, E + A * 4**-5
    assert R.cbs_power(e3, e4, 3, 4, alpha=5) == pytest.approx(E, abs=1e-12)


@pytest.mark.parametrize("alpha", [1.63, 9.0])
def test_hf_extrapolate_exp_is_exact_for_its_form(alpha):
    E, B, n = -76.06, 0.8, 4
    e = [E + B * (k + 1) * math.exp(-alpha * math.sqrt(k)) for k in (n, n + 1)]
    assert R.hf_extrapolate_exp(e[0], e[1], n, alpha) == pytest.approx(E, abs=1e-12)


def test_combination_with_parts_and_empirical():
    r = recipe(
        components={
            "small": {"method": "CCSD(T)", "basis": "6-31G(d)"},
            "large": {"method": "MP2", "basis": "G3MP2LargeXP"},
            "small_mp2": {"method": "MP2", "basis": "6-31G(d)"},
        },
        empirical={"hlc": {"kind": "constant"}},
        combine=(
            "delta = large - small_mp2\n"
            "corr = corr_part(large)\n"
            "E0 = small + delta + hlc\n"
        ),
    )
    results = r.combine(
        {
            "small": R.ComponentEnergy(-76.20, hf=-76.00),
            "large": R.ComponentEnergy(-76.35, hf=-76.05),
            "small_mp2": -76.19,
            "geometry": -76.4,
        },
        constants={"hlc": -0.01},
    )
    assert results["delta"] == pytest.approx(-0.16)
    assert results["corr"] == pytest.approx(-0.30)
    assert results["E0"] == pytest.approx(-76.20 - 0.16 - 0.01)


def test_part_of_a_component_without_hf_part():
    r = recipe(
        components={"c": {"method": "MP2", "basis": "x"}},
        combine="E0 = hf_part(c)\n",
    )
    with pytest.raises(R.RecipeError, match="has no HF part"):
        r.combine({"c": -1.0, "geometry": -1.0})


def test_thermal_without_scaling_or_H():
    r = recipe()
    t = r.thermal(-100.0, {"E": -99.0, "ZPE": 0.05})
    assert t["H0"] == -100.0 + 0.05 and "H" not in t


# -- the library --------------------------------------------------------------------

LOCAL = """\
recipe: {name}
version: {version}
geometry: {{level: "ORCA:DFT@PBE0/def2-SVP"}}
combine: |
  E0 = geometry
"""


def test_local_recipe_may_not_shadow_a_shipped_one(tmp_path):
    write(tmp_path, "a.yaml", LOCAL.format(name="r2SCAN-3c-HO", version=1))
    write(tmp_path, "b.yaml", LOCAL.format(name="r2SCAN-3c-HO", version=7))
    write(tmp_path, "c.yaml", LOCAL.format(name="mine", version=1))
    write(tmp_path, "d.yaml", "recipe: [broken\n")
    library = R.RecipeLibrary(user_directory=tmp_path)
    assert library.get("r2SCAN-3c-HO", 1).source == "shipped"
    assert library.get("r2SCAN-3c-HO").version == 7  # the newest
    assert library.get("mine").source == "user"
    assert len(library.problems) == 2
    assert any("already defined" in p for p in library.problems)
    assert any("not valid YAML" in p for p in library.problems)
    with pytest.raises(KeyError, match="No version 3"):
        library.get("mine", 3)
    with pytest.raises(KeyError, match="No recipe"):
        library.get("nothing")


# -- capabilities and fingerprints --------------------------------------------------


def test_capabilities():
    r = R.RecipeLibrary().get("revDSD-r2SCAN-3c-HO")
    assert r.requires == []
    have = ["ORCA:CMP@R2SCAN-3C", "orca:dft@revdsd-pbep86-d4_2021/def2-tzvpp"]
    assert R.check_capabilities(r, have) == []
    missing = R.check_capabilities(r, have[:1])
    assert missing == [
        "component 'revdsd': the model chemistry "
        "ORCA:DFT@REVDSD-PBEP86-D4_2021/def2-TZVPP"
    ]
    r2 = recipe(requires=["MP4(SDTQ)", "CCSD(T)"])
    assert R.check_capabilities(r2, capabilities=["ccsd(t)"]) == ["MP4(SDTQ)"]


def test_fingerprint():
    a = recipe(description="one")
    b = recipe(description="two", reference="somewhere")
    c = recipe(version=2)
    assert a.fingerprint() == b.fingerprint()
    assert a.fingerprint() != c.fingerprint()
    assert a.fingerprint({"x": "sha256:1"}) != a.fingerprint()
    assert a.fingerprint().startswith("sha256:")


def test_capabilities_with_a_matcher():
    r = R.RecipeLibrary().get("revDSD-r2SCAN-3c-HO")
    assert R.check_capabilities(r, lambda spec: spec.startswith("ORCA:")) == []
    assert len(R.check_capabilities(r, lambda spec: "CMP" in spec)) == 1


def test_model_chemistry_options(tmp_path):
    from composite_step import CompositeStep

    (tmp_path / "recipes").mkdir()
    write(tmp_path / "recipes", "mine.yaml", LOCAL.format(name="mine", version=2))
    options = CompositeStep.get_model_chemistry_options(data_path=[tmp_path])
    assert {
        "SEAMM:CMP@r2SCAN-3c-HO",
        "SEAMM:CMP@revDSD-r2SCAN-3c-HO",
        "SEAMM:CMP@mine",
    } <= set(options)
    default = options["SEAMM:CMP@revDSD-r2SCAN-3c-HO"]
    assert default["type"] == "CMP" and default["basis"] is None
    assert default["geometry_level"] == "ORCA:CMP@R2SCAN-3C"
    assert default["components"] == {
        "revdsd": "ORCA:DFT@REVDSD-PBEP86-D4_2021/def2-TZVPP"
    }
    assert default["zpe_scale"] == 0.9899
    assert options["SEAMM:CMP@mine"]["recipe_source"] == "job"
    assert CompositeStep.get_model_chemistry_options(periodic_only=True) == {}
    assert CompositeStep.get_model_chemistry_options(mdi_only=True) == {}


@pytest.mark.parametrize(
    "combine, message",
    [
        ("E0 = geometry * geometry\n", "scales as energy\\^2"),
        ("E0 = geometry + 1\n", "different kinds"),
        ("x = geometry / geometry\nE0 = x\n", "scales as energy\\^0"),
        (
            "E0 = scale(geometry, geometry)\n",
            "argument 2 of scale\\(\\) must be a number",
        ),
        ("E0 = scale(0.5, geometry)\n", "argument 1 of scale\\(\\) must be an energy"),
    ],
)
def test_nonlinear_combinations_are_refused(combine, message):
    """Every result must scale like an energy, so a combination gives the same
    answer in any unit (kJ/mol, E_h)."""
    with pytest.raises(R.RecipeError, match=message):
        recipe(combine=combine)


def test_unit_invariance_of_allowed_forms():
    """Ratios of energy differences are numbers (degree 0), so forms that fit an
    exponent from them, like Feller's three-point extrapolation, are allowed and
    give the same result in E_h and in kJ/mol."""
    r = recipe(
        components={
            "dz": {"method": "HF", "basis": "cc-pVDZ"},
            "tz": {"method": "HF", "basis": "cc-pVTZ"},
            "qz": {"method": "HF", "basis": "cc-pVQZ"},
        },
        combine=("E0 = (tz * tz - dz * qz) / (2 * tz - dz - qz)\n"),
    )
    e = {"dz": -76.0267, "tz": -76.0571, "qz": -76.0653, "geometry": -76.1}
    k = 2625.4996394798254
    in_hartree = r.combine(e)["E0"]
    in_kj = r.combine({n: v * k for n, v in e.items()})["E0"]
    assert in_kj == pytest.approx(in_hartree * k, rel=1e-12)
    # and it is Feller's limit: E(n) = E + B exp(-C n) exactly
    import math

    E, B, C = -76.0675, 0.9, 1.4
    exact = {f: E + B * math.exp(-C * n) for f, n in (("dz", 2), ("tz", 3), ("qz", 4))}
    exact["geometry"] = 0.0
    assert r.combine(exact)["E0"] == pytest.approx(E, abs=1e-10)
