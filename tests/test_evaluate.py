"""Tests of running recipes: the components through (fake) evaluators, the
combination, and the one-energy task delegation."""

import pytest

from composite_step import evaluate as E
from composite_step import recipes as R

LEVELS = {
    "A:QC@CCSD(T)/small": "A",
    "A:QC@MP2/large": "A",
    "A:QC@MP2/small": "A",
    "A:DFT@GEOM/b": "A",
}
AVAILABLE = {
    level: {"level": level, "step": step, "options": {}}
    for level, step in LEVELS.items()
}

# Fake energies (kJ/mol) per level and structure key
ENERGIES = {
    "A:QC@CCSD(T)/small": {"s1": -1000.0, "s2": -2000.0},
    "A:QC@MP2/large": {"s1": -1010.0, "s2": -2020.0},
    "A:QC@MP2/small": {"s1": -995.0, "s2": -1990.0},
    "A:DFT@GEOM/b": {"s1": -1003.0, "s2": -2003.0},
}


class FakeResult:
    def __init__(self, key, energy, ok=True, reason=None):
        self.key = key
        self.energy = energy
        self.ok = ok
        self.reason = reason
        self.data = {"calculation": {"fingerprint": f"fp:{key}"}} if ok else {}


class FakeEvaluator:
    calls = []

    def __init__(self, node, mc, properties, directory, name):
        assert properties == ("energy",)
        self.mc = mc
        self.submitted = []
        FakeEvaluator.calls.append((mc["level"], directory, name))

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False

    def submit(self, configuration, key):
        self.submitted.append(key)

    def results(self):
        for key in self.submitted:
            energy = ENERGIES[self.mc["level"]].get(key)
            if energy is None:
                yield FakeResult(key, None, ok=False, reason="SCF did not converge")
            else:
                yield FakeResult(key, energy)


def g4mp2_like():
    return R.Recipe.from_dict(
        {
            "recipe": "toy",
            "version": 1,
            "geometry": {"level": "A:DFT@GEOM/b"},
            "components": {
                "ccsdt": {"level": "A:QC@CCSD(T)/small"},
                "mp2_large": {"level": "A:QC@MP2/large"},
                "mp2_small": {"level": "A:QC@MP2/small"},
            },
            "empirical": {"hlc": {"kind": "constant"}},
            "combine": "E0 = ccsdt + (mp2_large - mp2_small) + hlc\n",
        }
    )


def test_energy_sources_and_single_source():
    library = R.RecipeLibrary()
    quick = library.get("r2SCAN-3c-HO")
    default = library.get("revDSD-r2SCAN-3c-HO")
    assert E.energy_sources(quick) == {"geometry"}
    assert E.single_source(quick) == "geometry"
    assert E.energy_sources(default) == {"revdsd"}
    assert E.single_source(default) == "revdsd"
    assert E.component_levels(default) == {
        "revdsd": "ORCA:DFT@REVDSD-PBEP86-D4_2021/def2-TZVPP"
    }
    toy = g4mp2_like()
    assert E.energy_sources(toy) == {"ccsdt", "mp2_large", "mp2_small"}
    assert E.single_source(toy) is None


def test_evaluate_recipe(tmp_path):
    FakeEvaluator.calls = []
    toy = g4mp2_like()

    def lookup(source, mc, key, configuration):
        if source == "ccsdt" and key == "s1":
            return -1000.0, "store"
        return None

    out = E.evaluate_recipe(
        None,
        toy,
        {"s1": object(), "s2": object()},
        directory=tmp_path,
        lookup=lookup,
        constants={"hlc": -0.01},
        available=AVAILABLE,
        evaluator_factory=FakeEvaluator,
    )
    assert list(out) == ["s1", "s2"]
    hlc = -0.01 * E.HARTREE_TO_KJ
    assert out["s1"].ok and out["s1"].E0 == pytest.approx(-1000 - 15 + hlc)
    assert out["s2"].E0 == pytest.approx(-2000 - 30 + hlc)
    assert out["s1"].components["ccsdt"]["how"] == "store"
    assert out["s2"].components["ccsdt"]["how"] == "computed"
    assert out["s2"].components["ccsdt"]["record"] == {"fingerprint": "fp:s2"}
    # One evaluator per component, each in its own directory
    assert sorted(level for level, _, _ in FakeEvaluator.calls) == sorted(
        ["A:QC@CCSD(T)/small", "A:QC@MP2/large", "A:QC@MP2/small"]
    )
    assert {d.name for _, d, _ in FakeEvaluator.calls} == {
        "ccsdt",
        "mp2_large",
        "mp2_small",
    }


def test_a_failed_component_fails_only_its_structure():
    ENERGIES["A:QC@MP2/large"].pop("s2")
    try:
        out = E.evaluate_recipe(
            None,
            g4mp2_like(),
            {"s1": object(), "s2": object()},
            constants={"hlc": 0.0},
            available=AVAILABLE,
            evaluator_factory=FakeEvaluator,
        )
    finally:
        ENERGIES["A:QC@MP2/large"]["s2"] = -2020.0
    assert out["s1"].ok
    assert not out["s2"].ok
    assert "mp2_large" in out["s2"].reason and "SCF" in out["s2"].reason


def test_missing_empirical_values():
    with pytest.raises(R.RecipeError, match="empirical values hlc"):
        E.evaluate_recipe(
            None,
            g4mp2_like(),
            {"s1": object()},
            available=AVAILABLE,
            evaluator_factory=FakeEvaluator,
        )


def test_single_component_recipe_is_exact():
    """A one-energy recipe gives the component's energy itself (no unit round
    trip): the default tier must reproduce the step's energy level exactly."""
    default = R.RecipeLibrary().get("revDSD-r2SCAN-3c-HO")
    level = "ORCA:DFT@REVDSD-PBEP86-D4_2021/def2-TZVPP"
    available = {level: {"level": level, "step": "ORCA", "options": {}}}
    ENERGIES[level] = {"s1": -405_398.55555321}
    out = E.evaluate_recipe(
        None,
        default,
        {"s1": object()},
        available=available,
        evaluator_factory=FakeEvaluator,
    )
    assert out["s1"].E0 == -405_398.55555321


class FakeProvider:
    @classmethod
    def get_task(cls, configuration, mc, **kwargs):
        return ("task", mc["level"], kwargs.get("key"))

    @classmethod
    def analyze_task(cls, result, mc, configuration, **kwargs):
        return {"energy": -123.25, "calculation": {"fingerprint": "fp"}}

    @classmethod
    def expected_fingerprints(cls, configuration, mc, **kwargs):
        return ["fp"]


@pytest.fixture
def delegated(monkeypatch):
    level = "ORCA:DFT@REVDSD-PBEP86-D4_2021/def2-TZVPP"
    monkeypatch.setattr(
        E, "resolve_level", lambda spec, available=None: {"level": spec, "step": "X"}
    )
    monkeypatch.setitem(E._PROVIDERS, "X", FakeProvider)
    mc = {
        "level": "SEAMM:CMP@revDSD-r2SCAN-3c-HO",
        "method": "revDSD-r2SCAN-3c-HO",
        "options": {"recipe": "revDSD-r2SCAN-3c-HO", "recipe_version": 1},
    }
    return mc, level


def test_delegated_task_route(delegated):
    from composite_step import CompositeStep

    mc, level = delegated
    task = CompositeStep.get_task(object(), mc, key="k1")
    assert task == ("task", level, "k1")
    data = CompositeStep.analyze_task(None, mc, object(), properties=("energy",))
    assert data["energy"] == -123.25
    assert data["calculation"] == {"fingerprint": "fp"}
    assert data["composite"]["component model chemistry"] == level
    assert data["composite"]["recipe"] == "revDSD-r2SCAN-3c-HO"
    assert CompositeStep.expected_fingerprints(object(), mc) == ["fp"]
    assert CompositeStep.can_run_task(object(), mc) is True


def test_no_gradients_yet(delegated):
    from seamm_exec import AnalysisError

    from composite_step import CompositeStep

    mc, _ = delegated
    with pytest.raises(AnalysisError, match="no gradients"):
        CompositeStep.analyze_task(
            None, mc, object(), properties=("energy", "gradients")
        )


def test_several_components_refuse_the_one_task_route(tmp_path, monkeypatch):
    from composite_step import CompositeStep

    (tmp_path / "recipes").mkdir()
    (tmp_path / "recipes" / "toy.yaml").write_text(
        "recipe: toy\nversion: 1\ngeometry: {level: 'A:DFT@GEOM/b'}\n"
        "components:\n  a: {level: 'A:QC@MP2/large'}\n  b: {level: 'A:QC@MP2/small'}\n"
        "combine: |\n  E0 = a - b\n"
    )
    monkeypatch.setattr(
        E,
        "_recipe_of",
        lambda mc, data_path=None: R.load_recipe(tmp_path / "recipes" / "toy.yaml"),
    )
    mc = {"level": "SEAMM:CMP@toy", "method": "toy", "options": {}}
    assert CompositeStep.can_run_task(object(), mc) is False
    with pytest.raises(R.RecipeError, match="combines several calculations"):
        CompositeStep.get_task(object(), mc, key="k")
