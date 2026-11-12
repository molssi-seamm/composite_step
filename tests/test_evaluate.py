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

    def submit(self, configuration, key, options=None):
        self.submitted.append(key)
        self.options = getattr(self, "options", {})
        self.options[key] = options

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
    record = data["calculation"]
    assert record["program"] == "SEAMM composite"
    assert record["components"]["revdsd"]["record"] == {"fingerprint": "fp"}
    assert data["composite"]["components"]["revdsd"]["model chemistry"] == level
    assert data["composite"]["recipe"] == "revDSD-r2SCAN-3c-HO"
    # The prediction is the composite record's fingerprint, as stored
    assert CompositeStep.expected_fingerprints(object(), mc) == [record["fingerprint"]]
    assert CompositeStep.can_run_task(object(), mc) is True
    assert CompositeStep.make_evaluator(None, mc) is None


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


# -- the evaluator of a multi-component recipe ------------------------------------


def two_point():
    """A recipe with a cheap and an expensive component, no empirical terms."""
    return R.Recipe.from_dict(
        {
            "recipe": "two",
            "version": 1,
            "geometry": {"level": "A:DFT@GEOM/b"},
            "components": {
                "ccsdt": {"level": "A:QC@CCSD(T)/small"},
                "mp2_large": {"level": "A:QC@MP2/large"},
                "mp2_small": {"level": "A:QC@MP2/small"},
            },
            "combine": "E0 = ccsdt + (mp2_large - mp2_small)\n",
        }
    )


@pytest.fixture
def composite(monkeypatch):
    FakeEvaluator.calls = []
    made = []

    def factory(node, mc, properties, directory, name):
        evaluator = FakeEvaluator(node, mc, properties, directory, name)
        made.append(evaluator)
        return evaluator

    def make(**kwargs):
        mc = {"level": "SEAMM:CMP@two", "method": "two", "options": {}}
        return E.CompositeEvaluator(
            None,
            mc,
            two_point(),
            properties=("energy",),
            directory="/job/3",
            available=AVAILABLE,
            evaluator_factory=factory,
            **kwargs,
        )

    return make, made


def test_composite_evaluator(composite):
    make, made = composite
    with make() as evaluator:
        evaluator.submit("conf1", key="s1", options={"charge": 0})
        evaluator.submit("conf2", key="s2")
        results = {r.key: r for r in evaluator.results()}
    # Cheapest first: MP2 before CCSD(T); each component in its own directory
    levels = [c[0] for c in FakeEvaluator.calls]
    assert levels.index("A:QC@CCSD(T)/small") > levels.index("A:QC@MP2/large")
    assert {str(c[1]) for c in FakeEvaluator.calls} == {
        "/job/3/ccsdt",
        "/job/3/mp2_large",
        "/job/3/mp2_small",
    }
    # Per-structure options reach every component
    assert all(e.options["s1"] == {"charge": 0} for e in made)
    r1 = results["s1"]
    assert r1.ok and r1.energy == pytest.approx(-1000.0 + (-1010.0 + 995.0))
    assert r1.data["energy"] == r1.energy
    assert r1.data["composite"]["recipe"] == "two"
    assert r1.data["calculation"]["components"]["ccsdt"]["fingerprint"] == "fp:s1"
    assert results["s2"].ok
    # Nothing is computed twice
    assert list(evaluator.results()) == []


def test_failed_component_skips_the_rest(composite, monkeypatch):
    make, made = composite
    monkeypatch.setitem(ENERGIES, "A:QC@MP2/small", {"s1": -995.0})
    with make() as evaluator:
        evaluator.submit("conf1", key="s1")
        evaluator.submit("conf2", key="s2")
        results = {r.key: r for r in evaluator.results()}
    assert not results["s2"].ok
    assert "mp2_small" in results["s2"].reason
    # s2 failed its cheap component and never went to CCSD(T)
    ccsdt = [e for e in made if e.mc["level"] == "A:QC@CCSD(T)/small"][0]
    assert ccsdt.submitted == ["s1"]


def test_same_form_on_both_routes(composite, delegated):
    from composite_step import CompositeStep

    make, _ = composite
    with make() as evaluator:
        evaluator.submit("conf1", key="s1")
        (result,) = [r for r in evaluator.results() if r.key == "s1"]
    mc, _ = delegated
    data = CompositeStep.analyze_task(None, mc, object(), properties=("energy",))
    assert set(result.data) == set(data) >= {"energy", "composite", "calculation"}
    assert set(result.data["composite"]) == set(data["composite"])
    assert set(result.data["calculation"]) == set(data["calculation"])


def test_composite_record_fingerprint():
    recipe = two_point()
    mc = {"level": "SEAMM:CMP@two"}

    def record(fp):
        components = {
            name: {"record": {"fingerprint": fp + name}, "model_chemistry": name}
            for name in recipe.components
        }
        return E.composite_record(recipe, mc, components)["fingerprint"]

    assert record("a") == record("a")
    assert record("a") != record("b")


def test_composite_evaluator_refuses_gradients(composite):
    make, _ = composite
    with pytest.raises(ValueError, match="energies only"):
        E.CompositeEvaluator(
            None,
            {"level": "SEAMM:CMP@two"},
            two_point(),
            properties=("energy", "gradients"),
        )


def test_structure_dependent_empirical_terms_are_refused():
    with pytest.raises(R.RecipeError, match="depends on the structure"):
        E.CompositeEvaluator(
            None, {"level": "SEAMM:CMP@toy"}, g4mp2_like(), properties=("energy",)
        )


def test_cancel(composite):
    make, _ = composite
    evaluator = make()
    evaluator.submit("conf1", key="s1")
    evaluator.cancel()
    (result,) = list(evaluator.results())
    assert not result.ok and result.reason == "cancelled"


def test_make_evaluator_for_several_components(monkeypatch):
    from composite_step import CompositeStep

    monkeypatch.setattr(E, "_recipe_of", lambda mc, data_path=None: two_point())
    mc = {"level": "SEAMM:CMP@two", "method": "two", "options": {}}
    evaluator = CompositeStep.make_evaluator(None, mc, properties=("energy",))
    assert isinstance(evaluator, E.CompositeEvaluator)


def test_a_user_recipe_is_found_from_its_path(tmp_path):
    (tmp_path / "mine.yaml").write_text(
        "recipe: mine\nversion: 2\ngeometry: {level: 'A:DFT@GEOM/b'}\n"
        "components:\n  a: {level: 'A:QC@MP2/large'}\n  b: {level: 'A:QC@MP2/small'}\n"
        "combine: |\n  E0 = a - b\n"
    )
    mc = {
        "level": "SEAMM:CMP@mine",
        "method": "mine",
        "options": {
            "recipe": "mine",
            "recipe_version": 2,
            "recipe_source": "user",
            "recipe_path": str(tmp_path / "mine.yaml"),
        },
    }
    recipe = E._recipe_of(mc)
    assert recipe.name == "mine" and recipe.version == 2
    with pytest.raises(KeyError):
        E._recipe_of(dict(mc, options=dict(mc["options"], recipe_version=3)))


def test_looked_up_components_keep_their_fingerprints(composite):
    """A component from a store carries its fingerprint into the composite
    record, so the record is the same as when it is computed."""
    make, _ = composite
    with make() as evaluator:
        evaluator.submit("conf1", key="s1")
        (computed,) = list(evaluator.results())

    def lookup(source, mc, key, configuration):
        energy = ENERGIES[mc["level"]][key]
        return energy, "store", {"fingerprint": f"fp:{key}"}

    with make(component_lookup=lookup) as evaluator:
        evaluator.submit("conf1", key="s1")
        (stored,) = list(evaluator.results())
    assert stored.energy == computed.energy
    assert stored.data["composite"]["components"]["ccsdt"]["how"] == "store"
    assert (
        stored.data["calculation"]["fingerprint"]
        == computed.data["calculation"]["fingerprint"]
    )
