"""End to end: a consumer's ``seamm_exec.Evaluator(...)`` for a multi-component
recipe becomes a CompositeEvaluator, whose components run as real (shell) tasks
of a fake program through seamm_exec's TaskSet, and are combined."""

import logging
import types

import pytest

seamm_exec = pytest.importorskip("seamm_exec")
if "__new__" not in vars(seamm_exec.Evaluator):
    pytest.skip("seamm_exec has no make_evaluator hook", allow_module_level=True)

from composite_step import CompositeStep  # noqa: E402
from composite_step import evaluate as E  # noqa: E402
from composite_step import recipes as R  # noqa: E402

#: The fake program's energy: the sum of the coordinates times a factor per method
FACTORS = {"HF/big": 1.0, "HF/small": 2.0, "CC/small": 3.0}


class FakeProgram:
    calls = []

    @classmethod
    def get_task(cls, configuration, model_chemistry, *, key, properties, options):
        data = seamm_exec.structure_data(configuration)
        factor = FACTORS[model_chemistry["method"]]
        energy = factor * float(data["coordinates"].sum()) + options.get("shift", 0.0)
        cls.calls.append((model_chemistry["method"], key))
        return seamm_exec.Task(
            key=key,
            program="fake",
            cmd=["echo", str(energy), ">", "energy.txt"],
            shell=True,
            return_files=["energy.txt"],
            config={"installation": "local"},
        )

    @classmethod
    def analyze_task(
        cls, result, model_chemistry, configuration, *, properties, options
    ):
        return {
            "energy": float(result.read("energy.txt")),
            "calculation": {"fingerprint": f"fp-{model_chemistry['method']}"},
        }


RECIPE = R.Recipe.from_dict(
    {
        "recipe": "toy-cbs",
        "version": 1,
        "geometry": {"level": "FAKE:QC@HF/small"},
        "components": {
            "hf_big": {"level": "FAKE:QC@HF/big"},
            "hf_small": {"level": "FAKE:QC@HF/small"},
            "cc_small": {"level": "FAKE:QC@CC/small"},
        },
        "combine": "E0 = cc_small + (hf_big - hf_small)\n",
    }
)
MC = {
    "level": "SEAMM:CMP@toy-cbs",
    "method": "toy-cbs",
    "step": "Composite",
    "options": {"recipe": "toy-cbs", "mdi_capable": False},
}


@pytest.fixture
def node(tmp_path, monkeypatch):
    providers = {"Composite": CompositeStep, "fake": FakeProgram}
    monkeypatch.setattr(E, "_recipe_of", lambda mc, data_path=None: RECIPE)

    def resolve(level, available=None):
        method = level.split("@")[1]
        return {"level": level, "method": method, "step": "fake", "options": {}}

    monkeypatch.setattr(E, "resolve_level", resolve)
    FakeProgram.calls = []
    return types.SimpleNamespace(
        directory=str(tmp_path / "step"),
        global_options={"root": str(tmp_path)},
        logger=logging.getLogger("test"),
        variable_exists=lambda name: False,
        flowchart=types.SimpleNamespace(
            executor=seamm_exec.Local(),
            plugin_manager=types.SimpleNamespace(get=lambda name: providers[name]),
            root_directory=str(tmp_path),
        ),
    )


def test_consumer_gets_composite_energies(node, tmp_path):
    geometries = {
        "w": seamm_exec.Geometry([8, 1, 1], [[0, 0, 0], [0.96, 0, 0], [0, 0.96, 0]]),
        "m": seamm_exec.Geometry([6], [[0.5, 0.25, 0.0]]),
    }
    with seamm_exec.Evaluator(
        node, MC, properties=("energy",), path="batch", name="cbs"
    ) as evaluator:
        assert isinstance(evaluator, E.CompositeEvaluator)
        for key, geometry in geometries.items():
            evaluator.submit(geometry, key=key, options={"shift": 0.5})
        results = {r.key: r for r in evaluator.results()}

    for key, geometry in geometries.items():
        s = float(seamm_exec.structure_data(geometry)["coordinates"].sum())
        # Each component: factor * s + shift; E0 = cc + (hf_big - hf_small)
        expected = (3.0 * s + 0.5) + ((1.0 * s + 0.5) - (2.0 * s + 0.5))
        assert results[key].ok, results[key].reason
        assert results[key].energy == pytest.approx(expected, abs=1e-6)
        record = results[key].data["calculation"]
        assert record["components"]["cc_small"]["fingerprint"] == "fp-CC/small"
    # Every component ran each structure once, as its own task set
    assert sorted(FakeProgram.calls) == sorted(
        (m, k) for m in FACTORS for k in geometries
    )
    for source in RECIPE.components:
        assert (tmp_path / "step" / source).is_dir()


def test_one_energy_recipes_keep_the_ordinary_evaluator(node, monkeypatch):
    single = R.Recipe.from_dict(
        {
            "recipe": "one",
            "version": 1,
            "geometry": {"level": "FAKE:QC@HF/small"},
            "components": {"cc": {"level": "FAKE:QC@CC/small"}},
            "combine": "E0 = cc\n",
        }
    )
    monkeypatch.setattr(E, "_recipe_of", lambda mc, data_path=None: single)
    evaluator = seamm_exec.Evaluator(node, MC, properties=("energy",), path="batch")
    assert isinstance(evaluator, seamm_exec.Evaluator)
