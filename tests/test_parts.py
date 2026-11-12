"""Composite phase 2 §2: the parts of the components' energies in recipes --
part(), hf_part and corr_part (correlation only, CABS singles apart), the derived
parts, the CABS-singles and reference checks, and parts through evaluate_recipe
(stored rows without parts are recomputed, never derived)."""

import pytest

from composite_step import evaluate as E
from composite_step import recipes as R


def _recipe(combine, components, **extra):
    return R.Recipe.from_dict(
        {
            "recipe": "p",
            "version": 1,
            "geometry": {"level": "A:DFT@GEOM/b"},
            "components": components,
            "combine": combine,
            **extra,
        }
    )


F12 = R.ComponentEnergy(
    total=-76.368957009,
    parts={
        "HF": -76.065298885,
        "CABS singles": -0.001072393,
        "CCSD correlation": -0.277268885,
        "F12 correction": -0.016427302,
        "(T)": -0.008889543,
    },
    factors={"(T) scale f": 1.059246828},
    reference="RHF",
)


def test_component_parts_and_derived_parts():
    assert F12.hf == F12.part("HF")
    # Correlation alone: E - HF - CABS singles
    assert F12.correlation == pytest.approx(-0.277268885 - 0.016427302 - 0.008889543)
    assert F12.part("CCSD-F12 correlation") == pytest.approx(-0.293696187)
    assert F12.part("HF+CABS") == pytest.approx(-76.066371278)
    assert F12.part("(T*)") == pytest.approx(1.059246828 * -0.008889543)
    with pytest.raises(R.RecipeError, match="no part 'DLPNO-\\(T0\\)'"):
        F12.part("DLPNO-(T0)")


def test_part_syntax_and_degree():
    r = _recipe(
        'E0 = hf_part(a) + corr_part(a) - part(a, "(T)") + part(a, "(T*)")\n',
        {"a": {"level": "A:QC@CCSD(T)/x"}},
    )
    assert r.needs_parts == {"a"}
    for bad in ("E0 = part(a)\n", "E0 = part(a, 3)\n", 'E0 = part("HF", a)\n'):
        with pytest.raises(R.RecipeError, match="part\\(\\) takes"):
            _recipe(bad, {"a": {"level": "A:QC@CCSD(T)/x"}})


def test_f12_split_must_account_for_cabs_singles():
    f12 = {"a": {"level": "A:QC@CCSD(T)-F12D_RI/cc-pVTZ-F12"}}
    with pytest.raises(R.RecipeError, match="leave out its CABS singles"):
        _recipe("E0 = hf_part(a) + 1.1 * corr_part(a)\n", f12)
    _recipe('E0 = hf_part(a) + part(a, "CABS singles") + 1.1 * corr_part(a)\n', f12)
    omitted = {"a": dict(f12["a"], omit=["CABS singles"])}
    _recipe("E0 = hf_part(a) + 1.1 * corr_part(a)\n", omitted)
    # Not F12: nothing to account for
    _recipe(
        "E0 = hf_part(a) + 1.1 * corr_part(a)\n", {"a": {"level": "A:QC@CCSD(T)/x"}}
    )


def test_mixed_references_are_refused_unless_allowed():
    comps = {"a": {"level": "A:QC@CCSD(T)/x"}, "b": {"level": "A:QC@CCSD(T)/y"}}
    energies = {
        "a": R.ComponentEnergy(-2.0, parts={"HF": -1.5}, reference="UHF"),
        "b": R.ComponentEnergy(-2.1, parts={"HF": -1.6}, reference="QRO"),
    }
    combine = "E0 = hf_part(a) + corr_part(b)\n"
    assert _recipe(combine, comps).combine(energies)["E0"] == pytest.approx(-2.0)
    combine = "E0 = 0.5 * (hf_part(a) + hf_part(b)) + corr_part(b)\n"
    with pytest.raises(R.RecipeError, match="different references"):
        _recipe(combine, comps).combine(energies)
    allowed = _recipe(combine, comps, allow_mixed_reference=True)
    assert allowed.combine(energies)["E0"] == pytest.approx(-1.55 - 0.5)


# -- parts through evaluate_recipe ---------------------------------------------------

LEVELS = ("A:QC@CCSD(T)/tz", "A:QC@CCSD(T)/qz")
#: E_h: (HF, CCSD correlation, (T)) per level and structure
VALUES = {
    "A:QC@CCSD(T)/tz": {"s1": (-76.05, -0.26, -0.008)},
    "A:QC@CCSD(T)/qz": {"s1": (-76.06, -0.28, -0.009)},
}
H = E.HARTREE_TO_KJ


class PartsEvaluator:
    calls = []
    without_parts = set()

    def __init__(self, node, mc, properties, directory, name):
        self.mc = mc
        self.keys = []

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False

    def submit(self, configuration, key, options=None):
        self.keys.append(key)

    def results(self):
        from seamm_exec.evaluator import EvaluatorResult

        level = self.mc["level"]
        for key in self.keys:
            PartsEvaluator.calls.append((level, key))
            hf, ccsd, t = VALUES[level][key]
            data = {"calculation": {"fingerprint": f"fp:{level}"}}
            if level not in PartsEvaluator.without_parts:
                data["energy parts"] = {
                    "units": "E_h",
                    "parts": {"HF": hf, "CCSD correlation": ccsd, "(T)": t},
                    "reference": "RHF",
                }
            yield EvaluatorResult(
                key=key, ok=True, energy=(hf + ccsd + t) * H, data=data
            )


CBS = _recipe(
    "E0 = hf_part(qz) + cbs_power(corr_part(tz), corr_part(qz), 3, 4)\n",
    {"tz": {"level": LEVELS[0]}, "qz": {"level": LEVELS[1]}},
)
AVAILABLE = {lv: {"level": lv, "step": "A", "options": {}} for lv in LEVELS}


@pytest.fixture
def parts():
    PartsEvaluator.calls = []
    PartsEvaluator.without_parts = set()
    return PartsEvaluator


def _expected():
    hf_tz, c_tz, t_tz = VALUES[LEVELS[0]]["s1"]
    hf_qz, c_qz, t_qz = VALUES[LEVELS[1]]["s1"]
    corr = R.cbs_power(c_tz + t_tz, c_qz + t_qz, 3, 4)
    return (hf_qz + corr) * H


def test_parts_in_evaluate_recipe(parts):
    out = E.evaluate_recipe(
        None,
        CBS,
        {"s1": "conf"},
        available=AVAILABLE,
        evaluator_factory=parts,
    )
    assert out["s1"].ok, out["s1"].reason
    assert out["s1"].E0 == pytest.approx(_expected(), rel=1e-12)


def test_a_component_without_parts_fails(parts):
    parts.without_parts = {LEVELS[0]}
    out = E.evaluate_recipe(
        None, CBS, {"s1": "conf"}, available=AVAILABLE, evaluator_factory=parts
    )
    assert not out["s1"].ok
    assert "gives no energy parts" in out["s1"].reason


def test_a_stored_row_without_parts_is_recomputed(parts):
    """The store has the tz energy but not its parts: tz is computed again;
    the qz row with parts is used."""

    def lookup(source, mc, key, configuration):
        hf, c, t = VALUES[mc["level"]][key]
        record = {"fingerprint": f"fp:{mc['level']}"}
        if source == "qz":
            record["energy parts"] = {
                "units": "E_h",
                "parts": {"HF": hf, "CCSD correlation": c, "(T)": t},
                "reference": "RHF",
            }
        return (hf + c + t) * H, "store", record

    out = E.evaluate_recipe(
        None,
        CBS,
        {"s1": "conf"},
        lookup=lookup,
        available=AVAILABLE,
        evaluator_factory=parts,
    )
    assert out["s1"].ok
    assert parts.calls == [(LEVELS[0], "s1")]
    components = out["s1"].components
    assert components["tz"]["how"] == "computed (store has no energy parts)"
    assert components["qz"]["how"] == "store"
    assert out["s1"].E0 == pytest.approx(_expected(), rel=1e-12)
