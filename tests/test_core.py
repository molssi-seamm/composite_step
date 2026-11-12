"""Composite phase 2 §3: frozen-core sets, relativistic Hamiltonians, the check
that a difference of two components changes one thing only, and the settings
each component runs with."""

import pytest

from composite_step import evaluate as E
from composite_step import recipes as R

SETS = {
    "valence": {"H": 0, "C": 1, "N": 1, "O": 1, "F": 1, "S": 5, "Cl": 5},
    "inner": {"H": 0, "C": 0, "N": 0, "O": 0, "F": 0, "S": 1, "Cl": 1},
}
F12 = "ORCA:QC@CCSD(T)-F12D_RI/cc-pVTZ-F12"
CV = "ORCA:QC@CCSD(T)/cc-pwCVTZ"


def tier_a(**changes):
    components = {
        "f12": {"level": F12, "frozen_core": "valence"},
        "cv_ae": {"level": CV, "frozen_core": "inner"},
        "cv_fc": {"level": CV, "frozen_core": "valence"},
        "rel": {
            "level": "ORCA:QC@CCSD(T)/cc-pVTZ-DK",
            "frozen_core": "valence",
            "relativistic": "DKH2",
            "nucleus": "point",
        },
        "nonrel": {"level": "ORCA:QC@CCSD(T)/cc-pVTZ", "frozen_core": "valence"},
    }
    for name, comp in changes.items():
        components[name] = comp
    return {
        "recipe": "tier-a-core",
        "version": 1,
        "geometry": {"level": "ORCA:CMP@R2SCAN-3C"},
        "frozen_core_definitions": SETS,
        "components": components,
        "combine": "E0 = f12 + (cv_ae - cv_fc) + (rel - nonrel)\n",
    }


def test_the_tier_a_core_recipe_is_valid():
    recipe = R.Recipe.from_dict(tier_a())
    assert not recipe.uses_default_frozen_core


@pytest.mark.parametrize(
    "change, message",
    [
        ({"cv_fc": {"level": CV, "frozen_core": "core"}}, "no frozen-core set 'core'"),
        (
            {"rel": {"level": "ORCA:QC@CCSD(T)/cc-pVTZ-DK", "relativistic": "DKH2"}},
            "nucleus model",
        ),
        (
            {"nonrel": {"level": "ORCA:QC@CCSD(T)/cc-pVTZ", "nucleus": "point"}},
            "without a relativistic",
        ),
        (
            {"rel": {"level": "x:QC@CCSD(T)/y", "relativistic": "DKH3"}},
            "relativistic must be one of",
        ),
    ],
)
def test_bad_settings_are_refused(change, message):
    with pytest.raises(R.RecipeError, match=message):
        R.Recipe.from_dict(tier_a(**change))


@pytest.mark.parametrize(
    "change, message",
    [
        # The core-valence pair also differing in basis
        (
            {"cv_fc": {"level": "ORCA:QC@CCSD(T)/cc-pwCVQZ", "frozen_core": "valence"}},
            "differ in basis, frozen core",
        ),
        # The relativistic pair with different frozen cores
        (
            {"nonrel": {"level": "ORCA:QC@CCSD(T)/cc-pVTZ", "frozen_core": "inner"}},
            "differ in basis, frozen core, relativistic",
        ),
        # Subtracting a component from the same calculation
        ({"cv_fc": {"level": CV, "frozen_core": "inner"}}, "same calculation"),
    ],
)
def test_a_difference_must_change_one_thing(change, message):
    with pytest.raises(R.RecipeError, match=message):
        R.Recipe.from_dict(tier_a(**change))


def test_a_shipped_recipe_must_give_its_frozen_cores(tmp_path):
    import yaml

    data = tier_a(f12={"level": F12, "frozen_core": "default"})
    data["combine"] = "E0 = f12 + (cv_ae - cv_fc)\n"
    (tmp_path / "x.yaml").write_text(yaml.safe_dump(data))
    library = R.RecipeLibrary(shipped=False)
    library._add_directory(tmp_path, "shipped")
    assert "tier-a-core" not in library.names()
    assert "must give its frozen cores explicitly" in library.problems[0]
    # A user's recipe may leave it to the program
    library = R.RecipeLibrary(shipped=False)
    library._add_directory(tmp_path, "user")
    assert library.get("tier-a-core").uses_default_frozen_core


# -- the settings a component runs with ---------------------------------------------

CAPABILITIES = {
    "frozen_core": ["default", "none", "per-element"],
    "relativistic": ["none", "DKH2", "X2C"],
    "nucleus": ["point", "finite"],
    "decontract": True,
    "energy parts": True,
    "reference": ["default", "UHF", "QRO"],
}


@pytest.fixture
def resolved(monkeypatch):
    def resolve(level, available=None):
        return {
            "level": level,
            "step": "ORCA",
            "options": {"capabilities": dict(CAPABILITIES)},
        }

    monkeypatch.setattr(E, "resolve_level", resolve)
    return resolve


def test_component_settings_and_labels(resolved):
    recipe = R.Recipe.from_dict(tier_a())
    ae = E.component_mc(recipe, "cv_ae")
    fc = E.component_mc(recipe, "cv_fc")
    assert ae["options"]["frozen_core"] == SETS["inner"]
    assert fc["options"]["frozen_core"] == SETS["valence"]
    # One level, two calculations: told apart by their labels
    inner, valence = R.frozen_core_key(SETS["inner"]), R.frozen_core_key(
        SETS["valence"]
    )
    assert ae["level"] == f"{CV}|fc={inner}" and fc["level"] == f"{CV}|fc={valence}"
    rel = E.component_mc(recipe, "rel")
    assert rel["options"]["relativistic"] == "DKH2"
    assert rel["options"]["nucleus"] == "point"
    assert rel["level"].endswith(f"|fc={valence}|DKH2(point)")


def test_labels_are_canonical(resolved):
    """The same counts under another name give the same label (one calculation,
    one key); the same name with other counts gives another."""
    renamed = tier_a()
    renamed["frozen_core_definitions"] = {
        "core1": dict(reversed(list(SETS["inner"].items()))),  # other order too
        "valence": SETS["valence"],
    }
    renamed["components"]["cv_ae"]["frozen_core"] = "core1"
    a = E.component_mc(R.Recipe.from_dict(tier_a()), "cv_ae")["level"]
    b = E.component_mc(R.Recipe.from_dict(renamed), "cv_ae")["level"]
    assert a == b
    changed = tier_a()
    changed["frozen_core_definitions"] = {
        "inner": dict(SETS["inner"], S=0),
        "valence": SETS["valence"],
    }
    c = E.component_mc(R.Recipe.from_dict(changed), "cv_ae")["level"]
    assert c != a
    assert R.frozen_core_key("none") == "none"


def test_capabilities_are_checked(monkeypatch):
    recipe = R.Recipe.from_dict(tier_a())

    def without_x2c(level, available=None):
        caps = dict(CAPABILITIES, relativistic=["none"])
        return {"level": level, "step": "ORCA", "options": {"capabilities": caps}}

    monkeypatch.setattr(E, "resolve_level", without_x2c)
    with pytest.raises(R.RecipeError, match="no DKH2 Hamiltonian"):
        E.component_mc(recipe, "rel")

    def says_nothing(level, available=None):
        return {"level": level, "step": "X", "options": {}}

    monkeypatch.setattr(E, "resolve_level", says_nothing)
    with pytest.raises(R.RecipeError, match="does not say what settings"):
        E.component_mc(recipe, "cv_ae")


def test_decontracted_relativistic_pair(resolved):
    """A relativistic correction with the same decontracted basis in both halves
    (where no DK basis exists): one change, the Hamiltonian."""
    data = tier_a(
        rel={
            "level": "ORCA:QC@CCSD(T)/cc-pVTZ",
            "decontract": True,
            "frozen_core": "valence",
            "relativistic": "DKH2",
            "nucleus": "point",
        },
        nonrel={
            "level": "ORCA:QC@CCSD(T)/cc-pVTZ",
            "decontract": True,
            "frozen_core": "valence",
        },
    )
    recipe = R.Recipe.from_dict(data)
    rel = E.component_mc(recipe, "rel")
    assert rel["options"]["decontract"] is True
    assert "|decon|DKH2(point)" in rel["level"]
    # Decontracting the relativistic half only changes the basis with the
    # Hamiltonian, as the cc-pVTZ-DK / cc-pVTZ pair does: allowed
    data["components"]["nonrel"].pop("decontract")
    R.Recipe.from_dict(data)
    # ...but not with another frozen core as well
    data["components"]["nonrel"]["frozen_core"] = "inner"
    with pytest.raises(R.RecipeError, match="differ in basis, frozen core"):
        R.Recipe.from_dict(data)


# -- the open-shell reference ---------------------------------------------------------


def test_open_shell_reference(resolved):
    """The recipe's open-shell reference goes to every correlated component, not to
    the geometry level, and is not in the label (closed shells are unchanged)."""
    data = tier_a()
    data["open_shell_reference"] = "QRO"
    recipe = R.Recipe.from_dict(data)
    plain = R.Recipe.from_dict(tier_a())
    for name in ("f12", "cv_ae", "rel"):
        mc = E.component_mc(recipe, name)
        assert mc["options"]["reference"] == "QRO"
        assert mc["level"] == E.component_mc(plain, name)["level"]
    # The geometry level, when combine uses it, keeps the program's default
    data["combine"] = "E0 = geometry + (cv_ae - cv_fc)\n"
    with_geometry = R.Recipe.from_dict(data)
    assert "reference" not in E.component_mc(with_geometry, R.GEOMETRY)["options"]
    # ... nor in the recipe's own fingerprint: the components' fingerprints carry it
    assert recipe.fingerprint() == plain.fingerprint()


def test_a_component_can_override_the_reference(resolved):
    data = tier_a()
    data["open_shell_reference"] = "QRO"
    data["components"]["f12"] = dict(data["components"]["f12"], reference="UHF")
    recipe = R.Recipe.from_dict(data)
    assert E.component_mc(recipe, "f12")["options"]["reference"] == "UHF"
    assert E.component_mc(recipe, "cv_ae")["options"]["reference"] == "QRO"


@pytest.mark.parametrize(
    "change, message",
    [
        ({"open_shell_reference": "ROHF"}, "open_shell_reference must be one of"),
        ({"cv_fc_reference": "RHF"}, "reference must be one of"),
        ({"cv_fc_reference": "UHF"}, "differ in frozen core, reference"),
    ],
)
def test_bad_references_are_refused(change, message):
    data = tier_a()
    data["open_shell_reference"] = "QRO"
    if "open_shell_reference" in change:
        data["open_shell_reference"] = change["open_shell_reference"]
    if "cv_fc_reference" in change:
        comp = data["components"]["cv_fc"]
        data["components"]["cv_fc"] = dict(comp, reference=change["cv_fc_reference"])
    with pytest.raises(R.RecipeError, match=message):
        R.Recipe.from_dict(data)


def test_a_program_without_qro_is_refused(monkeypatch):
    data = tier_a()
    data["open_shell_reference"] = "QRO"
    recipe = R.Recipe.from_dict(data)

    def uhf_only(level, available=None):
        caps = dict(CAPABILITIES, reference=["default", "UHF"])
        return {"level": level, "step": "ORCA", "options": {"capabilities": caps}}

    monkeypatch.setattr(E, "resolve_level", uhf_only)
    with pytest.raises(R.RecipeError, match="cannot use a QRO reference"):
        E.component_mc(recipe, "f12")


def test_the_shipped_tier_a_recipes_use_qro():
    library = R.RecipeLibrary()
    for name in library.names():
        if "+CV+" in name:
            assert library.get(name).data["open_shell_reference"] == "QRO"
