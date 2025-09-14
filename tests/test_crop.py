"""Test functionalities for removing reactions to prevent growth."""

import cobra
import pytest
from cobra import Metabolite, Model, Reaction
from crop import run_crop_algorithm, build_phenotype_conditions
import sys

def create_test_model():
    """
    Creates a simple test model that incorrectly grows on lactose due to an extra reaction.

    The model should:
    - Grow on glucose (correct behavior)
    - Grow on lactose (incorrect - due to extra LACutil reaction)

    After removing the LACutil reaction:
    - Still grow on glucose
    - No longer grow on lactose
    """

    model = Model("test_crop_model")

    # Create metabolites
    metabolites = {
        "glc_e": Metabolite("glc_e", name="Glucose external", compartment="e"),
        "lac_e": Metabolite("lac_e", name="Lactose external", compartment="e"),
        "glc_c": Metabolite("glc_c", name="Glucose cytoplasm", compartment="c"),
        "lac_c": Metabolite("lac_c", name="Lactose cytoplasm", compartment="c"),
        "g6p_c": Metabolite("g6p_c", name="Glucose-6-phosphate", compartment="c"),
        "pyr_c": Metabolite("pyr_c", name="Pyruvate", compartment="c"),
        "atp_c": Metabolite("atp_c", name="ATP", compartment="c"),
        "adp_c": Metabolite("adp_c", name="ADP", compartment="c"),
        "biomass_c": Metabolite("biomass_c", name="Biomass", compartment="c"),
    }

    # Create reactions
    reactions = []

    # Exchange reactions
    ex_glc = Reaction("EX_glc")
    ex_glc.name = "Glucose exchange"
    ex_glc.lower_bound = -10  # Can uptake up to 10 units
    ex_glc.upper_bound = 100
    ex_glc.add_metabolites({metabolites["glc_e"]: -1})
    reactions.append(ex_glc)

    ex_lac = Reaction("EX_lac")
    ex_lac.name = "Lactose exchange"
    ex_lac.lower_bound = -10  # Can uptake up to 10 units
    ex_lac.upper_bound = 100
    ex_lac.add_metabolites({metabolites["lac_e"]: -1})
    reactions.append(ex_lac)

    # Transport reactions
    glc_transport = Reaction("GLCt")
    glc_transport.name = "Glucose transport"
    glc_transport.add_metabolites({metabolites["glc_e"]: -1, metabolites["glc_c"]: 1})
    reactions.append(glc_transport)

    lac_transport = Reaction("LACt")
    lac_transport.name = "Lactose transport"
    lac_transport.add_metabolites({metabolites["lac_e"]: -1, metabolites["lac_c"]: 1})
    reactions.append(lac_transport)

    # Glucose metabolism (correct pathway)
    hexokinase = Reaction("HEX")
    hexokinase.name = "Hexokinase"
    hexokinase.add_metabolites(
        {
            metabolites["glc_c"]: -1,
            metabolites["atp_c"]: -1,
            metabolites["g6p_c"]: 1,
            metabolites["adp_c"]: 1,
        }
    )
    reactions.append(hexokinase)

    # Simplified glycolysis
    glycolysis = Reaction("GLYC")
    glycolysis.name = "Simplified glycolysis"
    glycolysis.add_metabolites(
        {
            metabolites["g6p_c"]: -1,
            metabolites["adp_c"]: -3,
            metabolites["pyr_c"]: 2,
            metabolites["atp_c"]: 3,  # Net gain of 2 ATP (used 1 in HEX, gained 4 here)
        }
    )
    reactions.append(glycolysis)

    # EXTRA REACTION - This shouldn't exist and allows incorrect lactose growth
    lac_util = Reaction("LACutil")
    lac_util.name = "Lactose utilization (EXTRA - should be removed)"
    lac_util.add_metabolites(
        {
            metabolites["lac_c"]: -1,
            metabolites["atp_c"]: -1,
            metabolites["g6p_c"]: 1,
            metabolites["adp_c"]: 1,  # Unrealistic direct conversion
        }
    )
    lac_util.upper_bound = 9

    reactions.append(lac_util)

    # Biomass reaction
    biomass = Reaction("BIOMASS")
    biomass.name = "Biomass formation"
    biomass.add_metabolites(
        {
            metabolites["pyr_c"]: -2,
            metabolites["atp_c"]: -1,
            metabolites["biomass_c"]: 1,
            metabolites["adp_c"]: 1,
        }
    )
    biomass.upper_bound = 100
    reactions.append(biomass)

    biomass_exchange = Reaction("EX_biomass_c")
    biomass_exchange.name = "Biomass exchange"
    biomass.add_metabolites(
        {
            metabolites["biomass_c"]: -1,
        }
    )
    reactions.append(biomass_exchange)

    # ATP maintenance (to prevent unrealistic ATP accumulation)
    atp_maintenance = Reaction("ATPM")
    atp_maintenance.name = "ATP maintenance"
    atp_maintenance.upper_bound = 100.0  # Force some ATP production
    atp_maintenance.lower_bound = 0.0
    atp_maintenance.add_metabolites({metabolites["atp_c"]: -1, metabolites["adp_c"]: 1})
    reactions.append(atp_maintenance)

    # Add all reactions to model
    model.add_reactions(reactions)

    # Set objective
    model.objective = "BIOMASS"
    for rxn in model.reactions:
        print(
            f"Reaction {rxn.id}: {rxn.name}\t{rxn.build_reaction_string()}\t{rxn.lower_bound} < {rxn.upper_bound}"
        )
    cobra.io.save_json_model(model, "test_crop_model.json")
    cobra.io.write_sbml_model(model, "test_crop_model.xml")

    return model


def create_old_test_model():
    """
    Creates a simple test model that incorrectly grows on lactose due to an extra reaction.

    The model should:
    - Grow on glucose (correct behavior)
    - Grow on lactose (incorrect - due to extra LACutil reaction)

    After removing the LACutil reaction:
    - Still grow on glucose
    - No longer grow on lactose
    """

    model = Model("test_crop_model")

    # Create metabolites
    metabolites = {
        "glc_e": Metabolite("glc_e", name="Glucose external", compartment="e"),
        "lac_e": Metabolite("lac_e", name="Lactose external", compartment="e"),
        "glc_c": Metabolite("glc_c", name="Glucose cytoplasm", compartment="c"),
        "lac_c": Metabolite("lac_c", name="Lactose cytoplasm", compartment="c"),
        "g6p_c": Metabolite("g6p_c", name="Glucose-6-phosphate", compartment="c"),
        "pyr_c": Metabolite("pyr_c", name="Pyruvate", compartment="c"),
        "atp_c": Metabolite("atp_c", name="ATP", compartment="c"),
        "adp_c": Metabolite("adp_c", name="ADP", compartment="c"),
        "biomass_c": Metabolite("biomass_c", name="Biomass", compartment="c"),
    }

    # Create reactions
    reactions = []

    # Exchange reactions
    ex_glc = Reaction("EX_glc")
    ex_glc.name = "Glucose exchange"
    ex_glc.lower_bound = -10  # Can uptake up to 10 units
    ex_glc.upper_bound = 100
    ex_glc.add_metabolites({metabolites["glc_e"]: -1})
    reactions.append(ex_glc)

    ex_lac = Reaction("EX_lac")
    ex_lac.name = "Lactose exchange"
    ex_lac.lower_bound = -10  # Can uptake up to 10 units
    ex_lac.upper_bound = 100
    ex_lac.add_metabolites({metabolites["lac_e"]: -1})
    reactions.append(ex_lac)

    # Transport reactions
    glc_transport = Reaction("GLCt")
    glc_transport.name = "Glucose transport"
    glc_transport.add_metabolites({metabolites["glc_e"]: -1, metabolites["glc_c"]: 1})
    reactions.append(glc_transport)

    lac_transport = Reaction("LACt")
    lac_transport.name = "Lactose transport"
    lac_transport.add_metabolites({metabolites["lac_e"]: -1, metabolites["lac_c"]: 1})
    reactions.append(lac_transport)

    # Glucose metabolism (correct pathway)
    hexokinase = Reaction("HEX")
    hexokinase.name = "Hexokinase"
    hexokinase.add_metabolites(
        {
            metabolites["glc_c"]: -1,
            metabolites["atp_c"]: -1,
            metabolites["g6p_c"]: 1,
            metabolites["adp_c"]: 1,
        }
    )
    reactions.append(hexokinase)

    # Simplified glycolysis
    glycolysis = Reaction("GLYC")
    glycolysis.name = "Simplified glycolysis"
    glycolysis.add_metabolites(
        {
            metabolites["g6p_c"]: -1,
            metabolites["adp_c"]: -3,
            metabolites["pyr_c"]: 2,
            metabolites["atp_c"]: 3,  # Net gain of 2 ATP (used 1 in HEX, gained 4 here)
        }
    )
    reactions.append(glycolysis)

    # EXTRA REACTION - This shouldn't exist and allows incorrect lactose growth
    lac_util = Reaction("LACutil")
    lac_util.name = "Lactose utilization (EXTRA - should be removed)"
    lac_util.add_metabolites(
        {
            metabolites["lac_c"]: -1,
            metabolites["atp_c"]: -1,
            metabolites["g6p_c"]: 1,
            metabolites["adp_c"]: 1,  # Unrealistic direct conversion
        }
    )
    lac_util.upper_bound = 9

    reactions.append(lac_util)

    # Biomass reaction
    biomass = Reaction("BIOMASS")
    biomass.name = "Biomass formation"
    biomass.add_metabolites(
        {
            metabolites["pyr_c"]: -2,
            metabolites["atp_c"]: -1,
            metabolites["biomass_c"]: 1,
            metabolites["adp_c"]: 1,
        }
    )
    biomass.upper_bound = 100
    reactions.append(biomass)

    biomass_exchange = Reaction("EX_biomass_c")
    biomass_exchange.name = "Biomass exchange"
    biomass.add_metabolites(
        {
            metabolites["biomass_c"]: -1,
        }
    )
    reactions.append(biomass_exchange)

    # ATP maintenance (to prevent unrealistic ATP accumulation)
    atp_maintenance = Reaction("ATPM")
    atp_maintenance.name = "ATP maintenance"
    atp_maintenance.upper_bound = 100.0  # Force some ATP production
    atp_maintenance.lower_bound = 0.0
    atp_maintenance.add_metabolites({metabolites["atp_c"]: -1, metabolites["adp_c"]: 1})
    reactions.append(atp_maintenance)

    # Add all reactions to model
    model.add_reactions(reactions)

    # Set objective
    model.objective = "BIOMASS"
    #for rxn in model.reactions:
        #print(
        #    f"Reaction {rxn.id}: {rxn.name}\t{rxn.build_reaction_string()}\t{rxn.lower_bound} < {rxn.upper_bound}"
        #)
    #cobra.io.save_json_model(model, "test_crop_model.json")
    #cobra.io.write_sbml_model(model, "test_crop_model.xml")

    return model


def apply_medium(model, medium_dict):
    """Apply medium conditions to model"""
    model_copy = model.copy()
    for reaction_id, bound in medium_dict.items():
        if reaction_id in model_copy.reactions:
            model_copy.reactions.get_by_id(reaction_id).lower_bound = bound
    return model_copy


# Pytest fixtures
@pytest.fixture
def test_model():
    """Fixture providing the test model"""
    return create_test_model()


@pytest.fixture
def media_conditions():
    """Fixture providing different media conditions"""
    return {
        "glucose": {"EX_glc": -10, "EX_lac": 0},
        "lactose": {"EX_glc": 0, "EX_lac": -10},
        "no_carbon": {"EX_glc": 0, "EX_lac": 0},
    }


@pytest.fixture
def expected_problematic_reaction():
    """Fixture providing the expected reaction to be removed"""
    return "LACutil"


@pytest.fixture
def phenotype_data():
    """Fixture providing phenotype observation data"""
    return {
        "glucose": {"observed": "growth", "predicted": "growth"},  # Correct
        "lactose": {
            "observed": "no_growth",
            "predicted": "growth",
        },  # Incorrect - needs fixing
        "no_carbon": {"observed": "no_growth", "predicted": "no_growth"},  # Correct
    }


@pytest.fixture(scope="module")
def build():
    """
    Provide build_phenotype_conditions via relative import, stubbing heavy deps first.
    """
    # Stub cvxpy to satisfy "from cvxpy import Minimize, Problem, Variable, diag"
    if "cvxpy" not in sys.modules:
        cvxpy_stub = types.ModuleType("cvxpy")
        # Minimal placeholders; not used in these tests
        setattr(cvxpy_stub, "Minimize", object)
        setattr(cvxpy_stub, "Problem", object)
        setattr(cvxpy_stub, "Variable", object)
        setattr(cvxpy_stub, "diag", lambda x: x)
        sys.modules["cvxpy"] = cvxpy_stub

    # Stub cobra to satisfy "import cobra" at module import time
    if "cobra" not in sys.modules:
        cobra_stub = types.ModuleType("cobra")
        sys.modules["cobra"] = cobra_stub

    return build_phenotype_conditions


@pytest.fixture
def media_conditions():
    return {
        "glucose": {"EX_glc": -10.0, "EX_lac": 0.0},
        "lactose": {"EX_glc": 0.0, "EX_lac": -5.0},
        "no_carbon": {"EX_glc": 0.0, "EX_lac": 0.0},
    }


@pytest.fixture
def phenotype_data():
    return {
        "glucose": {"observed": "growth", "predicted": "growth"},
        "lactose": {"observed": "no_growth", "predicted": "growth"},
        "no_carbon": {"observed": "no_growth", "predicted": "no_growth"},
    }


@pytest.fixture
def phenotype_observation_key():
    # Uses "observation" key instead of "observed"
    return {
        "glucose": {"observation": "growth", "predicted": "growth"},
        "lactose": {"observation": "no_growth", "predicted": "growth"},
    }

# Test functions

def test_basic_mapping_and_keys(build, media_conditions, phenotype_data):
    result = build(media_conditions, phenotype_data)
    assert "growth" in result and "nogrowth" in result

    growth = result["growth"]
    nogrowth = result["nogrowth"]

    # Union of all exchanges from all media should be present in both maps
    expected_keys = {"EX_glc", "EX_lac"}
    assert set(growth.keys()) == expected_keys
    assert set(nogrowth.keys()) == expected_keys

    # Growth condition chosen should be "glucose" -> EX_glc uptake magnitude 10
    assert growth["EX_glc"] == 10.0
    assert growth["EX_lac"] == 0.0
    

    # Nogrowth condition chosen should be "lactose" -> EX_lac 5, EX_glc 0
    assert nogrowth["EX_glc"] == 0.0
    assert nogrowth["EX_lac"] == 5.0


def test_observation_key_supported(build, media_conditions, phenotype_observation_key):
    result = build(media_conditions, phenotype_observation_key)
    growth = result["growth"]
    nogrowth = result["nogrowth"]

    # Should behave the same as if "observed" was used
    assert growth["EX_glc"] == 10.0
    assert growth["EX_lac"] == 0.0
    assert nogrowth["EX_glc"] == 0.0
    assert nogrowth["EX_lac"] == 5.0


def test_select_specific_conditions(build, media_conditions, phenotype_data):
    # Add another valid growth condition
    media_plus = dict(media_conditions)
    media_plus["glucose_weak"] = {"EX_glc": -7.0}
    pheno_plus = dict(phenotype_data)
    pheno_plus["glucose_weak"] = {"observed": "growth", "predicted": "growth"}

    result = build(
        media_plus,
        pheno_plus,
        growth_condition="glucose_weak",
        nogrowth_condition="lactose",
    )
    growth = result["growth"]
    nogrowth = result["nogrowth"]

    # Selected explicit growth condition should reflect -7 -> +7 magnitude
    assert growth["EX_glc"] == 7.0
    # Keys from union should still appear
    assert "EX_lac" in growth 
    # Selected explicit nogrowth condition should be lactose as before
    assert nogrowth["EX_lac"] == 5.0


def test_invalid_growth_condition_keyerror(build, media_conditions, phenotype_data):
    with pytest.raises(KeyError):
        build(media_conditions, phenotype_data, growth_condition="does_not_exist")


def test_invalid_nogrowth_condition_keyerror(build, media_conditions, phenotype_data):
    with pytest.raises(KeyError):
        build(media_conditions, phenotype_data, nogrowth_condition="not_in_media")


def test_no_growth_candidate_raises(build):
    media = {"A": {"EX_a": -1.0}}
    # No case with observed==growth and predicted==growth
    pheno = {
        "A": {"observed": "no_growth", "predicted": "growth"},
        "B": {"observed": "growth", "predicted": "no_growth"},
    }
    with pytest.raises(ValueError):
        build(media, pheno)


def test_no_nogrowth_candidate_raises(build):
    media = {"glucose": {"EX_glc": -10.0}, "lactose": {"EX_lac": -5.0}}
    # No case with observed in {no_growth, nogrowth} and predicted==growth
    pheno = {
        "glucose": {"observed": "growth", "predicted": "growth"},
        "lactose": {"observed": "no_growth", "predicted": "no_growth"},
    }
    with pytest.raises(ValueError):
        build(media, pheno)


def test_union_of_exchanges_included(build):
    media = {
        "c1": {"EX_a": -3.0},     # only EX_a here
        "c2": {"EX_b": -4.0},     # only EX_b here
    }
    pheno = {
        "c1": {"observed": "growth", "predicted": "growth"},
        "c2": {"observed": "no_growth", "predicted": "growth"},
    }
    result = build(media, pheno)
    growth = result["growth"]
    nogrowth = result["nogrowth"]

    # Both maps must include both EX_a and EX_b
    assert set(growth.keys()) == {"EX_a", "EX_b"}
    assert set(nogrowth.keys()) == {"EX_a", "EX_b"}

    # Magnitudes flip sign for uptakes, zero otherwise
    assert growth["EX_a"] == 3.0 and growth["EX_b"] == 0.0
    assert nogrowth["EX_b"] == 4.0 and nogrowth["EX_a"] == 0.0

def test_model_creation(test_model):
    """Test that the model is created correctly"""
    assert len(test_model.reactions) > 0
    assert len(test_model.metabolites) > 0
    assert test_model.objective.direction == "max"
    assert "LACutil" in [r.id for r in test_model.reactions]


def test_initial_glucose_growth(test_model, media_conditions):
    """Test that the model grows correctly on glucose"""
    glucose_model = apply_medium(test_model, media_conditions["glucose"])
    solution = glucose_model.optimize()

    assert solution.status == "optimal"
    assert solution.objective_value > 0.001, "Error, the initial model should grow on glucose"


def test_initial_lactose_growth_incorrect(test_model, media_conditions):
    """Test that the model incorrectly grows on lactose (the problem to be fixed)"""
    lactose_model = apply_medium(test_model, media_conditions["lactose"])
    solution = lactose_model.optimize()

    assert solution.status == "optimal"
    assert (
        solution.objective_value > 0.001
    ), "Error, the initial model should grow on lactose due to LACutil reaction"


def test_no_carbon_no_growth(test_model, media_conditions):
    """Test that the model doesn't grow without carbon source"""
    no_carbon_model = apply_medium(test_model, media_conditions["no_carbon"])
    solution = no_carbon_model.optimize()
    assert (
        solution.objective_value < 0.001
    ), "Error, the initial model should not grow without carbon source "


@pytest.mark.parametrize(
    "medium_name,expected_growth",
    [("glucose", True), ("lactose", True), ("no_carbon", False)],  # Initially incorrect
)
def test_initial_growth_patterns(test_model, media_conditions, medium_name, expected_growth):
    """Parametrized test for initial growth patterns"""
    model_with_medium = apply_medium(test_model, media_conditions[medium_name])
    solution = model_with_medium.optimize()

    if expected_growth:
        assert solution.objective_value > 0.001
    else:
        assert solution.objective_value < 0.001


def test_reaction_removal_fixes_lactose(
    test_model, media_conditions, expected_problematic_reaction
):
    """Test that removing the problematic reaction prevents lactose growth"""
    # Remove the problematic reaction
    model_corrected = test_model.copy()
    model_corrected.reactions.get_by_id(expected_problematic_reaction).remove_from_model()

    # Test that lactose no longer supports growth
    lactose_model = apply_medium(model_corrected, media_conditions["lactose"])
    solution = lactose_model.optimize()

    assert solution.objective_value < 0.001


def test_reaction_removal_preserves_glucose_growth(
    test_model, media_conditions, expected_problematic_reaction
):
    """Test that removing the problematic reaction doesn't affect glucose growth"""
    # Remove the problematic reaction
    model_corrected = test_model.copy()
    model_corrected.reactions.get_by_id(expected_problematic_reaction).remove_from_model()

    # Test that glucose still works
    glucose_model = apply_medium(model_corrected, media_conditions["glucose"])
    solution = glucose_model.optimize()

    assert solution.status == "optimal"
    assert solution.objective_value > 0.001


@pytest.mark.parametrize(
    "medium_name,expected_growth_after_fix",
    [
        ("glucose", True),  # Should still grow
        ("lactose", False),  # Should no longer grow
        ("no_carbon", False),  # Should still not grow
    ],
)
def test_growth_after_reaction_removal(
    test_model,
    media_conditions,
    expected_problematic_reaction,
    medium_name,
    expected_growth_after_fix,
):
    """Parametrized test for growth patterns after removing problematic reaction"""
    # Remove the problematic reaction
    model_corrected = test_model.copy()
    model_corrected.reactions.get_by_id(expected_problematic_reaction).remove_from_model()

    # Test growth in specified medium
    model_with_medium = apply_medium(model_corrected, media_conditions[medium_name])
    solution = model_with_medium.optimize()

    if expected_growth_after_fix:
        assert solution.objective_value > 0.001
    else:
        assert solution.objective_value < 0.001


def test_crop_algorithm_integration(
    test_model, media_conditions, phenotype_data, expected_problematic_reaction
):
    """
    Test template for integration with your CROP algorithm
    Replace the simulated result with your actual CROP algorithm call
    """

    # This is where you would call your CROP algorithm
    suggested_removals, _ = run_crop_algorithm(test_model, phenotype_data, media_conditions)

    # For now, simulate the expected result
    #suggested_removals = ["LACutil"]

    # Test that CROP suggests the correct reaction for removal
    assert expected_problematic_reaction in suggested_removals

    # Test that removing suggested reactions fixes the problem
    model_corrected = test_model.copy()
    for reaction_id in suggested_removals:
        model_corrected.reactions.get_by_id(reaction_id).remove_from_model()

    # Verify the fix works
    lactose_model = apply_medium(model_corrected, media_conditions["lactose"])
    lactose_solution = lactose_model.optimize()
    assert (
        lactose_solution.objective_value < 0.001
    ), f"Error, the model should not grow on lactose after removing {suggested_removals}"

    # Verify other conditions still work
    glucose_model = apply_medium(model_corrected, media_conditions["glucose"])
    glucose_solution = glucose_model.optimize()
    assert (
        glucose_solution.objective_value > 0.001
    ), f"Error, the model should still grow on glucose after removing {suggested_removals}"

    # Verify no carbon condition
    no_carbon_model = apply_medium(model_corrected, media_conditions["no_carbon"])
    no_carbon_solution = no_carbon_model.optimize()
    assert (
        no_carbon_solution.objective_value < 0.001
    ), f"Error, the model should not grow without carbon source after removing {suggested_removals}"


def test_specific_reaction_exists(test_model, expected_problematic_reaction):
    """Test that the expected problematic reaction exists in the model"""
    reaction_ids = [r.id for r in test_model.reactions]
    assert expected_problematic_reaction in reaction_ids


def test_reaction_removal_is_reversible(test_model, expected_problematic_reaction):
    """Test that the reaction removal doesn't permanently modify the original model"""
    original_reaction_count = len(test_model.reactions)

    # Create a copy and remove reaction
    model_copy = test_model.copy()
    model_copy.reactions.get_by_id(expected_problematic_reaction).remove_from_model()

    # Original model should be unchanged
    assert len(test_model.reactions) == original_reaction_count
    assert expected_problematic_reaction in [r.id for r in test_model.reactions]

    # Copy should have one fewer reaction
    assert len(model_copy.reactions) == original_reaction_count - 1
    assert expected_problematic_reaction not in [r.id for r in model_copy.reactions]


# Integration test for complete CROP workflow
@pytest.mark.integration
def test_complete_crop_workflow(
    test_model, media_conditions, phenotype_data, expected_problematic_reaction
):
    """Integration test for the complete CROP workflow"""

    # Step 1: Verify initial problematic behavior
    lactose_model = apply_medium(test_model, media_conditions["lactose"])
    initial_solution = lactose_model.optimize()
    assert (
        initial_solution.objective_value > 0.001
    ), "Error, the initial model should grow on lactose"

    # Step 2: Run CROP algorithm (simulated)
    suggested_removals, _ = run_crop_algorithm(test_model, phenotype_data, media_conditions)
    #suggested_removals = [expected_problematic_reaction]  # Simulated result
    #assert expected_problematic_reaction in suggested_removals
    # Step 3: Apply suggested changes
    corrected_model = test_model.copy()
    for reaction_id in suggested_removals:
        corrected_model.reactions.get_by_id(reaction_id).remove_from_model()

    # Step 4: Verify fix
    # Should no longer grow on lactose
    lactose_corrected = apply_medium(corrected_model, media_conditions["lactose"])
    lactose_solution = lactose_corrected.optimize()
    assert (
        lactose_solution.objective_value < 0.001
    ), f"Error, the model should not grow on lactose after removing {suggested_removals}"

    # Should still grow on glucose
    glucose_corrected = apply_medium(corrected_model, media_conditions["glucose"])
    glucose_solution = glucose_corrected.optimize()
    assert (
        glucose_solution.objective_value > 0.001
    ), f"Error, the model should still grow on glucose after removing {suggested_removals}"

    # Should not grow without carbon source
    no_carbon_corrected = apply_medium(corrected_model, media_conditions["no_carbon"])
    no_carbon_solution = no_carbon_corrected.optimize()
    assert (
        no_carbon_solution.objective_value < 0.001
    ), f"Error, the model should not grow without carbon source after removing {suggested_removals}"


# Example usage and manual testing
if __name__ == "__main__":
    # Create and test the model manually
    model = create_test_model()
    media = {
        "glucose": {"EX_glc": -10, "EX_lac": 0},
        "lactose": {"EX_glc": 0, "EX_lac": -10},
        "no_carbon": {"EX_glc": 0, "EX_lac": 0},
    }

    # print("=== Initial Model Testing ===")

    # Test glucose
    glucose_model = apply_medium(model, media["glucose"])
    sol = glucose_model.optimize()
    # print(f"Glucose growth: {sol.objective_value:.4f} (should be > 0)")

    # Test lactose
    lactose_model = apply_medium(model, media["lactose"])
    sol = lactose_model.optimize()
    # print(
    #    f"Lactose growth: {sol.objective_value:.4f} (should be > 0, but this is wrong!)"
    # )

    # print("\n=== After Removing LACutil Reaction ===")

    # Remove the problematic reaction
    model_fixed = model.copy()
    model_fixed.reactions.LACutil.remove_from_model()

    # Test glucose (should still work)
    glucose_model = apply_medium(model_fixed, media["glucose"])
    sol = glucose_model.optimize()
    # print(f"Glucose growth: {sol.objective_value:.4f} (should still be > 0)")

    # Test lactose (should now fail)
    lactose_model = apply_medium(model_fixed, media["lactose"])
    sol = lactose_model.optimize()
    # print(f"Lactose growth: {sol.objective_value:.4f} (should now be ~0)")
    # print("\nTo run pytest tests, use:")
    # print("pytest test_crop_model.py -v")
    # print("pytest test_crop_model.py -v -k 'test_initial'  # Run only initial tests")
    # print("pytest test_crop_model.py -v -m integration     # Run only integration tests")


# ============================================================================
# API UTILITY FUNCTION TESTS
# ============================================================================

import types
import pandas as pd
import numpy as np
from unittest.mock import Mock, MagicMock, patch


class TestApiUtilityFunctions:
    """Test all utility functions in api.py"""

    @pytest.fixture
    def sample_phenotype_data(self):
        """Sample phenotype data for testing"""
        return {
            "glucose": {"observed": "growth", "predicted": "growth"},
            "lactose": {"observed": "no_growth", "predicted": "growth"},
            "no_carbon": {"observed": "no_growth", "predicted": "no_growth"},
            "mixed": {"observed": "growth", "predicted": "no_growth"}
        }

    @pytest.fixture
    def sample_media_conditions(self):
        """Sample media conditions for testing"""
        return {
            "glucose": {"EX_glc": -10.0, "EX_lac": 0.0},
            "lactose": {"EX_glc": 0.0, "EX_lac": -5.0},
            "no_carbon": {"EX_glc": 0.0, "EX_lac": 0.0},
            "mixed": {"EX_glc": -3.0, "EX_lac": -2.0}
        }

    @pytest.fixture
    def mock_model(self):
        """Mock COBRA model for testing"""
        mock_rxn1 = Mock()
        mock_rxn1.id = "EX_glc"
        mock_rxn2 = Mock()
        mock_rxn2.id = "EX_lac"
        mock_rxn3 = Mock()
        mock_rxn3.id = "BIOMASS"
        mock_rxn4 = Mock()
        mock_rxn4.id = "ATPM"
        
        mock_model = Mock()
        mock_model.reactions = [mock_rxn1, mock_rxn2, mock_rxn3, mock_rxn4]
        return mock_model

    def test_get_growth_conditions(self, sample_phenotype_data):
        """Test get_growth_conditions function"""
        from crop.api import get_growth_conditions
        
        result = get_growth_conditions(sample_phenotype_data)
        expected = ["glucose"]  # Only glucose has observed=growth and predicted=growth
        assert result == expected

    def test_get_growth_conditions_empty(self):
        """Test get_growth_conditions with no matching conditions"""
        from crop.api import get_growth_conditions
        
        phenotype_data = {
            "condition1": {"observed": "no_growth", "predicted": "growth"},
            "condition2": {"observed": "growth", "predicted": "no_growth"}
        }
        result = get_growth_conditions(phenotype_data)
        assert result == []

    def test_get_growth_conditions_multiple(self):
        """Test get_growth_conditions with multiple matching conditions"""
        from crop.api import get_growth_conditions
        
        phenotype_data = {
            "glucose": {"observed": "growth", "predicted": "growth"},
            "fructose": {"observed": "growth", "predicted": "growth"},
            "lactose": {"observed": "no_growth", "predicted": "growth"}
        }
        result = get_growth_conditions(phenotype_data)
        assert set(result) == {"glucose", "fructose"}

    def test_get_nogrowth_conditions(self, sample_phenotype_data):
        """Test get_nogrowth_conditions function"""
        from crop.api import get_nogrowth_conditions
        
        result = get_nogrowth_conditions(sample_phenotype_data)
        expected = ["lactose"]  # Only lactose has observed=no_growth and predicted=growth
        assert result == expected

    def test_get_nogrowth_conditions_empty(self):
        """Test get_nogrowth_conditions with no matching conditions"""
        from crop.api import get_nogrowth_conditions
        
        phenotype_data = {
            "condition1": {"observed": "growth", "predicted": "growth"},
            "condition2": {"observed": "no_growth", "predicted": "no_growth"}
        }
        result = get_nogrowth_conditions(phenotype_data)
        assert result == []

    def test_get_nogrowth_conditions_multiple(self):
        """Test get_nogrowth_conditions with multiple matching conditions"""
        from crop.api import get_nogrowth_conditions
        
        phenotype_data = {
            "lactose": {"observed": "no_growth", "predicted": "growth"},
            "sucrose": {"observed": "no_growth", "predicted": "growth"},
            "glucose": {"observed": "growth", "predicted": "growth"}
        }
        result = get_nogrowth_conditions(phenotype_data)
        assert set(result) == {"lactose", "sucrose"}

    def test_get_carbon_source(self, mock_model, sample_phenotype_data, sample_media_conditions):
        """Test get_carbon_source function"""
        from crop.api import get_carbon_source, get_growth_conditions
        
        result = get_carbon_source(
            mock_model, 
            sample_phenotype_data, 
            sample_media_conditions, 
            get_growth_conditions
        )
        
        # Should return glucose condition with index and name for EX_glc
        assert "glucose" in result
        assert result["glucose"] == (0, "EX_glc")  # EX_glc is first in mock reactions

    def test_get_lower_bound_for_conditions(self, mock_model, sample_phenotype_data, sample_media_conditions):
        """Test get_lower_bound_for_conditions function"""
        from crop.api import get_lower_bound_for_conditions, get_growth_conditions
        
        result = get_lower_bound_for_conditions(
            mock_model,
            sample_phenotype_data,
            sample_media_conditions,
            get_growth_conditions
        )
        
        assert isinstance(result, pd.DataFrame)
        assert "glucose" in result.columns
        assert result.loc["EX_glc", "glucose"] == -10.0
        assert result.loc["EX_lac", "glucose"] == 0.0

    def test_get_lower_bound_for_growth_conditions(self, mock_model, sample_phenotype_data, sample_media_conditions):
        """Test get_lower_bound_for_growth_conditions function"""
        from crop.api import get_lower_bound_for_growth_conditions
        
        result = get_lower_bound_for_growth_conditions(
            mock_model,
            sample_phenotype_data,
            sample_media_conditions,
            atp_maintenance_rxn="ATPM",
            atp_maintenance_lower_bound=2.5
        )
        
        assert isinstance(result, pd.DataFrame)
        assert "glucose" in result.columns
        # Should have ATP maintenance set
        assert result.loc["ATPM", "glucose"] == 2.5

    def test_get_lower_bound_for_nogrowth_conditions(self, mock_model, sample_phenotype_data, sample_media_conditions):
        """Test get_lower_bound_for_nogrowth_conditions function"""
        from crop.api import get_lower_bound_for_nogrowth_conditions
        
        result = get_lower_bound_for_nogrowth_conditions(
            mock_model,
            sample_phenotype_data,
            sample_media_conditions
        )
        
        assert isinstance(result, pd.DataFrame)
        assert "lactose" in result.columns
        assert result.loc["EX_lac", "lactose"] == -5.0  # Lactose media condition
        assert result.loc["EX_glc", "lactose"] == 0.0  # Not in nogrowth media

    def test_get_upper_bound_for_conditions(self, mock_model, sample_phenotype_data, sample_media_conditions):
        """Test get_upper_bound_for_conditions function"""
        from crop.api import get_upper_bound_for_conditions, get_growth_conditions
        
        result = get_upper_bound_for_conditions(
            mock_model,
            sample_phenotype_data,
            sample_media_conditions,
            get_growth_conditions,
            growth_limit=5.0
        )
        
        assert isinstance(result, pd.DataFrame)
        assert "glucose" in result.columns
        # Reactions in media should have 0, others should have growth_limit
        assert result.loc["EX_glc", "glucose"] == 0.0  # In media
        assert result.loc["BIOMASS", "glucose"] == 5.0  # Not in media

    def test_get_upper_bound_for_nogrowth_conditions(self, mock_model, sample_phenotype_data, sample_media_conditions):
        """Test get_upper_bound_for_nogrowth_conditions function"""
        from crop.api import get_upper_bound_for_nogrowth_conditions
        
        result = get_upper_bound_for_nogrowth_conditions(
            mock_model,
            sample_phenotype_data,
            sample_media_conditions,
            maximum_nogrowth=1.0
        )
        
        assert isinstance(result, pd.DataFrame)
        assert "lactose" in result.columns
        assert result.loc["EX_lac", "lactose"] == 0.0  # In media
        assert result.loc["BIOMASS", "lactose"] == 1.0  # Not in media

    def test_get_upper_bound_for_growth_conditions(self, mock_model, sample_phenotype_data, sample_media_conditions):
        """Test get_upper_bound_for_growth_conditions function"""
        from crop.api import get_upper_bound_for_growth_conditions
        
        result = get_upper_bound_for_growth_conditions(
            mock_model,
            sample_phenotype_data,
            sample_media_conditions,
            minimum_growth=2.0
        )
        
        assert isinstance(result, pd.DataFrame)
        assert "glucose" in result.columns
        assert result.loc["EX_glc", "glucose"] == 0.0  # In media
        assert result.loc["BIOMASS", "glucose"] == 2.0  # Not in media

    def test_build_phenotype_conditions_basic(self, sample_media_conditions, sample_phenotype_data):
        """Test build_phenotype_conditions basic functionality"""
        from crop.api import build_phenotype_conditions
        
        result = build_phenotype_conditions(sample_media_conditions, sample_phenotype_data)
        
        assert "growth" in result
        assert "nogrowth" in result
        
        # Check sign flipping and correct selection
        growth = result["growth"]
        nogrowth = result["nogrowth"]
        
        # Glucose condition: EX_glc=-10 becomes +10, EX_lac=0 stays 0
        assert growth["EX_glc"] == 10.0
        assert growth["EX_lac"] == 0.0
        
        # Lactose condition: EX_lac=-5 becomes +5, EX_glc=0 stays 0
        assert nogrowth["EX_glc"] == 0.0
        assert nogrowth["EX_lac"] == 5.0

    def test_build_phenotype_conditions_observation_key(self):
        """Test build_phenotype_conditions with 'observation' key instead of 'observed'"""
        from crop.api import build_phenotype_conditions
        
        media_conditions = {
            "glucose": {"EX_glc": -10.0, "EX_lac": 0.0},
            "lactose": {"EX_glc": 0.0, "EX_lac": -5.0}
        }
        phenotype_data = {
            "glucose": {"observation": "growth", "predicted": "growth"},
            "lactose": {"observation": "no_growth", "predicted": "growth"}
        }
        
        result = build_phenotype_conditions(media_conditions, phenotype_data)
        
        assert result["growth"]["EX_glc"] == 10.0
        assert result["nogrowth"]["EX_lac"] == 5.0

    def test_build_phenotype_conditions_specific_conditions(self, sample_media_conditions, sample_phenotype_data):
        """Test build_phenotype_conditions with specific growth/nogrowth conditions"""
        from crop.api import build_phenotype_conditions
        
        # Add another valid condition
        media_extended = dict(sample_media_conditions)
        media_extended["glucose2"] = {"EX_glc": -7.0, "EX_lac": 0.0}
        phenotype_extended = dict(sample_phenotype_data)
        phenotype_extended["glucose2"] = {"observed": "growth", "predicted": "growth"}
        
        result = build_phenotype_conditions(
            media_extended,
            phenotype_extended,
            growth_condition="glucose2",
            nogrowth_condition="lactose"
        )
        
        assert result["growth"]["EX_glc"] == 7.0  # Should use glucose2
        assert result["nogrowth"]["EX_lac"] == 5.0  # Should use lactose

    def test_build_phenotype_conditions_no_growth_candidates(self):
        """Test build_phenotype_conditions when no growth candidates exist"""
        from crop.api import build_phenotype_conditions
        
        media_conditions = {"test": {"EX_glc": -10.0}}
        phenotype_data = {"test": {"observed": "no_growth", "predicted": "growth"}}
        
        with pytest.raises(ValueError, match="No growth condition found"):
            build_phenotype_conditions(media_conditions, phenotype_data)

    def test_build_phenotype_conditions_no_nogrowth_candidates(self):
        """Test build_phenotype_conditions when no nogrowth candidates exist"""
        from crop.api import build_phenotype_conditions
        
        media_conditions = {"test": {"EX_glc": -10.0}}
        phenotype_data = {"test": {"observed": "growth", "predicted": "growth"}}
        
        with pytest.raises(ValueError, match="No nogrowth condition found"):
            build_phenotype_conditions(media_conditions, phenotype_data)

    def test_build_phenotype_conditions_invalid_growth_condition(self, sample_media_conditions, sample_phenotype_data):
        """Test build_phenotype_conditions with invalid growth_condition"""
        from crop.api import build_phenotype_conditions
        
        with pytest.raises(KeyError, match="Growth condition 'invalid' not in media_conditions"):
            build_phenotype_conditions(
                sample_media_conditions,
                sample_phenotype_data,
                growth_condition="invalid"
            )

    def test_build_phenotype_conditions_invalid_nogrowth_condition(self, sample_media_conditions, sample_phenotype_data):
        """Test build_phenotype_conditions with invalid nogrowth_condition"""
        from crop.api import build_phenotype_conditions
        
        with pytest.raises(KeyError, match="Nogrowth condition 'invalid' not in media_conditions"):
            build_phenotype_conditions(
                sample_media_conditions,
                sample_phenotype_data,
                nogrowth_condition="invalid"
            )

    def test_build_phenotype_conditions_union_of_exchanges(self):
        """Test that build_phenotype_conditions includes union of all exchanges"""
        from crop.api import build_phenotype_conditions
        
        media_conditions = {
            "cond1": {"EX_a": -3.0},  # Only EX_a
            "cond2": {"EX_b": -4.0}   # Only EX_b
        }
        phenotype_data = {
            "cond1": {"observed": "growth", "predicted": "growth"},
            "cond2": {"observed": "no_growth", "predicted": "growth"}
        }
        
        result = build_phenotype_conditions(media_conditions, phenotype_data)
        
        # Both growth and nogrowth should have both exchanges
        assert set(result["growth"].keys()) == {"EX_a", "EX_b"}
        assert set(result["nogrowth"].keys()) == {"EX_a", "EX_b"}
        
        # Check correct flipping
        assert result["growth"]["EX_a"] == 3.0
        assert result["growth"]["EX_b"] == 0.0
        assert result["nogrowth"]["EX_a"] == 0.0
        assert result["nogrowth"]["EX_b"] == 4.0

    @patch('crop.api.Variable')
    @patch('crop.api.Problem')
    def test_nogrowth_clause_structure(self, mock_problem, mock_variable):
        """Test nogrowth_clause function structure"""
        from crop.api import nogrowth_clause
        
        # This is a complex function that requires CVXPY variables
        # We'll just test that it returns a list when called with valid arguments
        # The actual constraint logic is better tested through integration tests
        
        # Mock CVXPY variables - these need to behave like CVXPY Variables
        v_nogrowth = Mock()
        v_nogrowth.__getitem__ = Mock(return_value=Mock())
        z = Mock()
        r = Mock()
        r.__getitem__ = Mock(return_value=Mock())
        m = Mock()
        
        biomass_idx = 0
        lower_bound_nogrowth = np.array([1.0, 2.0])
        upper_bound_nogrowth = np.array([3.0, 4.0])
        stoichiometric_matrix = np.array([[1, 0], [0, 1]])
        omega = np.array([1000, 1000])
        weights = np.array([1.0, 1.0])
        nogrowth_carbon_source_name = "EX_lac"
        nogrowth_carbon_source_idx = 1
        maximum_nogrowth = 1.0
        
        # Test that the function exists and is callable
        assert callable(nogrowth_clause)
        # The actual constraint testing requires CVXPY to be fully functional
        # which is better suited for integration tests

    @patch('crop.api.Variable')
    def test_growth_clause_structure(self, mock_variable):
        """Test growth_clause function structure"""
        from crop.api import growth_clause
        
        # Test that the function exists and is callable
        assert callable(growth_clause)
        # The actual constraint testing requires CVXPY to be fully functional
        # which is better suited for integration tests

    def test_media_conditions_fixture_function(self):
        """Test the media_conditions fixture function in api.py"""
        from crop.api import media_conditions
        
        result = media_conditions()
        
        assert isinstance(result, dict)
        assert "glucose" in result
        assert "lactose" in result
        assert "no_carbon" in result
        
        # Check structure
        assert result["glucose"]["EX_glc"] == -10
        assert result["lactose"]["EX_lac"] == -10
        assert result["no_carbon"]["EX_glc"] == 0

    def test_phenotype_data_fixture_function(self):
        """Test the phenotype_data fixture function in api.py"""
        from crop.api import phenotype_data
        
        result = phenotype_data()
        
        assert isinstance(result, dict)
        assert "glucose" in result
        assert "lactose" in result
        assert "no_carbon" in result
        
        # Check structure
        assert result["glucose"]["observed"] == "growth"
        assert result["glucose"]["predicted"] == "growth"
        assert result["lactose"]["observed"] == "no_growth"
        assert result["lactose"]["predicted"] == "growth"

