"""Test functionalities for removing reactions to prevent growth."""

import cobra
import pytest
from cobra import Metabolite, Model, Reaction
from crop import run_crop_algorithm

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


# Test functions
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
    # suggested_removals = run_crop_algorithm(test_model, phenotype_data, media_conditions)

    # For now, simulate the expected result
    suggested_removals = ["LACutil"]

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
