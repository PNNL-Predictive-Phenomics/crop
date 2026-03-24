"""Oracle-style tests for the CROP reaction-removal algorithm."""

from itertools import combinations

from cobra import Metabolite, Model, Reaction

from crop import run_crop_algorithm


def create_unique_oracle_model():
    """Build a toy model with a single minimal correction set.

    The model has one correct growth condition and one incorrect no-growth condition.
    Lactose plus a helper metabolite should support growth through the legitimate pathway,
    while lactose alone incorrectly grows through a bypass reaction. Deleting the bypass is
    the only single-reaction removal that fixes the no-growth condition without breaking the
    growth condition.
    """

    model = Model("crop_oracle_model")

    lac_e = Metabolite("lac_e", compartment="e")
    lac_c = Metabolite("lac_c", compartment="c")
    helper_e = Metabolite("helper_e", compartment="e")
    helper_c = Metabolite("helper_c", compartment="c")
    pyr_c = Metabolite("pyr_c", compartment="c")

    ex_lac = Reaction("EX_lac")
    ex_lac.lower_bound = -10.0
    ex_lac.upper_bound = 100.0
    ex_lac.add_metabolites({lac_e: -1})

    ex_helper = Reaction("EX_helper")
    ex_helper.lower_bound = -10.0
    ex_helper.upper_bound = 100.0
    ex_helper.add_metabolites({helper_e: -1})

    lac_transport = Reaction("LACt")
    lac_transport.add_metabolites({lac_e: -1, lac_c: 1})

    helper_transport = Reaction("HELPt")
    helper_transport.add_metabolites({helper_e: -1, helper_c: 1})

    legitimate_lactose_use = Reaction("LAC_with_helper")
    legitimate_lactose_use.add_metabolites({lac_c: -1, helper_c: -1, pyr_c: 3})

    bypass_lactose_use = Reaction("LAC_bypass")
    bypass_lactose_use.add_metabolites({lac_c: -1, pyr_c: 3})

    biomass = Reaction("BIOMASS")
    biomass.lower_bound = 0.0
    biomass.upper_bound = 100.0
    biomass.add_metabolites({pyr_c: -1})

    model.add_reactions(
        [
            ex_lac,
            ex_helper,
            lac_transport,
            helper_transport,
            legitimate_lactose_use,
            bypass_lactose_use,
            biomass,
        ]
    )
    model.objective = "BIOMASS"
    return model


def apply_medium(model, medium):
    """Return a copy of model with exchange lower bounds set from medium."""

    model_copy = model.copy()
    for reaction_id, bound in medium.items():
        if reaction_id in model_copy.reactions:
            model_copy.reactions.get_by_id(reaction_id).lower_bound = bound
    return model_copy


def brute_force_minimal_removals(
    model,
    media_conditions,
    minimum_growth=2.0,
    maximum_nogrowth=1.0,
):
    """Enumerate all reaction-deletion subsets and return the minimal valid ones."""

    candidate_reactions = [reaction.id for reaction in model.reactions]
    growth_medium = media_conditions["lactose_with_helper"]
    nogrowth_medium = media_conditions["lactose_only"]

    for subset_size in range(len(candidate_reactions) + 1):
        valid_subsets = []
        for subset in combinations(candidate_reactions, subset_size):
            model_candidate = model.copy()
            for reaction_id in subset:
                if reaction_id in model_candidate.reactions:
                    model_candidate.reactions.get_by_id(reaction_id).remove_from_model()

            growth_flux = apply_medium(model_candidate, growth_medium).optimize().objective_value
            nogrowth_flux = apply_medium(model_candidate, nogrowth_medium).optimize().objective_value

            if growth_flux >= minimum_growth and nogrowth_flux <= maximum_nogrowth:
                valid_subsets.append(frozenset(subset))

        if valid_subsets:
            return set(valid_subsets)

    return set()


def create_unique_multi_condition_oracle_model():
    """Build a toy model whose only minimal multi-condition fix removes two bypasses."""

    model = Model("crop_multi_condition_oracle_model")

    lac_e = Metabolite("lac_e", compartment="e")
    lac_c = Metabolite("lac_c", compartment="c")
    lac_helper_e = Metabolite("lac_helper_e", compartment="e")
    lac_helper_c = Metabolite("lac_helper_c", compartment="c")
    sorb_e = Metabolite("sorb_e", compartment="e")
    sorb_c = Metabolite("sorb_c", compartment="c")
    sorb_helper_e = Metabolite("sorb_helper_e", compartment="e")
    sorb_helper_c = Metabolite("sorb_helper_c", compartment="c")
    g6p_e = Metabolite("g6p_e", compartment="e")
    g6p_c = Metabolite("g6p_c", compartment="c")
    pyr_c = Metabolite("pyr_c", compartment="c")

    ex_lac = Reaction("EX_lac")
    ex_lac.lower_bound = -10.0
    ex_lac.upper_bound = 100.0
    ex_lac.add_metabolites({lac_e: -1})

    ex_lac_helper = Reaction("EX_lac_helper")
    ex_lac_helper.lower_bound = -10.0
    ex_lac_helper.upper_bound = 100.0
    ex_lac_helper.add_metabolites({lac_helper_e: -1})

    ex_sorb = Reaction("EX_sorb")
    ex_sorb.lower_bound = -10.0
    ex_sorb.upper_bound = 100.0
    ex_sorb.add_metabolites({sorb_e: -1})

    ex_sorb_helper = Reaction("EX_sorb_helper")
    ex_sorb_helper.lower_bound = -10.0
    ex_sorb_helper.upper_bound = 100.0
    ex_sorb_helper.add_metabolites({sorb_helper_e: -1})

    ex_g6p = Reaction("EX_g6p")
    ex_g6p.lower_bound = -10.0
    ex_g6p.upper_bound = 100.0
    ex_g6p.add_metabolites({g6p_e: -1})

    lac_transport = Reaction("LACt")
    lac_transport.add_metabolites({lac_e: -1, lac_c: 1})

    lac_helper_transport = Reaction("LAC_HELPt")
    lac_helper_transport.add_metabolites({lac_helper_e: -1, lac_helper_c: 1})

    sorb_transport = Reaction("SORBt")
    sorb_transport.add_metabolites({sorb_e: -1, sorb_c: 1})

    sorb_helper_transport = Reaction("SORB_HELPt")
    sorb_helper_transport.add_metabolites({sorb_helper_e: -1, sorb_helper_c: 1})

    g6p_transport = Reaction("G6Pt")
    g6p_transport.add_metabolites({g6p_e: -1, g6p_c: 1})

    legitimate_lactose_use = Reaction("LAC_with_helper")
    legitimate_lactose_use.add_metabolites({lac_c: -1, lac_helper_c: -1, g6p_c: 1})

    legitimate_sorbose_use = Reaction("SORB_with_helper")
    legitimate_sorbose_use.add_metabolites({sorb_c: -1, sorb_helper_c: -1, g6p_c: 1})

    lactose_bypass = Reaction("LAC_bypass")
    lactose_bypass.add_metabolites({lac_c: -1, g6p_c: 1})

    sorbose_bypass = Reaction("SORB_bypass")
    sorbose_bypass.add_metabolites({sorb_c: -1, g6p_c: 1})

    glycolysis = Reaction("GLYC")
    glycolysis.add_metabolites({g6p_c: -1, pyr_c: 2})

    biomass = Reaction("BIOMASS")
    biomass.lower_bound = 0.0
    biomass.upper_bound = 100.0
    biomass.add_metabolites({pyr_c: -1})

    model.add_reactions(
        [
            ex_lac,
            ex_lac_helper,
            ex_sorb,
            ex_sorb_helper,
            ex_g6p,
            lac_transport,
            lac_helper_transport,
            sorb_transport,
            sorb_helper_transport,
            g6p_transport,
            legitimate_lactose_use,
            legitimate_sorbose_use,
            lactose_bypass,
            sorbose_bypass,
            glycolysis,
            biomass,
        ]
    )
    model.objective = "BIOMASS"
    return model


def multi_condition_oracle_media_conditions():
    """Return shared media conditions for the multi-condition oracle model."""

    return {
        "glucose_6_phosphate": {
            "EX_lac": 0.0,
            "EX_lac_helper": 0.0,
            "EX_sorb": 0.0,
            "EX_sorb_helper": 0.0,
            "EX_g6p": -10.0,
        },
        "lactose_with_helper": {
            "EX_lac": -10.0,
            "EX_lac_helper": -10.0,
            "EX_sorb": 0.0,
            "EX_sorb_helper": 0.0,
            "EX_g6p": 0.0,
        },
        "sorbose_with_helper": {
            "EX_lac": 0.0,
            "EX_lac_helper": 0.0,
            "EX_sorb": -10.0,
            "EX_sorb_helper": -10.0,
            "EX_g6p": 0.0,
        },
        "lactose_only": {
            "EX_lac": -10.0,
            "EX_lac_helper": 0.0,
            "EX_sorb": 0.0,
            "EX_sorb_helper": 0.0,
            "EX_g6p": 0.0,
        },
        "sorbose_only": {
            "EX_lac": 0.0,
            "EX_lac_helper": 0.0,
            "EX_sorb": -10.0,
            "EX_sorb_helper": 0.0,
            "EX_g6p": 0.0,
        },
        "no_carbon": {
            "EX_lac": 0.0,
            "EX_lac_helper": 0.0,
            "EX_sorb": 0.0,
            "EX_sorb_helper": 0.0,
            "EX_g6p": 0.0,
        },
    }


def brute_force_multi_condition_minimal_removals(
    model,
    media_conditions,
    expected_growth,
    growth_threshold=2.0,
    nogrowth_threshold=1.0,
):
    """Enumerate deletion sets and return the minimal ones satisfying all conditions."""

    candidate_reactions = [reaction.id for reaction in model.reactions]

    for subset_size in range(len(candidate_reactions) + 1):
        valid_subsets = []
        for subset in combinations(candidate_reactions, subset_size):
            model_candidate = model.copy()
            for reaction_id in subset:
                if reaction_id in model_candidate.reactions:
                    model_candidate.reactions.get_by_id(reaction_id).remove_from_model()

            matches_all_conditions = True
            for condition_name, should_grow in expected_growth.items():
                objective_value = apply_medium(
                    model_candidate, media_conditions[condition_name]
                ).optimize().objective_value
                if should_grow and objective_value < growth_threshold:
                    matches_all_conditions = False
                    break
                if not should_grow and objective_value > nogrowth_threshold:
                    matches_all_conditions = False
                    break

            if matches_all_conditions:
                valid_subsets.append(frozenset(subset))

        if valid_subsets:
            return set(valid_subsets)

    return set()


def test_oracle_model_has_unique_minimal_fix():
    """The toy model should have exactly one minimal correction set."""

    model = create_unique_oracle_model()
    media_conditions = {
        "lactose_with_helper": {"EX_lac": -10.0, "EX_helper": -10.0},
        "lactose_only": {"EX_lac": -10.0, "EX_helper": 0.0},
    }

    minimal_removals = brute_force_minimal_removals(model, media_conditions)

    assert minimal_removals == {frozenset({"LAC_bypass"})}


def test_crop_algorithm_matches_bruteforce_oracle():
    """CROP should match the unique minimal deletion set on the oracle model."""

    model = create_unique_oracle_model()
    media_conditions = {
        "lactose_with_helper": {"EX_lac": -10.0, "EX_helper": -10.0},
        "lactose_only": {"EX_lac": -10.0, "EX_helper": 0.0},
    }
    phenotype_data = {
        "lactose_with_helper": {"observed": "growth", "predicted": "growth"},
        "lactose_only": {"observed": "no_growth", "predicted": "growth"},
    }

    oracle_minimal_removals = brute_force_minimal_removals(model, media_conditions)
    suggested_removals, _ = run_crop_algorithm(
        model,
        phenotype_data,
        media_conditions,
        maximum_nogrowth=1.0,
        minimum_growth=2.0,
    )

    assert frozenset(suggested_removals) in oracle_minimal_removals

    corrected_model = model.copy()
    for reaction_id in suggested_removals:
        corrected_model.reactions.get_by_id(reaction_id).remove_from_model()

    growth_solution = apply_medium(
        corrected_model, media_conditions["lactose_with_helper"]
    ).optimize()
    nogrowth_solution = apply_medium(
        corrected_model, media_conditions["lactose_only"]
    ).optimize()

    assert growth_solution.objective_value >= 2.0
    assert nogrowth_solution.objective_value <= 1.0


def test_multi_condition_oracle_model_has_unique_minimal_fix():
    """The multi-condition oracle should require removing both bypass reactions."""

    model = create_unique_multi_condition_oracle_model()
    media_conditions = multi_condition_oracle_media_conditions()
    expected_growth = {
        "glucose_6_phosphate": True,
        "lactose_with_helper": True,
        "sorbose_with_helper": True,
        "lactose_only": False,
        "sorbose_only": False,
        "no_carbon": False,
    }

    minimal_removals = brute_force_multi_condition_minimal_removals(
        model, media_conditions, expected_growth
    )

    assert minimal_removals == {frozenset({"LAC_bypass", "SORB_bypass"})}


def test_crop_algorithm_matches_multicondition_oracle():
    """CROP should match the multi-condition oracle on the canonical mixed set."""

    model = create_unique_multi_condition_oracle_model()
    media_conditions = multi_condition_oracle_media_conditions()
    phenotype_data = {
        "glucose_6_phosphate": {"observed": "growth", "predicted": "growth"},
        "lactose_with_helper": {"observed": "growth", "predicted": "growth"},
        "sorbose_with_helper": {"observed": "growth", "predicted": "growth"},
        "lactose_only": {"observed": "no_growth", "predicted": "growth"},
        "sorbose_only": {"observed": "no_growth", "predicted": "growth"},
        "no_carbon": {"observed": "no_growth", "predicted": "no_growth"},
    }
    expected_growth = {
        "glucose_6_phosphate": True,
        "lactose_with_helper": True,
        "sorbose_with_helper": True,
        "lactose_only": False,
        "sorbose_only": False,
        "no_carbon": False,
    }

    oracle_minimal_removals = brute_force_multi_condition_minimal_removals(
        model, media_conditions, expected_growth
    )
    suggested_removals, _ = run_crop_algorithm(
        model,
        phenotype_data,
        media_conditions,
        maximum_nogrowth=1.0,
        minimum_growth=2.0,
    )

    assert frozenset(suggested_removals) in oracle_minimal_removals


def test_example_multigrowth_multinogrowth_returns_two_removals():
    """Example: multiple growth and multiple no-growth conditions."""

    model = create_unique_multi_condition_oracle_model()
    media_conditions = multi_condition_oracle_media_conditions()
    phenotype_data = {
        "glucose_6_phosphate": {"observed": "growth", "predicted": "growth"},
        "lactose_with_helper": {"observed": "growth", "predicted": "growth"},
        "sorbose_with_helper": {"observed": "growth", "predicted": "growth"},
        "lactose_only": {"observed": "no_growth", "predicted": "growth"},
        "sorbose_only": {"observed": "no_growth", "predicted": "growth"},
        "no_carbon": {"observed": "no_growth", "predicted": "no_growth"},
    }

    suggested_removals, _ = run_crop_algorithm(model, phenotype_data, media_conditions)

    assert suggested_removals == {"LAC_bypass", "SORB_bypass"}


def test_example_multigrowth_single_nogrowth_returns_single_removal():
    """Example: multiple growth conditions and one targeted no-growth condition."""

    model = create_unique_multi_condition_oracle_model()
    media_conditions = multi_condition_oracle_media_conditions()
    phenotype_data = {
        "glucose_6_phosphate": {"observed": "growth", "predicted": "growth"},
        "lactose_with_helper": {"observed": "growth", "predicted": "growth"},
        "sorbose_with_helper": {"observed": "growth", "predicted": "growth"},
        "lactose_only": {"observed": "no_growth", "predicted": "growth"},
        "sorbose_only": {"observed": "no_growth", "predicted": "no_growth"},
        "no_carbon": {"observed": "no_growth", "predicted": "no_growth"},
    }

    suggested_removals, _ = run_crop_algorithm(model, phenotype_data, media_conditions)

    assert suggested_removals == {"LAC_bypass"}


def test_example_mixed_outcomes_ignores_non_target_categories():
    """Example: mixed outcomes including false negatives and true no-growth cases."""

    model = create_unique_multi_condition_oracle_model()
    media_conditions = multi_condition_oracle_media_conditions()
    phenotype_data = {
        "glucose_6_phosphate": {"observed": "growth", "predicted": "growth"},
        "lactose_with_helper": {"observed": "growth", "predicted": "growth"},
        "sorbose_with_helper": {"observed": "growth", "predicted": "no_growth"},
        "lactose_only": {"observed": "no_growth", "predicted": "growth"},
        "sorbose_only": {"observed": "no_growth", "predicted": "growth"},
        "no_carbon": {"observed": "no_growth", "predicted": "no_growth"},
    }
    expected_growth = {
        "glucose_6_phosphate": True,
        "lactose_with_helper": True,
        "lactose_only": False,
        "sorbose_only": False,
        "no_carbon": False,
    }

    oracle_minimal_removals = brute_force_multi_condition_minimal_removals(
        model, media_conditions, expected_growth
    )

    suggested_removals, _ = run_crop_algorithm(model, phenotype_data, media_conditions)

    assert frozenset(suggested_removals) in oracle_minimal_removals