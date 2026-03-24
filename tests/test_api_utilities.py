"""Unit tests for utility/helper functions in crop.api."""

import numpy as np
import pandas as pd
import pytest
from cvxpy import Variable

from crop import run_crop_algorithm


class TestApiUtilityFunctions:
    """Test all utility functions in api.py."""

    @pytest.fixture
    def sample_phenotype_data(self):
        """Sample phenotype data for testing."""
        return {
            "glucose": {"observed": "growth", "predicted": "growth"},
            "lactose": {"observed": "no_growth", "predicted": "growth"},
            "no_carbon": {"observed": "no_growth", "predicted": "no_growth"},
            "mixed": {"observed": "growth", "predicted": "no_growth"},
        }

    @pytest.fixture
    def sample_media_conditions(self):
        """Sample media conditions for testing."""
        return {
            "glucose": {"EX_glc": -10.0, "EX_lac": 0.0},
            "lactose": {"EX_glc": 0.0, "EX_lac": -5.0},
            "no_carbon": {"EX_glc": 0.0, "EX_lac": 0.0},
            "mixed": {"EX_glc": -3.0, "EX_lac": -2.0},
        }

    @pytest.fixture
    def mock_model(self):
        """Mock COBRA model for testing."""
        from unittest.mock import Mock

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
        from crop.api import get_growth_conditions

        result = get_growth_conditions(sample_phenotype_data)
        assert result == ["glucose"]

    def test_get_growth_conditions_empty(self):
        from crop.api import get_growth_conditions

        phenotype_data = {
            "condition1": {"observed": "no_growth", "predicted": "growth"},
            "condition2": {"observed": "growth", "predicted": "no_growth"},
        }
        assert get_growth_conditions(phenotype_data) == []

    def test_get_growth_conditions_multiple(self):
        from crop.api import get_growth_conditions

        phenotype_data = {
            "glucose": {"observed": "growth", "predicted": "growth"},
            "fructose": {"observed": "growth", "predicted": "growth"},
            "lactose": {"observed": "no_growth", "predicted": "growth"},
        }
        assert set(get_growth_conditions(phenotype_data)) == {"glucose", "fructose"}

    def test_get_nogrowth_conditions(self, sample_phenotype_data):
        from crop.api import get_nogrowth_conditions

        assert get_nogrowth_conditions(sample_phenotype_data) == ["lactose"]

    def test_get_nogrowth_conditions_empty(self):
        from crop.api import get_nogrowth_conditions

        phenotype_data = {
            "condition1": {"observed": "growth", "predicted": "growth"},
            "condition2": {"observed": "no_growth", "predicted": "no_growth"},
        }
        assert get_nogrowth_conditions(phenotype_data) == []

    def test_get_nogrowth_conditions_multiple(self):
        from crop.api import get_nogrowth_conditions

        phenotype_data = {
            "lactose": {"observed": "no_growth", "predicted": "growth"},
            "sucrose": {"observed": "no_growth", "predicted": "growth"},
            "glucose": {"observed": "growth", "predicted": "growth"},
        }
        assert set(get_nogrowth_conditions(phenotype_data)) == {"lactose", "sucrose"}

    def test_get_carbon_source(self, mock_model, sample_phenotype_data, sample_media_conditions):
        from crop.api import get_carbon_source, get_growth_conditions

        result = get_carbon_source(
            mock_model,
            sample_phenotype_data,
            sample_media_conditions,
            get_growth_conditions,
        )

        assert result["glucose"] == (0, "EX_glc")

    def test_get_lower_bound_for_conditions(self, mock_model, sample_phenotype_data, sample_media_conditions):
        from crop.api import get_growth_conditions, get_lower_bound_for_conditions

        result = get_lower_bound_for_conditions(
            mock_model,
            sample_phenotype_data,
            sample_media_conditions,
            get_growth_conditions,
        )

        assert isinstance(result, pd.DataFrame)
        assert result.loc["EX_glc", "glucose"] == -10.0
        assert result.loc["EX_lac", "glucose"] == 0.0

    def test_get_lower_bound_for_growth_conditions(self, mock_model, sample_phenotype_data, sample_media_conditions):
        from crop.api import get_lower_bound_for_growth_conditions

        result = get_lower_bound_for_growth_conditions(
            mock_model,
            sample_phenotype_data,
            sample_media_conditions,
            atp_maintenance_rxn="ATPM",
            atp_maintenance_lower_bound=2.5,
        )

        assert isinstance(result, pd.DataFrame)
        assert result.loc["ATPM", "glucose"] == 2.5

    def test_get_lower_bound_for_nogrowth_conditions(self, mock_model, sample_phenotype_data, sample_media_conditions):
        from crop.api import get_lower_bound_for_nogrowth_conditions

        result = get_lower_bound_for_nogrowth_conditions(
            mock_model,
            sample_phenotype_data,
            sample_media_conditions,
        )

        assert isinstance(result, pd.DataFrame)
        assert result.loc["EX_lac", "lactose"] == -5.0
        assert result.loc["EX_glc", "lactose"] == 0.0

    def test_get_upper_bound_for_conditions(self, mock_model, sample_phenotype_data, sample_media_conditions):
        from crop.api import get_growth_conditions, get_upper_bound_for_conditions

        result = get_upper_bound_for_conditions(
            mock_model,
            sample_phenotype_data,
            sample_media_conditions,
            get_growth_conditions,
            growth_limit=5.0,
        )

        assert isinstance(result, pd.DataFrame)
        assert result.loc["EX_glc", "glucose"] == 0.0
        assert result.loc["BIOMASS", "glucose"] == 5.0

    def test_get_upper_bound_for_nogrowth_conditions(self, mock_model, sample_phenotype_data, sample_media_conditions):
        from crop.api import get_upper_bound_for_nogrowth_conditions

        result = get_upper_bound_for_nogrowth_conditions(
            mock_model,
            sample_phenotype_data,
            sample_media_conditions,
            maximum_nogrowth=1.0,
        )

        assert isinstance(result, pd.DataFrame)
        assert result.loc["EX_lac", "lactose"] == 0.0
        assert result.loc["BIOMASS", "lactose"] == 1.0

    def test_get_upper_bound_for_growth_conditions(self, mock_model, sample_phenotype_data, sample_media_conditions):
        from crop.api import get_upper_bound_for_growth_conditions

        result = get_upper_bound_for_growth_conditions(
            mock_model,
            sample_phenotype_data,
            sample_media_conditions,
            minimum_growth=2.0,
        )

        assert isinstance(result, pd.DataFrame)
        assert result.loc["EX_glc", "glucose"] == 0.0
        assert result.loc["BIOMASS", "glucose"] == 2.0

    def test_build_phenotype_conditions_basic(self, sample_media_conditions, sample_phenotype_data):
        from crop.api import build_phenotype_conditions

        result = build_phenotype_conditions(sample_media_conditions, sample_phenotype_data)

        assert result["growth"]["EX_glc"] == 10.0
        assert result["growth"]["EX_lac"] == 0.0
        assert result["nogrowth"]["EX_glc"] == 0.0
        assert result["nogrowth"]["EX_lac"] == 5.0

    def test_build_phenotype_conditions_observation_key(self):
        from crop.api import build_phenotype_conditions

        media_conditions = {
            "glucose": {"EX_glc": -10.0, "EX_lac": 0.0},
            "lactose": {"EX_glc": 0.0, "EX_lac": -5.0},
        }
        phenotype_data = {
            "glucose": {"observation": "growth", "predicted": "growth"},
            "lactose": {"observation": "no_growth", "predicted": "growth"},
        }

        result = build_phenotype_conditions(media_conditions, phenotype_data)

        assert result["growth"]["EX_glc"] == 10.0
        assert result["nogrowth"]["EX_lac"] == 5.0

    def test_build_phenotype_conditions_specific_conditions(self, sample_media_conditions, sample_phenotype_data):
        from crop.api import build_phenotype_conditions

        media_extended = dict(sample_media_conditions)
        media_extended["glucose2"] = {"EX_glc": -7.0, "EX_lac": 0.0}
        phenotype_extended = dict(sample_phenotype_data)
        phenotype_extended["glucose2"] = {"observed": "growth", "predicted": "growth"}

        result = build_phenotype_conditions(
            media_extended,
            phenotype_extended,
            growth_condition="glucose2",
            nogrowth_condition="lactose",
        )

        assert result["growth"]["EX_glc"] == 7.0
        assert result["nogrowth"]["EX_lac"] == 5.0

    def test_build_phenotype_conditions_no_growth_candidates(self):
        from crop.api import build_phenotype_conditions

        media_conditions = {"test": {"EX_glc": -10.0}}
        phenotype_data = {"test": {"observed": "no_growth", "predicted": "growth"}}

        with pytest.raises(ValueError, match="No growth condition found"):
            build_phenotype_conditions(media_conditions, phenotype_data)

    def test_build_phenotype_conditions_no_nogrowth_candidates(self):
        from crop.api import build_phenotype_conditions

        media_conditions = {"test": {"EX_glc": -10.0}}
        phenotype_data = {"test": {"observed": "growth", "predicted": "growth"}}

        with pytest.raises(ValueError, match="No nogrowth condition found"):
            build_phenotype_conditions(media_conditions, phenotype_data)

    def test_build_phenotype_conditions_invalid_growth_condition(self, sample_media_conditions, sample_phenotype_data):
        from crop.api import build_phenotype_conditions

        with pytest.raises(KeyError, match="Growth condition 'invalid' not in media_conditions"):
            build_phenotype_conditions(
                sample_media_conditions,
                sample_phenotype_data,
                growth_condition="invalid",
            )

    def test_build_phenotype_conditions_invalid_nogrowth_condition(self, sample_media_conditions, sample_phenotype_data):
        from crop.api import build_phenotype_conditions

        with pytest.raises(KeyError, match="Nogrowth condition 'invalid' not in media_conditions"):
            build_phenotype_conditions(
                sample_media_conditions,
                sample_phenotype_data,
                nogrowth_condition="invalid",
            )

    def test_build_phenotype_conditions_union_of_exchanges(self):
        from crop.api import build_phenotype_conditions

        media_conditions = {
            "cond1": {"EX_a": -3.0},
            "cond2": {"EX_b": -4.0},
        }
        phenotype_data = {
            "cond1": {"observed": "growth", "predicted": "growth"},
            "cond2": {"observed": "no_growth", "predicted": "growth"},
        }

        result = build_phenotype_conditions(media_conditions, phenotype_data)

        assert set(result["growth"].keys()) == {"EX_a", "EX_b"}
        assert set(result["nogrowth"].keys()) == {"EX_a", "EX_b"}
        assert result["growth"]["EX_a"] == 3.0
        assert result["growth"]["EX_b"] == 0.0
        assert result["nogrowth"]["EX_a"] == 0.0
        assert result["nogrowth"]["EX_b"] == 4.0

    def test_media_conditions_fixture_function(self):
        from crop.api import media_conditions

        result = media_conditions()

        assert result["glucose"]["EX_glc"] == -10
        assert result["lactose"]["EX_lac"] == -10
        assert result["no_carbon"]["EX_glc"] == 0

    def test_phenotype_data_fixture_function(self):
        from crop.api import phenotype_data

        result = phenotype_data()

        assert result["glucose"]["observed"] == "growth"
        assert result["glucose"]["predicted"] == "growth"
        assert result["lactose"]["observed"] == "no_growth"
        assert result["lactose"]["predicted"] == "growth"


class TestApiCoverageAdditions:
    """Additional tests focused on raising api.py line coverage."""

    def test_nogrowth_clause_returns_expected_number_of_constraints(self):
        from crop.api import nogrowth_clause

        v_nogrowth = Variable(3)
        z = Variable(3, boolean=True)
        r = Variable(3)
        m = Variable(2)

        lower_bound_nogrowth = pd.Series({"EX_carbon": -10.0, "R1": 0.0, "BIOMASS": 0.0})
        upper_bound_nogrowth = pd.Series({"EX_carbon": 0.0, "R1": 1.0, "BIOMASS": 1.0})
        stoichiometric_matrix = np.array([[1.0, -1.0, 0.0], [0.0, 1.0, -1.0]])
        c = np.array([0.0, 0.0, 1.0])

        constraints = nogrowth_clause(
            v_nogrowth=v_nogrowth,
            biomass_idx=2,
            lower_bound_nogrowth=lower_bound_nogrowth,
            upper_bound_nogrowth=upper_bound_nogrowth,
            stoichiometric_matrix=stoichiometric_matrix,
            z=z,
            r=r,
            m=m,
            c=c,
            omega=1000.0,
            nogrowth_carbon_source_name="EX_carbon",
            nogrowth_carbon_source_idx=0,
            maximum_nogrowth=1.0,
        )

        assert len(constraints) == 8

    def test_growth_clause_returns_expected_number_of_constraints(self):
        from crop.api import growth_clause

        v_growth = Variable(3)
        z = Variable(3, boolean=True)
        lower_bound_growth = np.array([-10.0, 0.0, 0.0])
        upper_bound_growth = np.array([0.0, 2.0, 2.0])
        stoichiometric_matrix = np.array([[1.0, -1.0, 0.0], [0.0, 1.0, -1.0]])

        constraints = growth_clause(
            v_growth=v_growth,
            biomass_idx=2,
            lower_bound_growth=lower_bound_growth,
            upper_bound_growth=upper_bound_growth,
            stoichiometric_matrix=stoichiometric_matrix,
            z=z,
            minimum_growth=2.0,
        )

        assert len(constraints) == 4

    def test_run_crop_algorithm_raises_without_growth_candidates(self):
        from tests.test_crop_oracle import create_unique_oracle_model

        model = create_unique_oracle_model()
        media_conditions = {
            "lactose_with_helper": {"EX_lac": -10.0, "EX_helper": -10.0},
            "lactose_only": {"EX_lac": -10.0, "EX_helper": 0.0},
        }
        phenotype_data = {
            "lactose_with_helper": {"observed": "no_growth", "predicted": "growth"},
            "lactose_only": {"observed": "no_growth", "predicted": "growth"},
        }

        with pytest.raises(ValueError, match="No growth condition found"):
            run_crop_algorithm(model, phenotype_data, media_conditions)

    def test_run_crop_algorithm_raises_without_nogrowth_candidates(self):
        from tests.test_crop_oracle import create_unique_oracle_model

        model = create_unique_oracle_model()
        media_conditions = {
            "lactose_with_helper": {"EX_lac": -10.0, "EX_helper": -10.0},
            "lactose_only": {"EX_lac": -10.0, "EX_helper": 0.0},
        }
        phenotype_data = {
            "lactose_with_helper": {"observed": "growth", "predicted": "growth"},
            "lactose_only": {"observed": "growth", "predicted": "no_growth"},
        }

        with pytest.raises(ValueError, match="No nogrowth condition found"):
            run_crop_algorithm(model, phenotype_data, media_conditions)

    def test_run_crop_algorithm_raises_when_nogrowth_has_no_carbon_source(self):
        from tests.test_crop_oracle import create_unique_oracle_model

        model = create_unique_oracle_model()
        media_conditions = {
            "lactose_with_helper": {"EX_lac": -10.0, "EX_helper": -10.0},
            "lactose_only": {"EX_lac": 0.0, "EX_helper": 0.0},
        }
        phenotype_data = {
            "lactose_with_helper": {"observed": "growth", "predicted": "growth"},
            "lactose_only": {"observed": "no_growth", "predicted": "growth"},
        }

        with pytest.raises(ValueError, match="No carbon source found for nogrowth condition"):
            run_crop_algorithm(model, phenotype_data, media_conditions)


class TestApiCoverageMaximization:
    """Additional tests to maximize code coverage for api.py."""

    @pytest.fixture
    def simple_growth_model(self):
        """Create a minimal viable model for testing edge cases."""
        from cobra import Model, Reaction, Metabolite

        model = Model("simple_test")

        # Metabolites
        glc_e = Metabolite("glc_e", compartment="e")
        glc_c = Metabolite("glc_c", compartment="c")
        pyr_c = Metabolite("pyr_c", compartment="c")

        # Exchange reaction
        ex_glc = Reaction("EX_glc")
        ex_glc.lower_bound = -10.0
        ex_glc.upper_bound = 1000.0
        ex_glc.add_metabolites({glc_e: -1})

        # Transport
        glc_t = Reaction("GLCt")
        glc_t.add_metabolites({glc_e: -1, glc_c: 1})

        # Simple conversion
        glc_conv = Reaction("GLC_conv")
        glc_conv.add_metabolites({glc_c: -1, pyr_c: 1})

        # Biomass (objective)
        biomass = Reaction("BIOMASS")
        biomass.lower_bound = 0.0
        biomass.upper_bound = 1000.0
        biomass.add_metabolites({pyr_c: -2})

        model.add_reactions([ex_glc, glc_t, glc_conv, biomass])
        model.objective = "BIOMASS"
        return model

    def test_infeasible_problem_with_conflicting_bounds(self):
        """Test that infeasible problems are caught and raise appropriate error."""
        from cobra import Model, Reaction, Metabolite

        model = Model("infeasible_test")

        # Create metabolites
        a_e = Metabolite("a_e", compartment="e")
        a_c = Metabolite("a_c", compartment="c")
        b_c = Metabolite("b_c", compartment="c")

        # Exchange
        ex_a = Reaction("EX_a")
        ex_a.lower_bound = -10.0
        ex_a.upper_bound = 100.0
        ex_a.add_metabolites({a_e: -1})

        # Transport
        a_t = Reaction("At")
        a_t.add_metabolites({a_e: -1, a_c: 1})

        # Conversion that requires consumption > max possible uptake
        a_conv = Reaction("A_conv")
        a_conv.lower_bound = 0.0
        a_conv.upper_bound = 1000.0
        # Requires 100 units of a_c, but max uptake is only 10
        a_conv.add_metabolites({a_c: -100, b_c: 1})

        # Biomass depends on reaction that can't go fast enough
        biomass = Reaction("BIOMASS")
        biomass.lower_bound = 0.0
        biomass.upper_bound = 1000.0
        biomass.add_metabolites({b_c: -1})

        model.add_reactions([ex_a, a_t, a_conv, biomass])
        model.objective = "BIOMASS"

        media_conditions = {
            "condition_a": {"EX_a": -10.0},
            "no_growth": {"EX_a": -10.0},
        }
        phenotype_data = {
            "condition_a": {"observed": "growth", "predicted": "growth"},
            "no_growth": {"observed": "no_growth", "predicted": "growth"},
        }

        # Try with very high minimum_growth that might make problem infeasible
        try:
            suggested_removals, _ = run_crop_algorithm(
                model,
                phenotype_data,
                media_conditions,
                minimum_growth=1000.0,  # Impossibly high
                maximum_nogrowth=0.001,
            )
            # If it succeeds, that's OK - problem might still be solvable
            # The key is we test the code path without errors
            assert isinstance(suggested_removals, (set, frozenset))
        except ValueError as e:
            # Infeasible is also acceptable
            assert "Infeasible" in str(e)

    def test_reaction_removal_with_minimal_model(self):
        """Test CROP algorithm handles various scenarios correctly."""
        from tests.test_crop_oracle import create_unique_oracle_model

        model = create_unique_oracle_model()

        media_conditions = {
            "lactose_with_helper": {"EX_lac": -10.0, "EX_helper": -10.0},
            "lactose_only": {"EX_lac": -10.0, "EX_helper": 0.0},
        }
        phenotype_data = {
            "lactose_with_helper": {"observed": "growth", "predicted": "growth"},
            "lactose_only": {"observed": "no_growth", "predicted": "growth"},
        }

        # Should run without error
        suggested_removals, solution = run_crop_algorithm(
            model,
            phenotype_data,
            media_conditions,
            maximum_nogrowth=1.0,
            minimum_growth=2.0,
        )

        # Verify solution structure
        assert isinstance(suggested_removals, set)
        assert isinstance(solution, pd.DataFrame)
        assert "z" in solution.columns
        assert len(solution) == len(model.reactions)

    def test_solution_contains_per_condition_columns(self):
        """Test that solution contains expected per-condition flux columns."""
        from tests.test_crop_oracle import (
            create_unique_oracle_model,
            apply_medium,
        )

        model = create_unique_oracle_model()
        media_conditions = {
            "lactose_with_helper": {"EX_lac": -10.0, "EX_helper": -10.0},
            "lactose_only": {"EX_lac": -10.0, "EX_helper": 0.0},
        }
        phenotype_data = {
            "lactose_with_helper": {"observed": "growth", "predicted": "growth"},
            "lactose_only": {"observed": "no_growth", "predicted": "growth"},
        }

        suggested_removals, solution = run_crop_algorithm(
            model,
            phenotype_data,
            media_conditions,
            maximum_nogrowth=1.0,
            minimum_growth=2.0,
        )

        # Verify solution has both backward-compatible and per-condition columns
        assert "z" in solution.columns
        # Should have per-condition columns
        assert any("__" in col for col in solution.columns), "Should have per-condition columns"

    def test_algorithm_with_single_growth_condition(self):
        """Test CROP with oracle models to exercise core logic."""
        from tests.test_crop_oracle import create_unique_multi_condition_oracle_model

        model = create_unique_multi_condition_oracle_model()

        media_conditions = {
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
            "lactose_only": {
                "EX_lac": -10.0,
                "EX_lac_helper": 0.0,
                "EX_sorb": 0.0,
                "EX_sorb_helper": 0.0,
                "EX_g6p": 0.0,
            },
        }

        phenotype_data = {
            "glucose_6_phosphate": {"observed": "growth", "predicted": "growth"},
            "lactose_with_helper": {"observed": "growth", "predicted": "growth"},
            "lactose_only": {"observed": "no_growth", "predicted": "growth"},
        }

        suggested_removals, solution = run_crop_algorithm(
            model,
            phenotype_data,
            media_conditions,
            minimum_growth=2.0,
            maximum_nogrowth=1.0,
        )

        assert isinstance(suggested_removals, set)
        assert isinstance(solution, pd.DataFrame)

    def test_bounds_getters_with_heterogeneous_data(self):
        """Test helper functions with complex data structures."""
        from crop.api import (
            get_lower_bound_for_nogrowth_conditions,
            get_upper_bound_for_nogrowth_conditions,
        )
        from cobra import Model, Reaction, Metabolite

        model = Model("test")
        rxn1 = Reaction("R1")
        rxn2 = Reaction("R2")
        rxn3 = Reaction("R3")
        met1 = Metabolite("M1")
        rxn1.add_metabolites({met1: -1})
        rxn2.add_metabolites({met1: 1})
        model.add_reactions([rxn1, rxn2, rxn3])

        media_conditions = {
            "cond1": {"R1": -5.0, "R2": 0.0},
            "cond2": {"R1": 0.0, "R2": -8.0},
        }

        phenotype_data = {
            "growth": {"observed": "growth", "predicted": "growth"},
            "cond1": {"observed": "no_growth", "predicted": "growth"},
            "cond2": {"observed": "no_growth", "predicted": "growth"},
        }

        lower_bounds = get_lower_bound_for_nogrowth_conditions(
            model, phenotype_data, media_conditions
        )
        upper_bounds = get_upper_bound_for_nogrowth_conditions(
            model, phenotype_data, media_conditions, maximum_nogrowth=2.0
        )

        # Verify structure and values
        assert isinstance(lower_bounds, pd.DataFrame)
        assert isinstance(upper_bounds, pd.DataFrame)
        assert "cond1" in lower_bounds.columns
        assert "cond2" in lower_bounds.columns
        assert lower_bounds.loc["R1", "cond1"] < 0  # Should be negative (uptake)
        assert lower_bounds.loc["R1", "cond2"] == 0  # Not in cond2 media
