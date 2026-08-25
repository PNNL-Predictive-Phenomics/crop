"""Regression tests for the single- versus multi-phenotype comparison."""

from scripts.compare_growmatch_crop import (
    build_results,
    compare_strategies,
    create_cross_coupled_model,
)


def test_joint_reconciliation_avoids_independent_solution_clobbering():
    models, removals = compare_strategies(create_cross_coupled_model())

    assert removals == {
        "GrowMatch-style local A": ["R_B"],
        "GrowMatch-style local B": ["R_A"],
        "GrowMatch-style union": ["R_A", "R_B"],
        "CROP joint": ["A_BYPASS_1", "A_BYPASS_2", "B_BYPASS_1", "B_BYPASS_2"],
    }

    matches = build_results(models).pivot(
        index="strategy",
        columns="condition",
        values="matches_observation",
    )
    assert matches.loc["GrowMatch-style independent"].to_dict() == {
        "growth_a": False,
        "nogrowth_a": True,
        "growth_b": False,
        "nogrowth_b": True,
    }
    assert matches.loc["CROP joint"].all()
