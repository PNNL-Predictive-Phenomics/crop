"""Tests for post-correction phenotype verification."""

from crop.verification import verify_phenotypes
from tests.test_crop_oracle import create_unique_oracle_model

MEDIA = {
    "lactose_with_helper": {"EX_lac": -10.0, "EX_helper": -10.0},
    "lactose_only": {"EX_lac": -10.0, "EX_helper": 0.0},
}
PHENOTYPES = {
    "lactose_with_helper": {"observed": "growth", "predicted": "growth"},
    "lactose_only": {"observed": "no_growth", "predicted": "growth"},
}


def test_corrected_model_reproduces_all_observed_phenotypes():
    model = create_unique_oracle_model()
    model.reactions.get_by_id("LAC_bypass").bounds = (0.0, 0.0)

    report = verify_phenotypes(model, PHENOTYPES, MEDIA)

    assert report.passed
    assert [result.passed for result in report.conditions] == [True, True]
    assert report.to_dict()["summary"] == {
        "total": 2,
        "passed": 2,
        "failed": 0,
        "overall_status": "pass",
    }


def test_uncorrected_model_reports_false_growth_without_mutating_model():
    model = create_unique_oracle_model()
    original_bounds = {reaction.id: reaction.bounds for reaction in model.reactions}

    report = verify_phenotypes(model, PHENOTYPES, MEDIA)

    failures = [result.condition for result in report.conditions if not result.passed]
    assert failures == ["lactose_only"]
    assert {reaction.id: reaction.bounds for reaction in model.reactions} == original_bounds
