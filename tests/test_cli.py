"""Command-line tests for post-correction phenotype verification."""

import json

from click.testing import CliRunner
from cobra.io import save_json_model, write_sbml_model

from crop.cli import main
from tests.test_crop_oracle import create_unique_oracle_model

MEDIA = {
    "lactose_with_helper": {"EX_lac": -10.0, "EX_helper": -10.0},
    "lactose_only": {"EX_lac": -10.0, "EX_helper": 0.0},
}
PHENOTYPES = {
    "lactose_with_helper": {"observed": "growth", "predicted": "growth"},
    "lactose_only": {"observed": "no_growth", "predicted": "growth"},
}


def _write_inputs(tmp_path, corrected=True, model_format="json"):
    model = create_unique_oracle_model()
    if corrected:
        model.reactions.get_by_id("LAC_bypass").bounds = (0.0, 0.0)
    model_path = tmp_path / f"model.{model_format}"
    if model_format == "json":
        save_json_model(model, str(model_path))
    else:
        write_sbml_model(model, str(model_path))

    phenotypes_path = tmp_path / "phenotypes.json"
    phenotypes_path.write_text(json.dumps(PHENOTYPES))
    media_path = tmp_path / "media.json"
    media_path.write_text(json.dumps(MEDIA))
    return model_path, phenotypes_path, media_path


def _invoke(tmp_path, corrected=True, model_format="json", extra_args=None):
    model_path, phenotypes_path, media_path = _write_inputs(tmp_path, corrected, model_format)
    arguments = [
        "verify",
        "--model",
        str(model_path),
        "--phenotypes",
        str(phenotypes_path),
        "--media",
        str(media_path),
    ]
    return CliRunner().invoke(main, arguments + (extra_args or []))


def test_verify_table_passes_for_corrected_json_model(tmp_path):
    result = _invoke(tmp_path)

    assert result.exit_code == 0, result.output
    assert "lactose_with_helper" in result.output
    assert "2/2 phenotypes passed" in result.output


def test_verify_returns_one_when_an_observation_fails(tmp_path):
    result = _invoke(tmp_path, corrected=False)

    assert result.exit_code == 1
    assert "lactose_only" in result.output
    assert "FAIL" in result.output


def test_verify_json_output_has_stable_summary(tmp_path):
    result = _invoke(tmp_path, extra_args=["--format", "json"])

    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["summary"]["overall_status"] == "pass"
    assert [condition["condition"] for condition in payload["conditions"]] == [
        "lactose_with_helper",
        "lactose_only",
    ]


def test_verify_writes_json_report_file(tmp_path):
    output_path = tmp_path / "reports" / "verification.json"

    result = _invoke(
        tmp_path,
        extra_args=["--format", "json", "--output", str(output_path)],
    )

    assert result.exit_code == 0, result.output
    assert result.output == ""
    assert json.loads(output_path.read_text())["summary"]["failed"] == 0


def test_verify_loads_corrected_sbml_model(tmp_path):
    result = _invoke(tmp_path, model_format="xml")

    assert result.exit_code == 0, result.output


def test_verify_invalid_thresholds_are_usage_error(tmp_path):
    result = _invoke(tmp_path, extra_args=["--minimum-growth", "1", "--maximum-nogrowth", "1"])

    assert result.exit_code == 2
    assert "maximum_nogrowth must be less than minimum_growth" in result.output
