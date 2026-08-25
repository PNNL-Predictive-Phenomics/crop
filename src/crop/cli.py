# -*- coding: utf-8 -*-

"""Command line interface for :mod:`crop`.

Why does this file exist, and why not put this in ``__main__``? You might be tempted to import things from ``__main__``
later, but that will cause problems--the code will get executed twice:

- When you run ``python3 -m crop`` python will execute``__main__.py`` as a script.
  That means there won't be any ``crop.__main__`` in ``sys.modules``.
- When you import __main__ it will get executed again (as a module) because
  there's no ``crop.__main__`` in ``sys.modules``.

.. seealso:: https://click.palletsprojects.com/en/8.1.x/setuptools/#setuptools-integration
"""

import json
import logging
from pathlib import Path
from typing import Dict

import click
from cobra.io import load_json_model, read_sbml_model

from .verification import VerificationReport, verify_phenotypes

__all__ = [
    "main",
]

logger = logging.getLogger(__name__)


@click.group()
@click.version_option()
def main():
    """CLI for crop."""


def _load_json_object(path: Path, label: str) -> Dict:
  try:
    with path.open() as file:
      value = json.load(file)
  except (OSError, json.JSONDecodeError) as error:
    raise click.UsageError(f"Unable to read {label} JSON '{path}': {error}") from error
  if not isinstance(value, dict):
    raise click.UsageError(f"{label.capitalize()} JSON must contain an object at its top level.")
  return value


def _load_model(path: Path):
  try:
    if path.suffix.lower() == ".json":
      return load_json_model(str(path))
    if path.suffix.lower() in {".xml", ".sbml"}:
      return read_sbml_model(str(path))
  except Exception as error:
    raise click.UsageError(f"Unable to load model '{path}': {error}") from error
  raise click.UsageError("Model must use a .json, .xml, or .sbml extension.")


def _format_table(report: VerificationReport) -> str:
  headers = ("Condition", "Observed", "Predicted", "FBA status", "Biomass", "Result")
  rows = []
  for result in report.conditions:
    biomass = "n/a" if result.biomass_flux is None else f"{result.biomass_flux:.6g}"
    rows.append(
      (
        result.condition,
        result.observed,
        result.predicted or "n/a",
        result.solver_status,
        biomass,
        "PASS" if result.passed else "FAIL",
      )
    )
  widths = [max(len(str(row[index])) for row in [headers] + rows) for index in range(len(headers))]

  def format_row(row):
    """Format one table row using the calculated column widths."""

    return "  ".join(str(value).ljust(widths[index]) for index, value in enumerate(row))

  separator = "  ".join("-" * width for width in widths)
  summary = f"{sum(result.passed for result in report.conditions)}/{len(report.conditions)} phenotypes passed"
  return "\n".join([format_row(headers), separator, *(format_row(row) for row in rows), summary])


@main.command("verify")
@click.option("--model", "model_path", required=True, type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.option(
  "--phenotypes", "phenotypes_path", required=True, type=click.Path(exists=True, dir_okay=False, path_type=Path)
)
@click.option("--media", "media_path", required=True, type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.option("--biomass-rxn")
@click.option("--minimum-growth", default=2.0, show_default=True, type=float)
@click.option("--maximum-nogrowth", default=1.0, show_default=True, type=float)
@click.option("--tolerance", default=1e-6, show_default=True, type=float)
@click.option("--solver", help="COBRA/optlang solver used for FBA.")
@click.option("--format", "output_format", type=click.Choice(["table", "json"]), default="table", show_default=True)
@click.option("--output", type=click.Path(dir_okay=False, path_type=Path))
@click.pass_context
def verify_command(
  context,
  model_path,
  phenotypes_path,
  media_path,
  biomass_rxn,
  minimum_growth,
  maximum_nogrowth,
  tolerance,
  solver,
  output_format,
  output,
):
  """Verify observed phenotypes on an already-modified model using FBA."""

  model = _load_model(model_path)
  phenotypes = _load_json_object(phenotypes_path, "phenotypes")
  media = _load_json_object(media_path, "media")
  try:
    report = verify_phenotypes(
      model,
      phenotypes,
      media,
      biomass_rxn=biomass_rxn,
      minimum_growth=minimum_growth,
      maximum_nogrowth=maximum_nogrowth,
      tolerance=tolerance,
      solver=solver,
    )
  except ValueError as error:
    raise click.UsageError(str(error)) from error

  rendered = json.dumps(report.to_dict(), indent=2) if output_format == "json" else _format_table(report)
  if output is None:
    click.echo(rendered)
  else:
    try:
      output.parent.mkdir(parents=True, exist_ok=True)
      output.write_text(f"{rendered}\n")
    except OSError as error:
      raise click.UsageError(f"Unable to write report '{output}': {error}") from error

  if not report.passed:
    context.exit(1)


if __name__ == "__main__":
    main()
