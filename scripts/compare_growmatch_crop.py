"""Compare single-phenotype reconciliation with joint multi-phenotype CROP.

This uses CROP's optimizer for both arms so that only the condition scope differs.
The independent arm represents GrowMatch's one-phenotype-at-a-time strategy; it
does not invoke the separate GrowMatch implementation or claim numerical parity
with its genome-scale model.
"""

import argparse
import json
from pathlib import Path
from typing import Dict, Iterable, Set

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from cobra import Metabolite, Model, Reaction

from crop import run_crop_algorithm, verify_phenotypes

MINIMUM_GROWTH = 2.0
MAXIMUM_NOGROWTH = 1.0

MEDIA = {
    "growth_a": {"EX_a": -10.0, "EX_helper_a": -10.0, "EX_b": 0.0, "EX_helper_b": 0.0},
    "nogrowth_a": {"EX_a": -10.0, "EX_helper_a": 0.0, "EX_b": 0.0, "EX_helper_b": 0.0},
    "growth_b": {"EX_a": 0.0, "EX_helper_a": 0.0, "EX_b": -10.0, "EX_helper_b": -10.0},
    "nogrowth_b": {"EX_a": 0.0, "EX_helper_a": 0.0, "EX_b": -10.0, "EX_helper_b": 0.0},
}

PHENOTYPES = {
    "growth_a": {"observed": "growth", "predicted": "growth"},
    "nogrowth_a": {"observed": "no_growth", "predicted": "growth"},
    "growth_b": {"observed": "growth", "predicted": "growth"},
    "nogrowth_b": {"observed": "no_growth", "predicted": "growth"},
}


def _reaction(reaction_id: str, metabolites: Dict[Metabolite, float]) -> Reaction:
    reaction = Reaction(reaction_id)
    reaction.lower_bound = 0.0
    reaction.upper_bound = 100.0
    reaction.add_metabolites(metabolites)
    return reaction


def create_cross_coupled_model() -> Model:
    """Build a model where local fixes conflict but a joint fix exists."""

    model = Model("single_vs_multi_phenotype")
    a_e = Metabolite("a_e", compartment="e")
    a_c = Metabolite("a_c", compartment="c")
    helper_a_e = Metabolite("helper_a_e", compartment="e")
    helper_a_c = Metabolite("helper_a_c", compartment="c")
    b_e = Metabolite("b_e", compartment="e")
    b_c = Metabolite("b_c", compartment="c")
    helper_b_e = Metabolite("helper_b_e", compartment="e")
    helper_b_c = Metabolite("helper_b_c", compartment="c")
    intermediate_a = Metabolite("intermediate_a", compartment="c")
    intermediate_b = Metabolite("intermediate_b", compartment="c")
    biomass_precursor = Metabolite("biomass_precursor", compartment="c")

    exchanges = [
        _reaction("EX_a", {a_e: -1}),
        _reaction("EX_helper_a", {helper_a_e: -1}),
        _reaction("EX_b", {b_e: -1}),
        _reaction("EX_helper_b", {helper_b_e: -1}),
    ]
    for exchange in exchanges:
        exchange.lower_bound = -10.0

    model.add_reactions(
        exchanges
        + [
            _reaction("A_TRANSPORT", {a_e: -1, a_c: 1}),
            _reaction("HELPER_A_TRANSPORT", {helper_a_e: -1, helper_a_c: 1}),
            _reaction("B_TRANSPORT", {b_e: -1, b_c: 1}),
            _reaction("HELPER_B_TRANSPORT", {helper_b_e: -1, helper_b_c: 1}),
            _reaction("A_LEGITIMATE", {a_c: -1, helper_a_c: -1, intermediate_a: 1}),
            _reaction("R_A", {intermediate_a: -1, biomass_precursor: 1}),
            _reaction("B_LEGITIMATE", {b_c: -1, helper_b_c: -1, intermediate_b: 1}),
            _reaction("R_B", {intermediate_b: -1, biomass_precursor: 1}),
            _reaction("A_BYPASS_1", {a_c: -1, intermediate_b: 1}),
            _reaction("A_BYPASS_2", {a_c: -1, intermediate_b: 1}),
            _reaction("B_BYPASS_1", {b_c: -1, intermediate_a: 1}),
            _reaction("B_BYPASS_2", {b_c: -1, intermediate_a: 1}),
            _reaction("BIOMASS", {biomass_precursor: -1}),
        ]
    )
    model.objective = "BIOMASS"
    return model


def apply_removals(model: Model, removals: Iterable[str]) -> Model:
    """Return a model copy with the selected reactions disabled."""

    corrected = model.copy()
    for reaction_id in removals:
        corrected.reactions.get_by_id(reaction_id).bounds = (0.0, 0.0)
    return corrected


def _run_crop(
    model: Model,
    phenotypes: Dict[str, Dict[str, str]],
    show_solver_output: bool = False,
) -> Set[str]:
    removals, _ = run_crop_algorithm(
        model,
        phenotypes,
        MEDIA,
        minimum_growth=MINIMUM_GROWTH,
        maximum_nogrowth=MAXIMUM_NOGROWTH,
        verbose=show_solver_output,
    )
    return removals


def compare_strategies(model: Model, show_solver_output: bool = False):
    """Run unmodified, independent, and joint reconciliation strategies."""

    local_a = _run_crop(
        model,
        {name: PHENOTYPES[name] for name in ("growth_a", "nogrowth_a")},
        show_solver_output,
    )
    local_b = _run_crop(
        model,
        {name: PHENOTYPES[name] for name in ("growth_b", "nogrowth_b")},
        show_solver_output,
    )
    independent_removals = local_a | local_b
    joint_removals = _run_crop(model, PHENOTYPES, show_solver_output)

    models = {
        "Unmodified": model,
        "GrowMatch-style independent": apply_removals(model, independent_removals),
        "CROP joint": apply_removals(model, joint_removals),
    }
    removal_sets = {
        "GrowMatch-style local A": sorted(local_a),
        "GrowMatch-style local B": sorted(local_b),
        "GrowMatch-style union": sorted(independent_removals),
        "CROP joint": sorted(joint_removals),
    }
    return models, removal_sets


def build_results(models: Dict[str, Model]) -> pd.DataFrame:
    """Independently verify every strategy under every phenotype condition."""

    rows = []
    for strategy, model in models.items():
        report = verify_phenotypes(
            model,
            PHENOTYPES,
            MEDIA,
            minimum_growth=MINIMUM_GROWTH,
            maximum_nogrowth=MAXIMUM_NOGROWTH,
        )
        for result in report.conditions:
            rows.append(
                {
                    "strategy": strategy,
                    "condition": result.condition,
                    "observed": result.observed,
                    "biomass_flux": result.biomass_flux or 0.0,
                    "matches_observation": result.passed,
                }
            )
    return pd.DataFrame(rows)


def plot_fluxes(results: pd.DataFrame, output_path: Path) -> None:
    """Plot FBA biomass flux by strategy and condition."""

    strategies = list(results["strategy"].drop_duplicates())
    conditions = list(PHENOTYPES)
    positions = np.arange(len(conditions))
    width = 0.24
    colors = ["#6c757d", "#d55e00", "#0072b2"]

    figure, axis = plt.subplots(figsize=(12, 5.5))
    for index, strategy in enumerate(strategies):
        subset = results[results["strategy"] == strategy].set_index("condition")
        bars = axis.bar(
            positions + (index - 1) * width,
            subset.loc[conditions, "biomass_flux"],
            width,
            label=strategy,
            color=colors[index],
        )
        axis.bar_label(bars, fmt="%.0f", padding=3)
    axis.axhline(MINIMUM_GROWTH, color="#009e73", linestyle="--", label="minimum growth")
    axis.axhline(MAXIMUM_NOGROWTH, color="#cc79a7", linestyle=":", label="maximum no-growth")
    axis.set_xticks(positions, conditions)
    axis.set_ylabel("FBA biomass flux")
    axis.set_title("Joint constraints preserve growth while suppressing false growth")
    axis.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    figure.tight_layout(rect=(0, 0, 0.82, 1))
    figure.savefig(output_path, dpi=180)
    plt.close(figure)


def plot_matches(results: pd.DataFrame, output_path: Path) -> None:
    """Plot whether each strategy reproduces each observed phenotype."""

    matrix = results.pivot(index="strategy", columns="condition", values="matches_observation")
    matrix = matrix.loc[
        ["Unmodified", "GrowMatch-style independent", "CROP joint"],
        list(PHENOTYPES),
    ]
    figure, axis = plt.subplots(figsize=(10, 3.5))
    axis.imshow(matrix.to_numpy(dtype=int), cmap="RdYlGn", vmin=0, vmax=1, aspect="auto")
    axis.set_xticks(range(len(matrix.columns)), matrix.columns)
    axis.set_yticks(range(len(matrix.index)), matrix.index)
    axis.set_title("Observed phenotype reproduced by FBA")
    for row in range(matrix.shape[0]):
        for column in range(matrix.shape[1]):
            value = bool(matrix.iloc[row, column])
            axis.text(
                column,
                row,
                "MATCH" if value else "DIFFER",
                ha="center",
                va="center",
                fontweight="bold",
            )
    figure.tight_layout()
    figure.savefig(output_path, dpi=180)
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("comparison_output"))
    parser.add_argument("--show-solver-output", action="store_true")
    arguments = parser.parse_args()
    arguments.output_dir.mkdir(parents=True, exist_ok=True)

    models, removal_sets = compare_strategies(
        create_cross_coupled_model(),
        show_solver_output=arguments.show_solver_output,
    )
    results = build_results(models)
    results.to_csv(arguments.output_dir / "phenotype_comparison.csv", index=False)
    (arguments.output_dir / "removals.json").write_text(json.dumps(removal_sets, indent=2) + "\n")
    plot_fluxes(results, arguments.output_dir / "biomass_fluxes.png")
    plot_matches(results, arguments.output_dir / "phenotype_matches.png")

    print(json.dumps(removal_sets, indent=2))
    print()
    print(results.to_string(index=False))
    print(f"\nWrote comparison artifacts to {arguments.output_dir}")


if __name__ == "__main__":
    main()
