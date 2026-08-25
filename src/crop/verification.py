"""Verify observed phenotypes on a corrected metabolic model."""

from dataclasses import asdict, dataclass
from math import isfinite
from typing import Dict, List, Optional

import cobra


@dataclass(frozen=True)
class PhenotypeVerification:
    """FBA result for one phenotype condition."""

    condition: str
    observed: str
    predicted: Optional[str]
    solver_status: str
    biomass_flux: Optional[float]
    threshold: float
    passed: bool
    message: str

    def to_dict(self) -> Dict[str, object]:
        """Return a JSON-serializable representation."""

        return asdict(self)


@dataclass(frozen=True)
class VerificationReport:
    """Post-correction verification results."""

    conditions: List[PhenotypeVerification]

    @property
    def passed(self) -> bool:
        """Return whether every observed phenotype was reproduced."""

        return all(result.passed for result in self.conditions)

    def to_dict(self) -> Dict[str, object]:
        """Return a JSON-serializable report."""

        passed_count = sum(result.passed for result in self.conditions)
        return {
            "summary": {
                "total": len(self.conditions),
                "passed": passed_count,
                "failed": len(self.conditions) - passed_count,
                "overall_status": "pass" if self.passed else "fail",
            },
            "conditions": [result.to_dict() for result in self.conditions],
        }


def _normalize_phenotype(value: str) -> str:
    normalized = value.strip().lower().replace("-", "_").replace(" ", "_")
    if normalized == "nogrowth":
        normalized = "no_growth"
    if normalized not in {"growth", "no_growth"}:
        raise ValueError(f"Unsupported phenotype '{value}'. Expected 'growth' or 'no_growth'.")
    return normalized


def verify_phenotypes(
    model: cobra.Model,
    phenotype_data: Dict[str, Dict[str, str]],
    media_conditions: Dict[str, Dict[str, float]],
    biomass_rxn: Optional[str] = None,
    minimum_growth: float = 2.0,
    maximum_nogrowth: float = 1.0,
    tolerance: float = 1e-6,
    solver: Optional[str] = None,
) -> VerificationReport:
    """Run FBA for every observed phenotype on an already-corrected model."""

    if minimum_growth < 0 or maximum_nogrowth < 0 or tolerance < 0:
        raise ValueError("Growth thresholds and tolerance must be non-negative.")
    if maximum_nogrowth >= minimum_growth:
        raise ValueError("maximum_nogrowth must be less than minimum_growth.")
    if not phenotype_data:
        raise ValueError("No phenotype conditions were provided.")

    results: List[PhenotypeVerification] = []
    for condition, phenotype in phenotype_data.items():
        if condition not in media_conditions:
            raise ValueError(f"No media definition found for condition '{condition}'.")

        observed_value = phenotype.get("observed", phenotype.get("observation"))
        if observed_value is None:
            raise ValueError(f"Condition '{condition}' has no observed phenotype.")
        observed = _normalize_phenotype(observed_value)
        predicted_value = phenotype.get("predicted")
        predicted = _normalize_phenotype(predicted_value) if predicted_value is not None else None

        condition_model = model.copy()
        if biomass_rxn is not None:
            if biomass_rxn not in condition_model.reactions:
                raise ValueError(f"Biomass reaction '{biomass_rxn}' was not found in the model.")
            condition_model.objective = biomass_rxn
        if solver is not None:
            try:
                condition_model.solver = solver
            except Exception as error:
                raise ValueError(f"Unable to use COBRA solver '{solver}': {error}") from error

        for exchange in condition_model.exchanges:
            exchange.lower_bound = 0.0
        for reaction_id, lower_bound in media_conditions[condition].items():
            if reaction_id not in condition_model.reactions:
                raise ValueError(
                    f"Media for condition '{condition}' references unknown reaction '{reaction_id}'."
                )
            condition_model.reactions.get_by_id(reaction_id).lower_bound = float(lower_bound)

        solution = condition_model.optimize()
        flux = None
        if solution.status == "optimal" and solution.objective_value is not None:
            candidate_flux = float(solution.objective_value)
            if isfinite(candidate_flux):
                flux = candidate_flux

        threshold = minimum_growth if observed == "growth" else maximum_nogrowth
        if observed == "growth":
            passed = flux is not None and flux >= minimum_growth - tolerance
            message = (
                "growth reproduced"
                if passed
                else f"biomass did not reach minimum growth {minimum_growth:g}"
            )
        else:
            passed = flux is None or flux <= maximum_nogrowth + tolerance
            message = (
                "no-growth reproduced"
                if passed
                else f"biomass exceeded maximum no-growth {maximum_nogrowth:g}"
            )

        results.append(
            PhenotypeVerification(
                condition=condition,
                observed=observed,
                predicted=predicted,
                solver_status=str(solution.status),
                biomass_flux=flux,
                threshold=threshold,
                passed=passed,
                message=message,
            )
        )

    return VerificationReport(results)
