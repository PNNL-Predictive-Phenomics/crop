# -*- coding: utf-8 -*-

"""Main code."""
import numpy as np
from cvxpy import Minimize, Problem, Variable, diag
import pandas as pd
import cobra
import numpy as np
from cvxpy import Minimize, Problem, Variable, diag


def run_crop_algorithm(
    test_crop_model, phenotype_data, media_conditions, biomass_rxn='BIOMASS', maximum_nogrowth=1, minimum_growth=2,
    atp_maintenance_rxn="ATPM", atp_maintenance_lower_bound=2 
):
    """
    Run the CROP algorithm on the test model with the given phenotype data and media conditions.
    """
    phenotype_conditions = build_phenotype_conditions(media_conditions, phenotype_data)
    # Set the media conditions
    S = cobra.util.array.create_stoichiometric_matrix(test_crop_model, array_type="DataFrame")
    nmets, nrxns = S.shape
    z = Variable(nrxns, boolean=True)
    r = Variable(nrxns)
    m = Variable(nmets)
    v_nogrowth = Variable(nrxns)
    v_growth = Variable(nrxns)
    omega = 1000
    weights = np.ones(nrxns)
    lower_bound_growth = pd.Series(
        {
            rxn.id: (
                -phenotype_conditions["growth"][rxn.id]
                if rxn.id in phenotype_conditions["growth"]
                else 0
            )
            for rxn in test_crop_model.reactions
        }
    )
    if atp_maintenance_rxn in lower_bound_growth.index:
        lower_bound_growth[atp_maintenance_rxn] = atp_maintenance_lower_bound
    lower_bound_nogrowth = pd.Series(
        {
            rxn.id: (
                -phenotype_conditions["nogrowth"][rxn.id]
                if rxn.id in phenotype_conditions["nogrowth"]
                else 0
            )
            for rxn in test_crop_model.reactions
        }
    )
    upper_bound_growth = pd.Series(
        {
            rxn.id: minimum_growth if rxn.id not in phenotype_conditions["growth"] else 0
            for rxn in test_crop_model.reactions
        }
    )
    upper_bound_nogrowth = pd.Series(
        {
            rxn.id: maximum_nogrowth if rxn.id not in phenotype_conditions["nogrowth"] else 0
            for rxn in test_crop_model.reactions
        }
    )
    nogrowth_carbon_source = [
        carbon_source for carbon_source in phenotype_conditions["nogrowth"] if carbon_source != 0
    ][0]
    # growth_carbon_source = [carbon_source for carbon_source in media_conditions['growth']
    #                         if carbon_source !=0][0]
    nogrowth_carbon_source_idx = lower_bound_nogrowth.index.get_loc(nogrowth_carbon_source)
    biomass_idx = lower_bound_nogrowth.index.get_loc(biomass_rxn)
    c = np.zeros(nrxns)
    c[biomass_idx] = 1

    # $$\begin{equation}\begin{array}{l}
    problem = Problem(
        Minimize(weights.T @ (1 - z)),[ # \min_z weights^T(1-z) \\  
            # \begin{array}{lll}
                v_nogrowth[biomass_idx] == lower_bound_nogrowth[nogrowth_carbon_source] * r[nogrowth_carbon_source_idx],  # & v_{biomass} = U_{nogrowth,lactose}r_{lactose} \\
                S.values @ v_nogrowth == 0,  # & Sv_{nogrowth}= 0 & \text{inner problem} \\
                diag(lower_bound_nogrowth.values) @ z <= v_nogrowth,
                v_nogrowth <= diag(upper_bound_nogrowth.values)@ z,  # & 0\leq v_i\leq U_{nogrowth,i}\cdot z_i & \text{every carbon source $i$ not in the nogrowth media has $U_{nogrowth,i} = 0$. \\ Every carbon source $j$ in the nogrowth media has $U_{nogrowth,j} > 0$} \\
                S.T.values @ m + c == r,  # & S^Tm +c = r \\
                r <= omega * (1 - z),  # & r_i\leq\Omega_i\cdot(1-z_i) & \text{for $i\neq$ lactose. This constraint ensures that lactose uptake is the only tight constraint in the model. Every other reaction with a tight constraint cannot be part of the model.}\\
                r[nogrowth_carbon_source_idx] <= 0,  # & r_{lactose} \geq 0 \\
            # \end{array}\\
            v_nogrowth[biomass_idx] <= maximum_nogrowth,  # v_{biomass} \leq\text{minimal growth} \\
            S.values @ v_growth == 0,  # Sw= 0 \\
            diag(lower_bound_growth.values) @ z <= v_growth,
            v_growth <= diag(upper_bound_growth.values) @ z,  # 0\leq w_i\leq U_{growth,i}\cdot z_i \\
            v_growth[biomass_idx] >= minimum_growth,  # w_{biomass} \geq \text{minimal growth} \\
        ] #+ non_removable_z ,# z_i = 1 \text{ for all non-removable reactions } i \\
    )  # w_{ATP} \geq \text{atp maintenance} \\
    # z\in \{0,1\} \\
    # \end{array}\end{equation}$$
    problem.solve()
    solution = pd.DataFrame(
        {"r": r.value, "z": z.value, "v_nogrowth": v_nogrowth.value, "v_growth": v_growth.value},
        index=S.columns,
    )
    if problem.status != "optimal":
        raise ValueError("Infeasible problem")
    suggested_removals = solution['z'].index[solution['z'].eq(0)].tolist()
    return suggested_removals, solution


def build_phenotype_conditions(
    media_conditions: dict,
    phenotype_data: dict,
    growth_condition: str | None = None,
    nogrowth_condition: str | None = None,
) -> dict:
    """
    Transform media_conditions and phenotype_data into phenotype_conditions.

    - media_conditions: mapping of condition -> {exchange_rxn: lower_bound}
      e.g., negative values allow uptake (COBRA convention).
    - phenotype_data: mapping of condition -> {"observed": "...", "predicted": "..."}.

    Selects:
      growth: observed == "growth" and predicted == "growth"
      nogrowth: observed == "no_growth" and predicted == "growth"

    Flips sign for uptakes: phenotype_conditions stores magnitudes (positive).
    Any non-uptake (>= 0) in media becomes 0.0 in phenotype_conditions.

    Optional growth_condition / nogrowth_condition let you pick specific names;
    otherwise the first match is used.
    """
    def obs_key(d: dict) -> str:
        # be tolerant if data uses "observation" instead of "observed"
        return d.get("observed", d.get("observation", ""))

    # Find candidate condition names
    growth_candidates = [
        name for name, ph in phenotype_data.items()
        if obs_key(ph) == "growth" and ph.get("predicted") == "growth" and name in media_conditions
    ]
    nogrowth_candidates = [
        name for name, ph in phenotype_data.items()
        if obs_key(ph) in ("no_growth", "nogrowth") and ph.get("predicted") == "growth" and name in media_conditions
    ]

    if growth_condition is None:
        if not growth_candidates:
            raise ValueError("No growth condition found where observed==growth and predicted==growth.")
        growth_condition = growth_candidates[0]
    elif growth_condition not in media_conditions:
        raise KeyError(f"Growth condition '{growth_condition}' not in media_conditions.")

    if nogrowth_condition is None:
        if not nogrowth_candidates:
            raise ValueError("No nogrowth condition found where observed==no_growth and predicted==growth.")
        nogrowth_condition = nogrowth_candidates[0]
    elif nogrowth_condition not in media_conditions:
        raise KeyError(f"Nogrowth condition '{nogrowth_condition}' not in media_conditions.")

    # Include all exchange rxns that appear in any media condition so keys are consistent
    all_exchanges = set().union(*[mc.keys() for mc in media_conditions.values()])

    def to_pheno(m: dict) -> dict:
        # Flip sign for allowed uptakes (negative in media -> positive magnitude here), else 0.0
        return {rxn: float(-m.get(rxn, 0)) if m.get(rxn, 0) < 0 else 0.0 for rxn in all_exchanges}

    growth_map = to_pheno(media_conditions[growth_condition])
    nogrowth_map = to_pheno(media_conditions[nogrowth_condition])

    return {"growth": growth_map, "nogrowth": nogrowth_map}

# Example with your data:
# phenotype_conditions = build_phenotype_conditions(media_conditions(), phenotype_data())
# print(phenotype_conditions)