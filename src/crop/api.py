# -*- coding: utf-8 -*-

"""Main code."""
from typing import Dict, Set, Tuple, List, Optional, Callable
import numpy as np
from cvxpy import Minimize, Problem, Variable, diag
import pandas as pd
import cobra
import numpy as np
from cvxpy import Minimize, Problem, Variable, diag


def run_crop_algorithm(
    test_crop_model: cobra.Model,
    phenotype_data: Dict[str, Dict[str, str]],
    media_conditions: Dict[str, Dict[str, float]],
    biomass_rxn: str = 'BIOMASS',
    maximum_nogrowth: float = 1.0,
    minimum_growth: float = 2.0,
    atp_maintenance_rxn: str = "ATPM",
    atp_maintenance_lower_bound: float = 2.0,
    solver: str = 'SCIPY'
) -> Tuple[Set[str], pd.DataFrame]:
    """
    Run the CROP algorithm to identify reactions to remove for fixing phenotype predictions.
    
    This implements a single-level MILP formulation that solves the bi-level optimization problem:
    
    .. math::
        \\begin{equation}\\begin{array}{l}
        \\min_z weights^T(1-z) \\\\
        \\begin{array}{lll}
        & v_{biomass} = U_{nogrowth,carbon}r_{carbon} \\\\
        & Sv_{nogrowth}= 0 & \\text{inner problem} \\\\
        & L_{nogrowth,i}\\cdot z_i\\leq v_i\\leq U_{nogrowth,i}\\cdot z_i \\\\
        & S^Tm +c = r \\\\
        & r_i\\leq\\Omega_i\\cdot(1-z_i) \\\\
        & r_{carbon} \\leq 0 \\\\
        \\end{array}\\\\
        v_{biomass} \\leq\\text{maximum nogrowth} \\\\
        Sw= 0 \\\\
        L_{growth,i}\\cdot z_i \\leq w_i\\leq U_{growth,i}\\cdot z_i \\\\
        w_{biomass} \\geq \\text{minimum growth} \\\\
        w_{ATP} \\geq \\text{atp maintenance} \\\\
        z\\in \\{0,1\\} \\\\
        \\end{array}\\end{equation}
    
    Parameters
    ----------
    test_crop_model : cobra.Model
        The metabolic model to analyze. Should contain reactions with IDs matching
        those in media_conditions.
    phenotype_data : Dict[str, Dict[str, str]]
        Phenotype observations and predictions for each condition. Format:
        {condition_name: {"observed": "growth"|"no_growth", "predicted": "growth"|"no_growth"}}
        The algorithm identifies conditions where observed="no_growth" and predicted="growth"
        as the nogrowth condition, and where both are "growth" as the growth condition.
    media_conditions : Dict[str, Dict[str, float]]
        Media compositions for each condition. Format:
        {condition_name: {exchange_rxn_id: lower_bound}}
        Negative values indicate allowed uptake (COBRA convention).
        Example: {"glucose": {"EX_glc": -10.0, "EX_lac": 0.0}}
    biomass_rxn : str, default='BIOMASS'
        Reaction ID of the biomass objective function in the model.
        Used as the growth indicator: $v_{biomass}$ and $w_{biomass}$.
    maximum_nogrowth : float, default=1.0
        Maximum allowed biomass flux in the nogrowth condition.
        Corresponds to $\\text{maximum nogrowth}$ in the formulation.
        Should be set below the threshold that defines growth.
    minimum_growth : float, default=2.0
        Minimum required biomass flux in the growth condition.
        Corresponds to $\\text{minimum growth}$ in the formulation.
        Should be set above the threshold that defines growth.
    atp_maintenance_rxn : str, default="ATPM"
        Reaction ID for ATP maintenance. If present, its lower bound will be
        set to atp_maintenance_lower_bound in growth conditions.
        Corresponds to $w_{ATP} \\geq \\text{atp maintenance}$ in the formulation.
    atp_maintenance_lower_bound : float, default=2.0
        Minimum ATP maintenance flux required in growth conditions.
        Only applied if atp_maintenance_rxn exists in the model.
    solver : str, default='SCIPY'
        CVXPY solver to use. Options include 'SCIPY', 'GUROBI', 'CPLEX', etc.
        Some problems may require commercial solvers for better performance.
    
    Returns
    -------
    Tuple[Set[str], pd.DataFrame]
        A tuple containing:
        - Set of reaction IDs suggested for removal (where z_i = 0)
        - DataFrame with solution details including columns:
          * 'r': dual variable values (reduced costs)
          * 'z': binary inclusion variables (1 = keep, 0 = remove)
          * 'v_nogrowth': flux values in nogrowth condition
          * 'v_growth': flux values in growth condition
    
    Raises
    ------
    ValueError
        If the optimization problem is infeasible or no valid growth/nogrowth
        conditions are found in phenotype_data.
    
    Notes
    -----
    The algorithm works by:
    1. Identifying a nogrowth condition where the model incorrectly predicts growth
    2. Setting up a bi-level MILP that finds reactions whose removal would:
       - Prevent growth in the nogrowth condition (biomass ≤ maximum_nogrowth)
       - Maintain growth in the growth condition (biomass ≥ minimum_growth)
    3. Minimizing the number of reactions removed (weighted by evidence)
    
    The variable z_i ∈ {0,1} determines whether reaction i is included:
    - z_i = 1: reaction is kept in the model
    - z_i = 0: reaction should be removed
    
    The dual variable r_i represents the shadow price on flux bounds:
    - r_i > 0: upper bound is constraining
    - r_i < 0: lower bound is constraining
    - r_i ≈ 0: reaction is not at a bound
    
    Examples
    --------
    >>> phenotype_data = {
    ...     "glucose": {"observed": "growth", "predicted": "growth"},
    ...     "lactose": {"observed": "no_growth", "predicted": "growth"}
    ... }
    >>> media_conditions = {
    ...     "glucose": {"EX_glc": -10.0, "EX_lac": 0.0},
    ...     "lactose": {"EX_glc": 0.0, "EX_lac": -10.0}
    ... }
    >>> removals, solution = run_crop_algorithm(model, phenotype_data, media_conditions)
    >>> print(f"Suggested removals: {removals}")
    """
    np.random.seed(42)  # For reproducibility
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
    nogrowth_carbon_sources = [
        carbon_source for carbon_source, uptake_rate in phenotype_conditions["nogrowth"].items()
        if uptake_rate != 0
    ]
    print(f"Phenotype conditions: {phenotype_conditions}")
    print(f"Media conditions: {media_conditions}")
    print(f"Phenotype data: {phenotype_data}")
    print(f"nogrowth carbon sources: {nogrowth_carbon_sources}")
    nogrowth_carbon_source = sorted(nogrowth_carbon_sources)[0]
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
    problem.solve(verbose=True, solver=solver)
    solution = pd.DataFrame(
        {"r": r.value, "z": z.value, "v_nogrowth": v_nogrowth.value, "v_growth": v_growth.value},
        index=S.columns,
    )
    if problem.status != "optimal":
        raise ValueError("Infeasible problem")
    zero_z = np.isclose(solution['z'], 0)
    print(f"Zero Z: ")
    print(solution['z'])
    suggested_removals = set(solution['z'].index[zero_z].tolist())
    return suggested_removals, solution


def build_phenotype_conditions(
    media_conditions: Dict[str, Dict[str, float]],
    phenotype_data: Dict[str, Dict[str, str]],
    growth_condition: Optional[str] = None,
    nogrowth_condition: Optional[str] = None,
) -> Dict[str, Dict[str, float]]:
    """
    Transform media conditions and phenotype data into standardized phenotype conditions.
    
    This function processes media compositions and phenotype observations to create
    standardized uptake rates for growth and nogrowth conditions. It converts COBRA's
    negative-uptake convention to positive magnitudes for use in the CROP algorithm.
    
    Parameters
    ----------
    media_conditions : Dict[str, Dict[str, float]]
        Media compositions for each condition. Format:
        {condition_name: {exchange_rxn_id: lower_bound}}
        Negative values indicate allowed uptake (COBRA convention).
        Non-negative values mean the nutrient is absent or blocked.
        Example: {"glucose": {"EX_glc": -10.0, "EX_lac": 0.0}}
    phenotype_data : Dict[str, Dict[str, str]]
        Phenotype observations and predictions. Format:
        {condition_name: {"observed": "growth"|"no_growth", "predicted": "growth"|"no_growth"}}
        The function identifies:
        - growth condition: observed="growth" AND predicted="growth" (correct prediction)
        - nogrowth condition: observed="no_growth" AND predicted="growth" (incorrect prediction to fix)
    growth_condition : Optional[str], default=None
        Specific condition name to use as the growth condition.
        If None, automatically selects the first condition where
        observed="growth" and predicted="growth".
    nogrowth_condition : Optional[str], default=None
        Specific condition name to use as the nogrowth condition.
        If None, automatically selects the first condition where
        observed="no_growth" and predicted="growth".
    
    Returns
    -------
    Dict[str, Dict[str, float]]
        Standardized phenotype conditions with format:
        {"growth": {exchange_rxn: uptake_rate}, "nogrowth": {exchange_rxn: uptake_rate}}
        where uptake_rate is:
        - positive magnitude for allowed uptakes (flipped from negative in media_conditions)
        - 0.0 for blocked or absent nutrients
        All exchange reactions from all media conditions are included for consistency.
    
    Raises
    ------
    ValueError
        If no valid growth condition is found (observed="growth" and predicted="growth")
        or no valid nogrowth condition is found (observed="no_growth" and predicted="growth").
    KeyError
        If a specified growth_condition or nogrowth_condition is not in media_conditions.
    
    Notes
    -----
    The function:
    1. Identifies candidate growth/nogrowth conditions based on phenotype data
    2. Selects specific conditions (first match or user-specified)
    3. Converts media bounds to positive uptake magnitudes
    4. Ensures all exchange reactions appear in both growth and nogrowth dictionaries
    
    The COBRA convention uses negative lower bounds for uptake:
    - EX_glc = -10.0: allows up to 10 mmol/gDW/h glucose uptake
    - EX_glc = 0.0: no glucose uptake allowed
    
    This function flips the sign so uptake rates are positive:
    - EX_glc: 10.0 (in phenotype_conditions) corresponds to -10.0 (in media_conditions)
    
    Examples
    --------
    >>> media = {
    ...     "glucose": {"EX_glc": -10.0, "EX_lac": 0.0},
    ...     "lactose": {"EX_glc": 0.0, "EX_lac": -5.0}
    ... }
    >>> phenotypes = {
    ...     "glucose": {"observed": "growth", "predicted": "growth"},
    ...     "lactose": {"observed": "no_growth", "predicted": "growth"}
    ... }
    >>> result = build_phenotype_conditions(media, phenotypes)
    >>> print(result)
    {'growth': {'EX_glc': 10.0, 'EX_lac': 0.0},
     'nogrowth': {'EX_glc': 0.0, 'EX_lac': 5.0}}
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


def get_carbon_source(
    model: cobra.Model,
    phenotype_data: Dict[str, Dict[str, str]],
    media_conditions: Dict[str, Dict[str, float]],
    condition_function: Callable[[Dict[str, Dict[str, str]]], List[str]]
) -> Dict[str, Tuple[int, str]]:
    """
    Identify the carbon source exchange reaction for specified conditions.
    
    For each condition returned by condition_function, this identifies the primary
    carbon source (the exchange reaction with a negative lower bound in the media).
    
    Parameters
    ----------
    model : cobra.Model
        The metabolic model containing the reactions.
    phenotype_data : Dict[str, Dict[str, str]]
        Phenotype observations and predictions for each condition.
    media_conditions : Dict[str, Dict[str, float]]
        Media compositions. Format: {condition: {exchange_rxn: lower_bound}}
    condition_function : callable
        Function that takes phenotype_data and returns a list of condition names.
        Typically get_growth_conditions or get_nogrowth_conditions.
    
    Returns
    -------
    Dict[str, Tuple[int, str]]
        Mapping of condition names to (index, reaction_id) tuples where:
        - index: integer position of the reaction in model.reactions
        - reaction_id: string identifier of the carbon source exchange reaction
    
    Notes
    -----
    This function assumes the first exchange reaction with a negative bound
    in each media condition is the primary carbon source.
    
    Examples
    --------
    >>> media = {"glucose": {"EX_glc": -10.0, "EX_lac": 0.0}}
    >>> result = get_carbon_source(model, phenotypes, media, get_growth_conditions)
    >>> result
    {'glucose': (5, 'EX_glc')}  # EX_glc is at index 5 in model.reactions
    """
    return {
        condition: [
            ([rxn.id for rxn in model.reactions].index(carbon_source), carbon_source) for carbon_source in media_conditions[condition]
        ][0] for condition in condition_function(phenotype_data)
    }

def get_lower_bound_for_conditions(
    model: cobra.Model,
    phenotype_data: Dict[str, Dict[str, str]],
    media_conditions: Dict[str, Dict[str, float]],
    condition_function: Callable[[Dict[str, Dict[str, str]]], List[str]]
) -> pd.DataFrame:
    """
    Construct lower bound matrix for specified conditions.
    
    Creates a DataFrame where each column represents a condition and each row
    represents a reaction, with values being the lower bounds (typically negative
    for exchange reactions allowing uptake).
    
    Parameters
    ----------
    model : cobra.Model
        The metabolic model containing the reactions.
    phenotype_data : Dict[str, Dict[str, str]]
        Phenotype observations and predictions for each condition.
    media_conditions : Dict[str, Dict[str, float]]
        Media compositions. Format: {condition: {exchange_rxn: lower_bound}}
        Negative values indicate allowed uptake.
    condition_function : callable
        Function that takes phenotype_data and returns a list of condition names.
        Typically get_growth_conditions or get_nogrowth_conditions.
    
    Returns
    -------
    pd.DataFrame
        DataFrame with reactions as rows and conditions as columns.
        Values are lower bounds from media_conditions, defaulting to 0
        for reactions not specified in the media.
    
    Examples
    --------
    >>> lb = get_lower_bound_for_conditions(model, phenotypes, media, get_growth_conditions)
    >>> lb['glucose']['EX_glc']
    -10.0
    >>> lb['glucose']['BIOMASS']
    0.0
    """
    return pd.DataFrame({
        condition: pd.Series(
            {
                rxn.id: (
                    media_conditions[condition][rxn.id]
                    if rxn.id in media_conditions[condition]
                    else 0
                )
                for rxn in model.reactions
            }
        )
        for condition in condition_function(phenotype_data)
    })

def get_lower_bound_for_nogrowth_conditions(
    model: cobra.Model,
    phenotype_data: Dict[str, Dict[str, str]],
    media_conditions: Dict[str, Dict[str, float]]
) -> pd.DataFrame:
    """
    Get lower bound matrix for nogrowth conditions.
    
    Convenience wrapper around get_lower_bound_for_conditions that automatically
    filters to conditions where observed="no_growth" and predicted="growth".
    
    Parameters
    ----------
    model : cobra.Model
        The metabolic model containing the reactions.
    phenotype_data : Dict[str, Dict[str, str]]
        Phenotype observations and predictions for each condition.
    media_conditions : Dict[str, Dict[str, float]]
        Media compositions with negative values indicating allowed uptake.
    
    Returns
    -------
    pd.DataFrame
        DataFrame with reactions as rows and nogrowth conditions as columns.
        Values are lower bounds from media_conditions.
    
    See Also
    --------
    get_lower_bound_for_conditions : General function for any condition set
    get_nogrowth_conditions : Function used to filter conditions
    """
    return get_lower_bound_for_conditions(model, phenotype_data, media_conditions, get_nogrowth_conditions)
    

def get_lower_bound_for_growth_conditions(
    model: cobra.Model,
    phenotype_data: Dict[str, Dict[str, str]],
    media_conditions: Dict[str, Dict[str, float]],
    atp_maintenance_rxn: str,
    atp_maintenance_lower_bound: float
) -> pd.DataFrame:
    """
    Get lower bound matrix for growth conditions with ATP maintenance.
    
    Like get_lower_bound_for_nogrowth_conditions but also sets a minimum ATP
    maintenance requirement if the ATP maintenance reaction exists.
    
    Parameters
    ----------
    model : cobra.Model
        The metabolic model containing the reactions.
    phenotype_data : Dict[str, Dict[str, str]]
        Phenotype observations and predictions for each condition.
    media_conditions : Dict[str, Dict[str, float]]
        Media compositions with negative values indicating allowed uptake.
    atp_maintenance_rxn : str
        Reaction ID for ATP maintenance (e.g., "ATPM").
    atp_maintenance_lower_bound : float
        Minimum ATP maintenance flux required in growth conditions.
    
    Returns
    -------
    pd.DataFrame
        DataFrame with reactions as rows and growth conditions as columns.
        The ATP maintenance reaction (if present) has its lower bound set
        to atp_maintenance_lower_bound in all growth conditions.
    
    Notes
    -----
    ATP maintenance represents the energy requirement for cellular processes
    beyond growth. Setting a minimum ensures the model maintains realistic
    energetic constraints.
    
    See Also
    --------
    get_lower_bound_for_conditions : General function for any condition set
    get_growth_conditions : Function used to filter conditions
    """
    lower_bound_growth =  get_lower_bound_for_conditions(model, phenotype_data, media_conditions, get_growth_conditions)
    for growth_condition in get_growth_conditions(phenotype_data):
        if atp_maintenance_rxn in lower_bound_growth.index:
            lower_bound_growth.loc[atp_maintenance_rxn] = atp_maintenance_lower_bound
    return lower_bound_growth

def get_upper_bound_for_conditions(
    model: cobra.Model,
    phenotype_data: Dict[str, Dict[str, str]],
    media_conditions: Dict[str, Dict[str, float]],
    condition_function: Callable[[Dict[str, Dict[str, str]]], List[str]],
    growth_limit: float
) -> pd.DataFrame:
    """
    Construct upper bound matrix for specified conditions.
    
    Creates a DataFrame where reactions not in the media have their upper bound
    set to growth_limit (allowing internal fluxes), while reactions in the media
    (exchange reactions) have upper bound 0 (preventing secretion).
    
    Parameters
    ----------
    model : cobra.Model
        The metabolic model containing the reactions.
    phenotype_data : Dict[str, Dict[str, str]]
        Phenotype observations and predictions for each condition.
    media_conditions : Dict[str, Dict[str, float]]
        Media compositions. Reactions present in media get upper bound 0.
    condition_function : Callable[[Dict[str, Dict[str, str]]], List[str]]
        Function that takes phenotype_data and returns condition names.
    growth_limit : float
        Upper bound for reactions not in the media (internal reactions).
        For nogrowth conditions, typically maximum_nogrowth.
        For growth conditions, typically minimum_growth.
    
    Returns
    -------
    pd.DataFrame
        DataFrame with reactions as rows and conditions as columns.
        Exchange reactions in media have upper bound 0 (no secretion).
        Other reactions have upper bound growth_limit (allow metabolism).
    
    Notes
    -----
    The asymmetric bounds (negative lower, zero upper) on exchange reactions
    enforce the biological constraint that cells can uptake but not secrete
    nutrients from a defined medium.
    """
    return pd.DataFrame({
        condition: pd.Series(
            {
                rxn.id: (
                    growth_limit
                    if rxn.id not in media_conditions[condition]
                    else 0
                )
                for rxn in model.reactions
            }
        )
        for condition in condition_function(phenotype_data)
    })

def get_upper_bound_for_nogrowth_conditions(
    model: cobra.Model,
    phenotype_data: Dict[str, Dict[str, str]],
    media_conditions: Dict[str, Dict[str, float]],
    maximum_nogrowth: float
) -> pd.DataFrame:
    """
    Get upper bound matrix for nogrowth conditions.
    
    Convenience wrapper that sets upper bounds to maximum_nogrowth for
    conditions where the model incorrectly predicts growth.
    
    Parameters
    ----------
    model : cobra.Model
        The metabolic model containing the reactions.
    phenotype_data : Dict[str, Dict[str, str]]
        Phenotype observations and predictions for each condition.
    media_conditions : Dict[str, Dict[str, float]]
        Media compositions.
    maximum_nogrowth : float
        Upper bound for internal reactions in nogrowth conditions.
        Should be below the growth threshold.
    
    Returns
    -------
    pd.DataFrame
        DataFrame with upper bounds for nogrowth conditions.
    
    See Also
    --------
    get_upper_bound_for_conditions : General function
    """
    return get_upper_bound_for_conditions(model, phenotype_data, media_conditions, get_nogrowth_conditions, maximum_nogrowth)

def get_upper_bound_for_growth_conditions(
    model: cobra.Model,
    phenotype_data: Dict[str, Dict[str, str]],
    media_conditions: Dict[str, Dict[str, float]],
    minimum_growth: float
) -> pd.DataFrame:
    """
    Get upper bound matrix for growth conditions.
    
    Convenience wrapper that sets upper bounds to minimum_growth for
    conditions where the model correctly predicts growth.
    
    Parameters
    ----------
    model : cobra.Model
        The metabolic model containing the reactions.
    phenotype_data : Dict[str, Dict[str, str]]
        Phenotype observations and predictions for each condition.
    media_conditions : Dict[str, Dict[str, float]]
        Media compositions.
    minimum_growth : float
        Upper bound for internal reactions in growth conditions.
        Should be at or above the growth threshold.
    
    Returns
    -------
    pd.DataFrame
        DataFrame with upper bounds for growth conditions.
    
    See Also
    --------
    get_upper_bound_for_conditions : General function
    """
    return get_upper_bound_for_conditions(model, phenotype_data, media_conditions, get_growth_conditions, minimum_growth)

def get_growth_conditions(phenotype_data: Dict[str, Dict[str, str]]) -> List[str]:
    """
    Extract condition names where both observed and predicted phenotypes are growth.
    
    Identifies conditions where the model correctly predicts growth. These are
    used as positive controls - the CROP algorithm must maintain growth in
    these conditions after removing problematic reactions.
    
    Parameters
    ----------
    phenotype_data : Dict[str, Dict[str, str]]
        Phenotype observations and predictions for each condition. Format:
        {condition_name: {"observed": "growth"|"no_growth", "predicted": "growth"|"no_growth"}}
    
    Returns
    -------
    List[str]
        List of condition names where observed="growth" AND predicted="growth".
        These represent correctly predicted growth conditions.
    
    Examples
    --------
    >>> phenotype_data = {
    ...     "glucose": {"observed": "growth", "predicted": "growth"},
    ...     "lactose": {"observed": "no_growth", "predicted": "growth"},
    ...     "minimal": {"observed": "no_growth", "predicted": "no_growth"}
    ... }
    >>> get_growth_conditions(phenotype_data)
    ['glucose']
    """
    return [media_condition for media_condition, phenotype in phenotype_data.items()
            if phenotype["observed"] == "growth" and phenotype["predicted"] == "growth"]

def get_nogrowth_conditions(phenotype_data: Dict[str, Dict[str, str]]) -> List[str]:
    """
    Extract condition names where observed is no-growth but predicted is growth.
    
    Identifies conditions where the model incorrectly predicts growth. These are
    the problematic conditions that the CROP algorithm attempts to fix by removing
    reactions that enable false-positive growth predictions.
    
    Parameters
    ----------
    phenotype_data : Dict[str, Dict[str, str]]
        Phenotype observations and predictions for each condition. Format:
        {condition_name: {"observed": "growth"|"no_growth", "predicted": "growth"|"no_growth"}}
    
    Returns
    -------
    List[str]
        List of condition names where observed="no_growth" AND predicted="growth".
        These represent incorrectly predicted growth conditions that need fixing.
    
    Examples
    --------
    >>> phenotype_data = {
    ...     "glucose": {"observed": "growth", "predicted": "growth"},
    ...     "lactose": {"observed": "no_growth", "predicted": "growth"},
    ...     "minimal": {"observed": "no_growth", "predicted": "no_growth"}
    ... }
    >>> get_nogrowth_conditions(phenotype_data)
    ['lactose']
    """
    return [media_condition for media_condition, phenotype in phenotype_data.items()
            if phenotype["observed"] == "no_growth" and phenotype["predicted"] == "growth"]

def phenotype_data() -> Dict[str, Dict[str, str]]:
    """
    Example phenotype observation data for testing.
    
    Returns
    -------
    Dict[str, Dict[str, str]]
        Sample phenotype data with three conditions:
        - "glucose": correctly predicted growth
        - "lactose": incorrectly predicted growth (needs fixing)
        - "no_carbon": correctly predicted no growth
    
    Notes
    -----
    This is a fixture function used for testing and demonstrations.
    Real applications should provide actual experimental observations.
    """
    return {
        "glucose": {"observed": "growth", "predicted": "growth"},  # Correct
        "lactose": {
            "observed": "no_growth",
            "predicted": "growth",
        },  # Incorrect - needs fixing
        "no_carbon": {"observed": "no_growth", "predicted": "no_growth"},  # Correct
    }

def media_conditions() -> Dict[str, Dict[str, float]]:
    """
    Example media composition data for testing.
    
    Returns
    -------
    Dict[str, Dict[str, float]]
        Sample media conditions with three scenarios:
        - "glucose": 10 mmol/gDW/h glucose uptake allowed
        - "lactose": 10 mmol/gDW/h lactose uptake allowed
        - "no_carbon": no carbon source available
        
        Negative values indicate allowed uptake (COBRA convention).
    
    Notes
    -----
    This is a fixture function used for testing and demonstrations.
    Real applications should provide actual media compositions.
    """
    return {
        "glucose": {"EX_glc": -10, "EX_lac": 0},
        "lactose": {"EX_glc": 0, "EX_lac": -10},
        "no_carbon": {"EX_glc": 0, "EX_lac": 0},
    }

def nogrowth_clause(
    v_nogrowth: Variable,
    biomass_idx: int,
    lower_bound_nogrowth: np.ndarray,
    upper_bound_nogrowth: np.ndarray,
    stoichiometric_matrix: np.ndarray,
    z: Variable,
    r: Variable,
    m: Variable,
    c: np.ndarray,
    omega: float,
    nogrowth_carbon_source_name: str,
    nogrowth_carbon_source_idx: int,
    maximum_nogrowth: float
) -> List:
    """
    Generate constraints that enforce no-growth behavior in nogrowth conditions.
    
    These constraints form the "inner problem" of the bi-level optimization,
    ensuring that the model cannot achieve significant growth in the nogrowth condition
    after removing the identified reactions.
    
    Mathematical formulation:
    
    .. math::
        v_{biomass} = U_{nogrowth,carbon} \\cdot r_{carbon} \\\\
        S \\cdot v_{nogrowth} = 0 \\\\
        L_{nogrowth,i} \\cdot z_i \\leq v_i \\leq U_{nogrowth,i} \\cdot z_i \\\\
        S^T m + c = r \\\\
        r_i \\leq \\Omega \\cdot (1 - z_i) \\\\
        r_{carbon} \\leq 0 \\\\
        v_{biomass} \\leq \\text{maximum nogrowth}
    
    Parameters
    ----------
    v_nogrowth : Variable
        CVXPY variable representing reaction fluxes in the nogrowth condition.
        Corresponds to $v_{nogrowth}$ in the formulation.
        Shape: (n_reactions,)
    biomass_idx : int
        Index of the biomass reaction in the reaction list.
        Used to access $v_{biomass}$ from v_nogrowth.
    lower_bound_nogrowth : np.ndarray
        Lower bounds on reaction fluxes in nogrowth condition.
        Corresponds to $L_{nogrowth}$ in the formulation.
        Typically negative for exchange reactions (uptake allowed), 0 otherwise.
        Shape: (n_reactions,)
    upper_bound_nogrowth : np.ndarray
        Upper bounds on reaction fluxes in nogrowth condition.
        Corresponds to $U_{nogrowth}$ in the formulation.
        Set to maximum_nogrowth for reactions not in media, 0 for absent nutrients.
        Shape: (n_reactions,)
    stoichiometric_matrix : np.ndarray
        Stoichiometric matrix $S$ of the metabolic network.
        S[i,j] is the coefficient of metabolite i in reaction j.
        Shape: (n_metabolites, n_reactions)
    z : Variable
        CVXPY binary variable indicating reaction inclusion (1) or removal (0).
        Corresponds to $z$ in the formulation.
        Shape: (n_reactions,)
    r : Variable
        CVXPY variable for dual variables (reduced costs) associated with flux bounds.
        Corresponds to $r$ in the formulation.
        Interpretation: r_i > 0 means upper bound is active, r_i < 0 means lower bound is active.
        Shape: (n_reactions,)
    m : Variable
        CVXPY variable for dual variables associated with mass balance constraints.
        Corresponds to $m$ in the formulation (metabolite shadow prices).
        Shape: (n_metabolites,)
    c : np.ndarray
        Objective function coefficients (typically 1 for biomass, 0 elsewhere).
        Corresponds to $c$ in the formulation.
        Shape: (n_reactions,)
    omega : float
        Large positive constant (typically 1000) serving as upper bound on |r_i|.
        Corresponds to $\\Omega$ in the formulation.
        The constraint $r_i \\leq \\Omega(1-z_i)$ forces r_i = 0 when z_i = 1.
    nogrowth_carbon_source_name : str
        Name/ID of the primary carbon source exchange reaction in nogrowth condition.
        Used to index into lower_bound_nogrowth array.
    nogrowth_carbon_source_idx : int
        Integer index of the carbon source reaction in the reaction list.
        Used to access elements of r and v_nogrowth.
    maximum_nogrowth : float
        Maximum allowed biomass flux in nogrowth condition.
        Corresponds to $\\text{maximum nogrowth}$ in the formulation.
        Typically set below the threshold defining growth (e.g., 1.0).
    
    Returns
    -------
    List
        List of CVXPY constraints that enforce the nogrowth condition behavior.
        These constraints ensure the model cannot grow significantly in the
        nogrowth condition after removing reactions (z_i = 0).
    
    Notes
    -----
    The key insight is linking biomass flux to the dual variable of the carbon source:
    $v_{biomass} = U_{nogrowth,carbon} \\cdot r_{carbon}$
    
    This relationship, combined with $r_{carbon} \\leq 0$ and the dual feasibility
    constraints, ensures that when the carbon source cannot support growth, neither
    can any other pathway.
    
    The constraint $r_i \\leq \\Omega(1-z_i)$ is crucial:
    - When z_i = 1 (reaction included): r_i ≤ 0, so upper bound cannot be tight
    - When z_i = 0 (reaction removed): r_i ≤ Ω, essentially unbounded
    
    This prevents reactions from being at their upper bound when included, which
    would indicate they're constraining growth.
    """

    return [
        v_nogrowth[biomass_idx] == lower_bound_nogrowth[nogrowth_carbon_source_name] * r[nogrowth_carbon_source_idx],  # & v_{biomass} = U_{nogrowth,lactose}r_{lactose} \\
        stoichiometric_matrix @ v_nogrowth == 0,  # & Sv_{nogrowth}= 0 & \text{inner problem} \\
        diag(lower_bound_nogrowth) @ z <= v_nogrowth,
        v_nogrowth <= diag(upper_bound_nogrowth) @ z,  # & 0\leq v_i\leq U_{nogrowth,i}\cdot z_i & \text{every carbon source $i$ not in the nogrowth media has $U_{nogrowth,i} = 0$. \\ Every carbon source $j$ in the nogrowth media has $U_{nogrowth,j} > 0$} \\
        stoichiometric_matrix.T @ m + c == r,  # & S^Tm +c = r \\
        r <= omega * (1 - z),  # & r_i\leq\Omega_i\cdot(1-z_i) & \text{for $i\neq$ lactose. This constraint ensures that lactose uptake is the only tight constraint in the model. Every other reaction with a tight constraint cannot be part of the model.}\\
        r[nogrowth_carbon_source_idx] <= 0,  # & r_{lactose_uptake} \leq 0 \\
        v_nogrowth[biomass_idx] <= maximum_nogrowth,  # v_{biomass} \leq\text{minimal growth} \\
    ]

def growth_clause(
    v_growth: Variable,
    biomass_idx: int,
    lower_bound_growth: np.ndarray,
    upper_bound_growth: np.ndarray,
    stoichiometric_matrix: np.ndarray,
    z: Variable,
    minimum_growth: float
) -> List:
    """
    Generate constraints that enforce growth behavior in growth conditions.
    
    These constraints ensure that the model can still achieve sufficient growth
    in the designated growth condition after removing the identified reactions.
    
    Mathematical formulation:
    
    .. math::
        S \\cdot w_{growth} = 0 \\\\
        L_{growth,i} \\cdot z_i \\leq w_i \\leq U_{growth,i} \\cdot z_i \\\\
        w_{biomass} \\geq \\text{minimum growth}
    
    Parameters
    ----------
    v_growth : Variable
        CVXPY variable representing reaction fluxes in the growth condition.
        Corresponds to $w_{growth}$ (or $w$) in the formulation.
        Shape: (n_reactions,)
    biomass_idx : int
        Index of the biomass reaction in the reaction list.
        Used to access $w_{biomass}$ from v_growth.
    lower_bound_growth : np.ndarray
        Lower bounds on reaction fluxes in growth condition.
        Corresponds to $L_{growth}$ in the formulation.
        Typically negative for exchange reactions providing nutrients, 0 otherwise.
        May include ATP maintenance requirements.
        Shape: (n_reactions,)
    upper_bound_growth : np.ndarray
        Upper bounds on reaction fluxes in growth condition.
        Corresponds to $U_{growth}$ in the formulation.
        Set to minimum_growth for reactions not supplying the carbon source, 0 for absent nutrients.
        Shape: (n_reactions,)
    stoichiometric_matrix : np.ndarray
        Stoichiometric matrix $S$ of the metabolic network.
        S[i,j] is the coefficient of metabolite i in reaction j.
        Shape: (n_metabolites, n_reactions)
    z : Variable
        CVXPY binary variable indicating reaction inclusion (1) or removal (0).
        Corresponds to $z$ in the formulation.
        When z_i = 0, both bounds on reaction i become 0, effectively removing it.
        Shape: (n_reactions,)
    minimum_growth : float
        Minimum required biomass flux to constitute growth.
        Corresponds to $\\text{minimum growth}$ in the formulation.
        Typically set above the threshold defining growth (e.g., 2.0).
    
    Returns
    -------
    List
        List of CVXPY constraints that enforce the growth condition behavior.
        These constraints ensure the model can still grow after reaction removal.
    
    Notes
    -----
    The growth constraints are simpler than nogrowth constraints because we only need
    to verify that a feasible flux distribution exists that achieves minimum growth.
    
    Key differences from nogrowth_clause:
    - No dual variables (r, m) needed - this is a primal feasibility check
    - Direct constraint on biomass flux: $w_{biomass} \\geq \\text{minimum growth}$
    - Bounds are gated by z: when z_i = 0, both bounds become 0 (reaction removed)
    
    The constraint $L_{growth,i} \\cdot z_i \\leq w_i \\leq U_{growth,i} \\cdot z_i$
    effectively implements:
    - When z_i = 1 (included): normal flux bounds apply
    - When z_i = 0 (removed): w_i must be 0 (no flux through removed reaction)
    
    Examples
    --------
    For a reaction with bounds [-10, 1000]:
    - If z = 1: -10 ≤ w ≤ 1000 (normal operation)
    - If z = 0: 0 ≤ w ≤ 0 (forced to zero, effectively removed)
    """
    return [
        stoichiometric_matrix @ v_growth == 0,  # Sw= 0 \\
        diag(lower_bound_growth) @ z <= v_growth,  # & L_{growth,i} \cdot z_i \leq w_i \\
        v_growth <= diag(upper_bound_growth) @ z,  # 0\leq w_i\leq U_{growth,i}\cdot z_i \\
        v_growth[biomass_idx] >= minimum_growth  # w_{biomass} \geq \text{minimal growth} 
    ]