"""This script generates an example metabolic model in SBML format using the CobraPy library.
It creates a simple model with a few metabolites and reactions, writes it to an SBML file,
and validates the model."""
from pprint import pprint

import cobra
import pandas as pd
from cobra.io import validate_sbml_model, write_sbml_model
from cobra.util.array import create_stoichiometric_matrix

pd.options.display.float_format = "{:.0f}".format

mets = ["A", "B", "C", "D", "E", "F", "A_ext", "C_ext", "D_ext", "E_ext", "F_ext"]

internal_mets = [met for met in mets if "ext" not in met]
rxns = ["R_{}".format(i) for i in range(1, 11)]

data = {
    "R_1": pd.Series({"A": 1, "A_ext": -1}),
    "R_2": pd.Series({"A": -1, "B": 1}),
    "R_3": pd.Series({"A": -1, "C": 1}),
    "R_4": pd.Series({"B": -1, "D": 2, "E": -1}),
    "R_5": pd.Series({"E": 1, "E_ext": -1}),
    "R_6": pd.Series({"B": -2, "C": 1, "F": 1}),
    "R_7": pd.Series({"C": -1, "D": 1}),
    "R_8": pd.Series({"D": -1, "D_ext": 1}),
    "R_9": pd.Series({"F": -1, "F_ext": 1}),
    "R_10": pd.Series({"C": -1, "C_ext": 1}),
}

model = cobra.Model("example_model")

A = cobra.Metabolite("A_c", compartment="c")
A_e = cobra.Metabolite("A_e", compartment="e")
B = cobra.Metabolite("B_c", compartment="c")
C = cobra.Metabolite("C_c", compartment="c")
C_e = cobra.Metabolite("C_e", compartment="e")
D = cobra.Metabolite("D_c", compartment="c")
D_e = cobra.Metabolite("D_e", compartment="e")
E = cobra.Metabolite("E_c", compartment="c")
E_e = cobra.Metabolite("E_e", compartment="e")
F = cobra.Metabolite("F_c", compartment="c")
F_e = cobra.Metabolite("F_e", compartment="e")

model.add_metabolites([A, B, C, D, E, F, A_e, C_e, D_e, E_e, F_e])

model.add_boundary(model.metabolites.get_by_id("A_e"), type="exchange")
model.add_boundary(model.metabolites.get_by_id("C_e"), type="exchange")
model.add_boundary(model.metabolites.get_by_id("D_e"), type="exchange")
model.add_boundary(model.metabolites.get_by_id("E_e"), type="exchange")
model.add_boundary(model.metabolites.get_by_id("F_e"), type="exchange")

ex1 = cobra.Reaction("EX_A")
ex1.add_metabolites(
    {
        A_e: -1.0,
        A: 1.0,
    }
)

ex2 = cobra.Reaction("EX_C")
ex2.add_metabolites(
    {
        C: -1.0,
        C_e: 1.0,
    }
)

ex3 = cobra.Reaction("EX_D")
ex3.add_metabolites({D: -1.0, D_e: 1.0})

ex4 = cobra.Reaction("EX_E")
ex4.add_metabolites(
    {
        E_e: -1.0,
        E: 1.0,
    }
)

ex5 = cobra.Reaction("SINK_F")
ex5.add_metabolites(
    {
        F: -1.0,
        F_e: 1.0,
    }
)

r1 = cobra.Reaction("R_A_to_B")
r1.add_metabolites(
    {
        A: -1.0,
        B: 1.0,
    }
)

r2 = cobra.Reaction("R_BE_to_D")
r2.add_metabolites(
    {
        B: -1.0,
        E: -1.0,
        D: 2.0,
    }
)

r3 = cobra.Reaction("R_A_to_C")
r3.add_metabolites(
    {
        A: -1.0,
        C: 1.0,
    }
)

r4 = cobra.Reaction("R_BC_to_F")
r4.add_metabolites({B: -2.0, C: 1.0, F: 1})

r5 = cobra.Reaction("R_C_to_D")
r5.add_metabolites(
    {
        C: -1.0,
        D: 1.0,
    }
)

model.add_reactions([ex1, ex2, ex3, ex4, ex5, r1, r2, r3, r4, r5])
model.objective = model.reactions.get_by_id("EX_D")

write_sbml_model(model, filename="example_model.xml")
report = validate_sbml_model(filename="example_model.xml")
pprint(report)


stoich_matrix = create_stoichiometric_matrix(model, "DataFrame")

fullS = pd.DataFrame(data, columns=rxns, index=mets, dtype="int64").fillna(0)

order = [
    "EX_A",
    "R_A_to_B",
    "R_A_to_C",
    "R_BE_to_D",
    "EX_E",
    "R_BC_to_F",
    "R_C_to_D",
    "EX_D",
    "SINK_F",
    "EX_C",
]

stoich_matrix = stoich_matrix[order]

(stoich_matrix.values == fullS.values).all()
print(stoich_matrix)
