"""Test functionalities of gap filling."""

import pathlib

import cobra.io

path = pathlib.Path(__file__).parent


def load_model():
    return cobra.io.read_sbml_model(path.joinpath("example_model.xml").__str__())


def test_e_is_no_growth():
    # E  nutrient condition is no growth
    model = load_model()
    media = {"EX_E_e": 100}
    model.medium = media
    obj_func = model.slim_optimize()
    assert obj_func == 0.0


def test_a_alone_is_growth():
    # A alone is growth
    model = cobra.io.read_sbml_model(path.joinpath("example_model.xml").__str__())
    media = {"EX_A_e": 100}
    model.medium = media
    obj_func = model.slim_optimize()
    assert obj_func == 100.0


def test_a_and_v3_ko_is_growth():
    # A + v3 knockout is no-growth
    model = cobra.io.read_sbml_model(path.joinpath("example_model.xml").__str__())
    media = {"EX_A_e": 100}
    model.remove_reactions(model.reactions.get_by_id("R_A_to_C"))
    model.medium = media
    obj_func = model.slim_optimize()

    assert obj_func == 50.0


def test_a_and_e_and_v3_ko_is_growth():
    # A + E + v3 knockout is growth
    model = cobra.io.read_sbml_model(path.joinpath("example_model.xml").__str__())
    media = {"EX_A_e": 100, "EX_E_e": 100}
    model.remove_reactions(model.reactions.get_by_id("R_A_to_C"))
    model.medium = media
    obj_func = model.slim_optimize()
    assert obj_func > 50.0
