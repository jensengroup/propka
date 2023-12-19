import propka.lib as m

import argparse
import pytest


def test_parse_res_string():
    assert m.parse_res_string("C:123") == ("C", 123, " ")
    assert m.parse_res_string("C:123B") == ("C", 123, "B")
    assert m.parse_res_string("ABC:123x") == ("ABC", 123, "x")
    with pytest.raises(ValueError):
        m.parse_res_string("C:B123")
    with pytest.raises(ValueError):
        m.parse_res_string("123B")
    with pytest.raises(ValueError):
        m.parse_res_string("C:123:B")


def test_parse_res_list():
    assert m.parse_res_list("C:123") == [("C", 123, " ")]
    assert m.parse_res_list("ABC:123,D:4,F:56X") == [
        ("ABC", 123, " "),
        ("D", 4, " "),
        ("F", 56, "X"),
    ]
    with pytest.raises(argparse.ArgumentTypeError):
        m.parse_res_list("C:B123")
