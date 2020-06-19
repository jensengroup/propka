"""
Hybrid36 PDB-like file format
=============================

`hybrid36`_ is an alternative PDB format that can encode larger atom
numbers. This module provides the :func:`decode` functon to parse the
atom numbers in hybrid36 "PDB" files.

.. _hybrid36: http://cci.lbl.gov/hybrid_36/

"""
import string


_HYBRID36_UPPER_CHARS = set(string.ascii_uppercase)
_HYBRID36_LOWER_CHARS = set(string.ascii_lowercase)
_HYBRID36_DIGITS = set(string.digits)
_HYBRID36_UPPER_SET = _HYBRID36_UPPER_CHARS | _HYBRID36_DIGITS
_HYBRID36_LOWER_SET = _HYBRID36_LOWER_CHARS | _HYBRID36_DIGITS


def decode(input_string):
    """Convert an input string of a number in hybrid-36 format to an integer.

    Args:
        input_string:  input string
    Returns:
        integer
    """
    value_error_message = "invalid literal for hybrid-36 conversion: '{0:s}'"

    original_input_string = input_string
    input_string = input_string.strip()

    # Manually handle negative sign.
    if input_string.startswith("-"):
        sign = -1
        input_string = input_string[1:]
    else:
        sign = 1

    if len(input_string) == 0:
        raise ValueError(value_error_message.format(input_string))

    # See http://cci.lbl.gov/hybrid_36/ for documentation on the format.

    num_chars = len(input_string)
    first_char = input_string[0]

    if first_char in _HYBRID36_DIGITS:
        return sign * int(input_string)
    elif first_char in _HYBRID36_UPPER_CHARS:
        reference = - (10 * 36 ** (num_chars - 1) - 10 ** num_chars)
        _hybrid36_set = _HYBRID36_UPPER_SET
    elif first_char in _HYBRID36_LOWER_CHARS:
        reference = (16 * 36 ** (num_chars - 1) + 10 ** num_chars)
        _hybrid36_set = _HYBRID36_LOWER_SET
    else:
        raise ValueError(value_error_message.format(original_input_string))

    # Check the validity of the input string: ASCII characters should be
    # either all uppercase or all lowercase.
    for char in input_string[1:]:
        if char not in _hybrid36_set:
            raise ValueError(value_error_message.format(original_input_string))

    # Convert with the int function.
    return sign * (int(input_string, 36) + reference)
