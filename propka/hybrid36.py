import string

_hybrid36_upper_chars = set(string.ascii_uppercase)
_hybrid36_lower_chars = set(string.ascii_lowercase)
_hybrid36_digits = set(string.digits)
_hybrid36_upper_set = _hybrid36_upper_chars | _hybrid36_digits
_hybrid36_lower_set = _hybrid36_lower_chars | _hybrid36_digits

def decode(input_string):
    """
    Convert an input string of a number in hybrid-36 format to an integer.

    """
    value_error_message = "invalid literal for hybrid-36 conversion: '%s'"

    original_input_string = input_string
    input_string = input_string.strip()

    # Manually handle negative sign.
    if input_string.startswith("-"):
        sign = -1
        input_string = input_string[1:]
    else:
        sign = 1

    if not len(input_string):
        raise ValueError(value_error_message % input_string)

    # See http://cci.lbl.gov/hybrid_36/ for documentation on the format.

    num_chars = len(input_string)
    first_char = input_string[0]

    if first_char in _hybrid36_digits:
        return sign * int(input_string)
    elif first_char in _hybrid36_upper_chars:
        reference = - (10 * 36 ** (num_chars - 1) - 10 ** num_chars)
        _hybrid36_set = _hybrid36_upper_set
    elif first_char in _hybrid36_lower_chars:
        reference = (16 * 36 ** (num_chars - 1) + 10 ** num_chars)
        _hybrid36_set = _hybrid36_lower_set
    else:
        raise ValueError(value_error_message % original_input_string)

    # Check the validity of the input string: ASCII characters should be
    # either all uppercase or all lowercase.
    for c in input_string[1:]:
        if c not in _hybrid36_set:
            raise ValueError(value_error_message % original_input_string)

    # Convert with the int function.
    return sign * (int(input_string, 36) + reference)
