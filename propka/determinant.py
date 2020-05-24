"""Holds the Determinant class

TODO - it is confusing to have both `determinant.py` and `determinants.py`.
Should these be merged?
"""


class Determinant:
    """Determinant class.

    Appears to be a container for storing information and values about
    groups that interact to influence titration states.

    TODO - figure out what this class does.
    """

    def __init__(self, group, value):
        """Initialize the object.

        Args:
            group:  group associated with Determinant object
            value:  value to assign to group
        """
        self.group = group
        self.label = group.label
        self.value = value

    def add(self, value):
        """Increment determinant value.

        Args:
            value:  value to add to determinant
        """
        self.value += value

    def __str__(self):
        return '%s: %8.2f' % (self.label, self.value)
