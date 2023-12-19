"""
Determinant
===========

Provides the :class:`Determinant` class.

.. TODO::

   It is confusing to have both `determinant.py` and `determinants.py`.
   Should these be merged?

.. SeeAlso::
   - :mod:`propka.determinants`
   - :mod:`propka.iterative`

"""

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from propka.group import Group


class Determinant:
    """Determinant class.

    Appears to be a container for storing information and values about
    groups that interact to influence titration states.

    TODO - figure out what this class does.
    """

    def __init__(self, group: "Group", value: float):
        """Initialize the object.

        Args:
            group:  group associated with Determinant object
            value:  value to assign to group
        """
        self.group = group
        self.label = group.label
        self.value = value

    def add(self, value: float):
        """Increment determinant value.

        Args:
            value:  value to add to determinant
        """
        self.value += value

    def __str__(self):
        return '{0:s}: {1:8.2f}'.format(self.label, self.value)
