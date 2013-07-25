
from __future__ import division
from __future__ import print_function

class Determinant:
    """
        Determinant class - set up for later structurization
    """

    def __init__(self, group, value):
        """
        Contructer of determinant object - simple, but helps in creating structure!
        """
        self.group = group
        self.label = group.label
        self.value = value

        return

    def add(self, value):
        """
        adding a value to determinant
        """
        self.value += value

        return
    
    def __str__(self):
        return '%s: %8.2f'%(self.label,self.value)
