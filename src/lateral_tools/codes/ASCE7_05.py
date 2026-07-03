"""
A library for carrying out calculations to ASCE7-05

*** INCOMPLETE ***

"This Standard provides requirements for
general structural design and includes
means for determining dead, live, soil,
flood, wind, snow, rain, atmospheric
ice, and earthquake loads, as well as
their combinations, which are suitable
for inclusion in building codes and
other documents."

Previous editions are:
* ASCE 7-02
* ASCE 7-98
* ASCE 7-95
* ASCE 7-93
* ASCE 7-88
"""

import doctest

# F_a and F_v taken from ASCE7-05, tables 11.4-1 & 11.4-2
F_a_headers = (0.25, 0.5, 0.75, 1.0, 1.25, 1.5)
F_v_headers = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
F_a_dict = {
    "A": (0.8, 0.8, 0.8, 0.8, 0.8, 0.8),
    "B": (1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
    "C": (1.2, 1.2, 1.1, 1.0, 1.0, 1.0),
    "D": (1.6, 1.4, 1.2, 1.1, 1.0, 1.0),
    "E": (2.5, 1.7, 1.2, 0.9, 0.9, 0.9),
    "F": (None, None, None, None, None, None),
}
F_v_dict = {
    "A": (0.8, 0.8, 0.8, 0.8, 0.8, 0.8),
    "B": (1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
    "C": (1.7, 1.6, 1.5, 1.4, 1.3, 1.3),
    "D": (2.4, 2.0, 1.8, 1.6, 1.5, 1.5),
    "E": (3.5, 3.2, 2.8, 2.4, 2.4, 2.4),
    "F": (None, None, None, None, None, None),
}


# These are values from an unknown source
# F_a and F_v taken from Section 11.4 of ASCE7 ??
# F_a_headers = (0.25,0.5,0.75,1.0,1.25,1.5)
# F_v_headers = (0.1,0.2,0.3,0.4,0.5,0.6)
# F_a_dict = {'A': (0.8,0.8,0.8,0.8,0.8,0.8),
#            'B': (0.9,0.9,0.9,0.9,0.9,0.9),
#            'C': (1.3,1.3,1.2,1.2,1.2,1.2),
#            'D': (1.6,1.4,1.2,1.1,1.0,1.0),
#            'E': (2.4,1.7,1.3,None,None,None),
#            'F': (None,None,None,None,None,None)}
# F_v_dict = {'A': (0.8,0.8,0.8,0.8,0.8,0.8),
#            'B': (0.8,0.8,0.8,0.8,0.8,0.8),
#            'C': (1.5,1.5,1.5,1.5,1.5,1.4),
#            'D': (2.4,2.2,2.0,1.9,1.8,1.7),
#            'E': (4.2,3.3,2.8,2.4,2.2,2.0),
#            'F': (None,None,None,None,None,None)}


def main():  # noqa
    pass


if __name__ == "__main__":
    doctest.testmod()
