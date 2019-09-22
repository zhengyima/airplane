#       __DUBINSAIRPLANEFUNCTIONS__
#       A small set of helping functions
#
#       Authors: 
#       Kostas Alexis (konstantinos.alexis@mavt.ethz.ch)

import builtins
    
def my_max(a, b):
    # just a simple max function between 2 values
    if a >=b:
        return a
    else:
        return b

def my_min(a, b=0, nargout=0):
    # just a simple min function between 2 values
    if a <= b:
        return a
    else:
        return b
