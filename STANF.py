from sympy import pi

from kinematics import SymbolicJointValues, Kinematics
from simplifier import Simplifier
from transforms import DHSet

def computeSTANF(simp = Simplifier("None")):
    q = SymbolicJointValues(3)
    h1 = 0.5
    h2 = 0.2
    h3 = 0.5
    l  = 0.2
    DHParameters = [DHSet(q[0]+pi/2,         h1, 0,  pi/2),
                    DHSet(q[1]     ,          l, 0, -pi/2),
                    DHSet(    -pi/2, q[2]+h2+h3, 0,     0)]
    STANF = Kinematics(DHParameters, q, simp)

    return STANF, None