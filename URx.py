from sympy import Symbol, pi
from sympy.abc import g
from sympy.matrices import Matrix

from kinematics import SymbolicJointValues, Kinematics
from dynamics import Dynamics, Body
from simplifier import Simplifier
from transforms import DHSet
from symbolic import SymbolicConstantVector

def computeURxRobot(simp = Simplifier("None")):
    q = SymbolicJointValues(6)
    a2 = Symbol('a2', real=True)
    a3 = Symbol('a3', real=True)
    d1 = Symbol('d1', real=True)
    d4 = Symbol('d4', real=True)
    d5 = Symbol('d5', real=True)
    d6 = Symbol('d6', real=True)

    m = SymbolicConstantVector("m", 6)
    x = SymbolicConstantVector("x", 6)
    y = SymbolicConstantVector("y", 6)
    z = SymbolicConstantVector("z", 6)

    # DH Parameters from https://www.universal-robots.com/articles/ur/application-installation/dh-parameters-for-calculations-of-kinematics-and-dynamics/
    #                       theta,  d,  a,   alpha
    DHParameters = [DHSet( q[0], d1,  0,  pi/2.0),
                    DHSet( q[1],  0, a2,       0),
                    DHSet( q[2],  0, a3,       0),
                    DHSet( q[3], d4,  0,  pi/2.0),
                    DHSet( q[4], d5,  0, -pi/2.0),
                    DHSet( q[5], d6,  0,       0)]

    URx = Kinematics(DHParameters, q, simp)

    Bodies = [ Body(CoM=Matrix([x[i], y[i], z[i]]), Mass=m[i]) for i in range(6)]
    Gravitation = Matrix([0,0, g])

    DynamicURx = Dynamics(URx, Bodies, Gravitation, simp)
    
    return URx, DynamicURx