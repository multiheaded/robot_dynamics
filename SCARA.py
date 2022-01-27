from sympy import Symbol, pi
from sympy.abc import g
from sympy.matrices import Matrix

from kinematics import SymbolicJointValues, Kinematics
from dynamics import Dynamics, CylindricBody
from simplifier import Simplifier
from transforms import DHSet
from symbolic import SymbolicConstantVector

def computeSCARA(simp = Simplifier("None")):
    q = SymbolicJointValues(2)
    l1 = Symbol('l_1', real=True)
    m1 = Symbol('m_1', real=True)
    r1 = Symbol('r_1', real=True)
    l2 = Symbol('l_2', real=True)
    m2 = Symbol('m_2', real=True)
    r2 = Symbol('r_2', real=True)
    DHParameters = [DHSet(q[0], 0, l1, 0),
                    DHSet(q[1], 0, l2, 0)]
    
    SCARA = Kinematics(DHParameters, q, simp)
    
    Bodies = [ CylindricBody(m1, l1, r1), CylindricBody(m2, l2, r2) ]
    Gravitation = Matrix([0,g,0])
    
    DynamicSCARA = Dynamics(SCARA, Bodies, Gravitation, simp)
    
    return SCARA, DynamicSCARA