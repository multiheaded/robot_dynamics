import numpy as np

from sympy import Symbol, Function, pretty, latex, diff, Eq, dotprint, Add, simplify, trigsimp, symbols, shape, pi
from sympy.abc import theta, alpha, d, a, t, g
from sympy.core.function import Subs
from sympy.matrices import Matrix, BlockMatrix
from sympy.solvers.solveset import linear_eq_to_matrix

from kinematics import SymbolicJointValues, Kinematics
from dynamics import CylindricBody, Dynamics, Body
from transforms import DHSet, Jacobian

#, eye, zeros, ones, diag, GramSchmidt
#from sympy.solvers import solve

from sympy.interactive.printing import init_printing
init_printing(use_unicode=True, wrap_line=False, pretty_print=True, latex_mode=True)


#first = DHSet(theta = 0, d=d, a=a, alpha=alpha)
#DHMat = transforms.DH(theta, d, a, alpha)
#DHMat = DHMat.subs(theta,first.theta)

""" q = SymbolicJointValues(3)
h1 = 0.5
h2 = 0.2
h3 = 0.5
l  = 0.2
DHParameters = [  DHSet(q[0]+pi/2, h1, 0, pi/2),
                  DHSet(q[1], l, 0, -pi/2),
                  DHSet(-pi/2, q[2]+h2+h3, 0, 0) ]
STANF = KinematicChain(DHParameters, q)
print(latex(STANF.forwardKinematics()))
print(latex(STANF.SymbolicJacobian)) """

q = SymbolicJointValues(2)
l1 = Symbol('l1', real=True)
m1 = Symbol('m1', real=True)
r1 = Symbol('r1', real=True)
l2 = Symbol('l2', real=True)
m2 = Symbol('m2', real=True)
r2 = Symbol('r2', real=True)
DHParameters = [  DHSet(q[0], 0, l1, 0),
                  DHSet(q[1], 0, l2, 0) ]
SCARA = Kinematics(DHParameters, q)
#Bodies = [  Body(CoM=Matrix([-l1/2, 0, 0])),
#            Body(CoM=Matrix([-l2/2, 0, 0])) ]
Bodies = [ CylindricBody(m1, l1, r1), CylindricBody(m2, l2, r2) ]
Gravitation = Matrix([0,g,0])
DynamicSCARA = Dynamics(SCARA, Bodies, Gravitation)

print(latex(simplify(DynamicSCARA.Grav)))


#print(latex(SCARA.forwardKinematics()))
#print(latex(SCARA.SymbolicJacobian))