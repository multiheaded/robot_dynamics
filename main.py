import numpy as np

from sympy import Symbol, Function, pretty, latex, diff, Eq, dotprint, Add, simplify, trigsimp, symbols, shape, pi
from sympy.abc import theta, alpha, d, a, t
from sympy.core.function import Subs
from sympy.matrices import Matrix, BlockMatrix
from sympy.solvers.solveset import linear_eq_to_matrix

from kinematics import SymbolicJointValues, KinematicChain
from transforms import DHSet, Jacobian

#, eye, zeros, ones, diag, GramSchmidt
#from sympy.solvers import solve

from sympy.interactive.printing import init_printing
init_printing(use_unicode=True, wrap_line=False, pretty_print=True, latex_mode=True)


#first = DHSet(theta = 0, d=d, a=a, alpha=alpha)
#DHMat = transforms.DH(theta, d, a, alpha)
#DHMat = DHMat.subs(theta,first.theta)

q = SymbolicJointValues(3)

h1 = 0.5
h2 = 0.2
h3 = 0.5
l  = 0.2

DHParameters = [  DHSet(q[0]+pi/2, h1, 0, pi/2),
                  DHSet(q[1], l, 0, -pi/2),
                  DHSet(-pi/2, q[2]+h2+h3, 0, 0) ]

STANF = KinematicChain(DHParameters, q)

print(latex(STANF.forwardKinematics()))
print(latex(STANF.SymbolicJacobian))