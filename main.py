import numpy as np
import time
from sympy import Symbol, Function, pretty, latex, diff, Eq, dotprint, Add, simplify, trigsimp, symbols, shape, pi, expand
from sympy.abc import theta, alpha, d, a, t, g
from sympy.core.function import Subs
from sympy.matrices import Matrix, BlockMatrix
from sympy.solvers.solveset import linear_eq_to_matrix

from kinematics import SymbolicJointValues, Kinematics
from dynamics import CylindricBody, Dynamics, Body
from simplifier import Simplifier
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

"""q = SymbolicJointValues(2)
l1 = Symbol('l1', real=True)
m1 = Symbol('m1', real=True)
r1 = Symbol('r1', real=True)
l2 = Symbol('l2', real=True)
m2 = Symbol('m2', real=True)
r2 = Symbol('r2', real=True)
DHParameters = [  DHSet(q[0], 0, l1, 0),
                  DHSet(q[1], 0, l2, 0) ]

tic = time.perf_counter()
SCARA = Kinematics(DHParameters, q)
toc = time.perf_counter()
print(f"Kinematics done. Took {toc - tic:0.4f} seconds")
print(SCARA.forward(1,2))

#Bodies = [  Body(CoM=Matrix([-l1/2, 0, 0])),
#            Body(CoM=Matrix([-l2/2, 0, 0])) ]
Bodies = [ CylindricBody(m1, l1, r1), CylindricBody(m2, l2, r2) ]
Gravitation = Matrix([0,g,0])

tic = time.perf_counter()
DynamicSCARA = Dynamics(SCARA, Bodies, Gravitation)
toc = time.perf_counter()
print(f"Dynamics done. Took {toc - tic:0.4f} seconds")

tic = time.perf_counter()
T_inv = simplify(DynamicSCARA.Theta.inv())
toc = time.perf_counter()
print(f"Forward dynamics done. Took {toc - tic:0.4f} seconds")

tic = time.perf_counter()
L = simplify(DynamicSCARA.Theta.cholesky(False))
toc = time.perf_counter()
print(f"Decomposition done. Took {toc - tic:0.4f} seconds")

tic = time.perf_counter()
testmat = simplify(L*L.transpose())
toc = time.perf_counter()
print(f"LL* done. Took {toc - tic:0.4f} seconds")

print(latex(simplify(DynamicSCARA.Theta)))

print("")

print(latex(simplify(trigsimp(expand(testmat[1,1])))))


#print(latex(SCARA.forwardKinematics()))
#print(latex(SCARA.SymbolicJacobian)) """

def SymbolicConstantVector(name, N):
  return Matrix([Symbol('%s%d'%(name, i+1), real=True) for i in range(N)])

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

simp = Simplifier("None")

# DH Parameters from https://www.universal-robots.com/articles/ur/application-installation/dh-parameters-for-calculations-of-kinematics-and-dynamics/
#                       theta,  d,  a,   alpha
DHParameters = [  DHSet( q[0], d1,  0,  pi/2.0),
                  DHSet( q[1],  0, a2,       0),
                  DHSet( q[2],  0, a3,       0),
                  DHSet( q[3], d4,  0,  pi/2.0),
                  DHSet( q[4], d5,  0, -pi/2.0),
                  DHSet( q[5], d6,  0,       0) ]

print("Computing kinematics")
tic = time.perf_counter()
SCARA = Kinematics(DHParameters, q)
toc = time.perf_counter()
print(f"Kinematics done. Took {toc - tic:0.4f} seconds")

Bodies = [ Body(CoM=Matrix([x[i], y[i], z[i]]), Mass=m[i]) for i in range(6)]
Gravitation = Matrix([0,g,0])

tic = time.perf_counter()
DynamicSCARA = Dynamics(SCARA, Bodies, Gravitation)
toc = time.perf_counter()
print(f"Dynamics done. Took {toc - tic:0.4f} seconds")

print("Inverting Theta")
tic = time.perf_counter()
T_inv = simplify(DynamicSCARA.Theta.inv())
toc = time.perf_counter()
print(f"Forward dynamics done. Took {toc - tic:0.4f} seconds")

tic = time.perf_counter()
L = simplify(DynamicSCARA.Theta.cholesky(False))
toc = time.perf_counter()
print(f"Decomposition done. Took {toc - tic:0.4f} seconds")

tic = time.perf_counter()
testmat = simplify(L*L.transpose())
toc = time.perf_counter()
print(f"LL* done. Took {toc - tic:0.4f} seconds")

print("")
print("")
print("")
print("")
print("T")
print(latex(DynamicSCARA.Theta))

print("")
print("")
print("")
print("")
print("CORZEN")
print(latex(DynamicSCARA.CORZEN))

print("")
print("")
print("")
print("")
print("GRAV")
print(latex(DynamicSCARA.Grav))

print("")
print("")
print("")
print("")
print("T_inv")
print(latex(T_inv))

print("")
print("")
print("")
print("")
print("L")
print(latex(L))

print("")
print("")
print("")
print("")
print("testmat (should be equal to T)")
print(latex(testmat))

