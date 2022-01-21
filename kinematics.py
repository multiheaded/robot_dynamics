from sympy import Function
from sympy.abc import t
from sympy.matrices import Matrix, eye

from simplifier import Simplifier

import time

from transforms import DH, Jacobian

def SymbolicJointValues(N):
  return Matrix([Function('q%d'%(i+1), real=True)(t) for i in range(N)])

class Kinematics:
  DHParameter = None
  CachedPartialTransforms = None
  SymbolicJacobian = None
  JointValueSymbols = None
  simpl = None
  
  def __init__(self, DHParms, jointSyms, simp = Simplifier("None")):
    self.simpl = simp
    self.JointValueSymbols = jointSyms
    self.DHParameter = DHParms
    print(" Filling transformation cache")
    tic = time.perf_counter()
    self.fillTransformationCache()
    toc = time.perf_counter()
    print(f" Filling transformation cache done. Took {toc - tic:0.4f} seconds")
    print(" Computing simplified symbolic Jacobian")
    tic = time.perf_counter()
    self.SymbolicJacobian = self.computeSymbolicJacobian()
    toc = time.perf_counter()
    print(f" Computing simplified symbolic Jacobian done. Took {toc - tic:0.4f} seconds")

  def countJoints(self):
    return len(self.DHParameter)
  
  def forward(self, a, b):
    T_a_b = eye(4)
    for i in range(a,b):
      T_a_b = T_a_b * self.simpl.execute(DH(self.DHParameter[i]))
    return T_a_b

  def forwardKinematics(self):
    return self.CachedPartialTransforms[self.countJoints()-1]

  def fillTransformationCache(self):
    expanded = [self.forward(0, i+1) for i in range(self.countJoints())]
    self.CachedPartialTransforms = [self.simpl.execute(expanded[i]) for i in range(self.countJoints())]

  def computeSymbolicJacobian(self):
    tic = time.perf_counter()
    jac = Jacobian(self.forwardKinematics(), self.JointValueSymbols)
    toc = time.perf_counter()
    print(f"  Computing symbolic Jacobian done. Took {toc - tic:0.4f} seconds")
    tic = time.perf_counter()
    result = self.simpl.execute(jac)
    toc = time.perf_counter()
    print(f"  Simplifying symbolic Jacobian done. Took {toc - tic:0.4f} seconds")

    return result
