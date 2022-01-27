from sympy import Function
from sympy.abc import t
from sympy.matrices import Matrix, eye

from simplifier import Simplifier

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
    self.fillTransformationCache()
    self.SymbolicJacobian = self.computeSymbolicJacobian()

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
    jac = Jacobian(self.forwardKinematics(), self.JointValueSymbols)
    result = self.simpl.execute(jac)

    return result
