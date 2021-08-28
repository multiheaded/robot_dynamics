from sympy import Function
from sympy.abc import t
from sympy.matrices import Matrix
from sympy.simplify.simplify import simplify

from transforms import DH, Jacobian

def SymbolicJointValues(N):
  return Matrix([Function('q%d'%(i+1), real=True)(t) for i in range(N)])

class Kinematics:
  DHParameter = None
  CachedPartialTransforms = None
  SymbolicJacobian = None
  JointValueSymbols = None
  
  def __init__(self, DHParms, jointSyms):
    self.JointValueSymbols = jointSyms
    self.DHParameter = DHParms
    self.fillTransformationCache()
    self.SymbolicJacobian = self.computeSymbolicJacobian()

  def countJoints(self):
    return len(self.DHParameter)
  
  def forward(self, a, b):
    T_a_b = Matrix.eye(4)
    for i in range(b-a):
      T_a_b = T_a_b * DH(self.DHParameter[i])
    return simplify(T_a_b)

  def forwardKinematics(self):
    return self.CachedPartialTransforms[self.countJoints()-1]

  def fillTransformationCache(self):
    self.CachedPartialTransforms = [self.forward(0, i+1) for i in range(self.countJoints())]

  def computeSymbolicJacobian(self):
    return simplify(Jacobian(self.forwardKinematics(), self.JointValueSymbols))
