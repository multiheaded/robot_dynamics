from sympy import Symbol
from sympy.matrices import Matrix

def SymbolicConstantVector(name, N):
  return Matrix([Symbol('%s%d'%(name, i+1), real=True) for i in range(N)])