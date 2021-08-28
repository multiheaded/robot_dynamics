from sympy import cos, sin, diff, symbols, shape
from sympy.abc import t
from sympy.matrices import Matrix, BlockMatrix
from sympy.solvers.solveset import linear_eq_to_matrix

from dataclasses import dataclass

@dataclass
class DHSet:
    theta:  float = 0.0
    d:      float = 0.0
    a:      float = 0.0
    alpha:  float = 0.0

def rotx(var):
    return Matrix([[1, 0, 0, 0], [0, cos(var), -sin(var), 0], [0, sin(var), cos(var), 0], [0, 0, 0, 1]])

def roty(var):
    return Matrix([[cos(var), 0, sin(var), 0], [0, 1, 0, 0], [-sin(var), 0, cos(var), 0], [0, 0, 0, 1]])

def rotz(var):
    return Matrix([[cos(var), -sin(var), 0, 0], [sin(var), cos(var), 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])

def trans(x, y, z):
    return Matrix([[1, 0, 0, x], [0, 1, 0, y], [0, 0, 1, z], [0, 0, 0, 1]])

def DH(dhSet):
    return rotz(dhSet.theta)*trans(0,0,dhSet.d)*trans(dhSet.a,0,0)*rotx(dhSet.alpha)

def Jacobian(mat, syms):
    # get the position vector via extract(rowlist, collist), MATLAB: mat(1:3, 4)
    p = mat.extract([0, 1, 2], [3])

    # use the position vector to derive the jacobian using built-in functions, MATLAB: jacobian(p, syms);
    Jv = p.jacobian(syms)

    # get the rotational submatrix
    R = mat.extract([0, 1, 2], [0, 1, 2])
    # transpose to get the inverse
    R_inv = R.transpose()
    # get the skew symmetric angular velocity matrix
    O = diff(R,t)*R_inv
    # extract the x, y, and z components
    w = Matrix([O[2,1], O[0,2], O[1,0]])

    rows = shape(syms)[0]
    dsyms = diff(syms, t)
    # create a list of temporary symbols
    tmpq = symbols('tmpq0:%d'%rows, real=True)
    # make a list of tuples containing temporary replacements, e.g. [(diff(q1,t),tmpq1), ...]
    replacements = list(zip(dsyms,tmpq))
    tmpw = w.subs(replacements)

    Jw, b = linear_eq_to_matrix(tmpw, tmpq)
    
    return BlockMatrix([[Jv],[Jw]]).as_explicit()