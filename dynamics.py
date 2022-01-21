from dataclasses import dataclass

import time

from sympy.abc import t, M
from sympy.core.facts import TautologyDetected
from sympy.sets.ordinals import OmegaPower
from sympy.utilities.decorator import threaded_factory
from kinematics import Kinematics
from sympy import Symbol, latex, diff
from sympy.matrices import Matrix, BlockMatrix
from transforms import trans, Jacobian

from simplifier import Simplifier

@dataclass
class Body:
    CoM:        Matrix  = Matrix.zeros(3,1)
    Mass:       Symbol  = Symbol('m', real = True)
    Inertia:    Matrix  = Matrix.eye(3)

def Steiner(inertia, mass, displacement):
    a = displacement
    Idisp = mass*Matrix([[ a[1]**2+a[2]**2,      -a[0]*a[1],      -a[0]*a[2]],
                         [      -a[0]*a[1], a[0]**2+a[2]**2,      -a[1]*a[2]],
                         [      -a[0]*a[2],      -a[1]*a[2], a[0]**2+a[1]**2]])
    return (inertia+Idisp)

def CylindricBody(mass, length, radius):
    return Body(
            CoM = Matrix([-length/2, 0, 0]),
            Mass = mass,
            Inertia = mass/12*Matrix([[3*radius**2+length**2, 0, 0], [0, 3*radius**2+length**2, 0], [0, 0, 6*length**2]])
            #Steiner(mass/12*Matrix([[3*radius**2+length**2, 0, 0], [0, 3*radius**2+length**2, 0], [0, 0, 6*length**2]]), mass, Matrix([-length/2,0,0]))
        )

class Dynamics:
    Kinematics  = None
    Bodies      = None
    Gravitation = None
    Js          = None
    M           = None
    Omega       = None
    TauE        = None

    Theta       = None
    CORZEN      = None
    Grav        = None
    simpl  = None
    def __init__(self, kinematics, bodies, gravitation, simp = Simplifier("None")):
        self.Kinematics = kinematics
        self.Bodies = bodies
        self.Gravitation = gravitation
        self.simpl = simp

        assert kinematics.countJoints() == len(bodies), 'A body must be supplied for each link in the kinematics object.'

        tic = time.perf_counter()
        Js = self.computeCoMJacobian()
        self.Js = Js
        toc = time.perf_counter()
        print(f"  Computing CoM Jacobian done. Took {toc - tic:0.4f} seconds")

        tic = time.perf_counter()
        M  = self.computeMassMatrix()
        self.M = M
        toc = time.perf_counter()
        print(f"  Computing mass matrix done. Took {toc - tic:0.4f} seconds")

        tic = time.perf_counter()
        O = self.computeAngularBodyVelocity()
        self.Omega = O
        toc = time.perf_counter()
        print(f"  Computing angular velocities done. Took {toc - tic:0.4f} seconds")

        tic = time.perf_counter()
        tauE = self.computeGravitationalForce()
        self.TauE = tauE
        toc = time.perf_counter()
        print(f"  Computing gravitational forces done. Took {toc - tic:0.4f} seconds")

        toc = time.perf_counter()
        print(f" Computing partial matrices done. Took {toc - tic:0.4f} seconds")

        print(" Computing inertia")
        tic = time.perf_counter()
        I = Js.transpose()*M*Js
        toc = time.perf_counter()
        print(f"  Computing Theta done. Took {toc - tic:0.4f} seconds")
        tic = time.perf_counter()
        self.Theta = self.simpl.execute( I )
        toc = time.perf_counter()
        print(f"  Simplifying Theta done. Took {toc - tic:0.4f} seconds")

        print(" Computing coriolis and centrifucal forces")
        tic = time.perf_counter()
        CORZEN = Js.transpose()*( M*diff(Js,t) + O*M*Js) 
        toc = time.perf_counter()
        print(f"  Computing CORZEN done. Took {toc - tic:0.4f} seconds")

        tic = time.perf_counter()
        self.CORZEN = self.simpl.execute( CORZEN )
        toc = time.perf_counter()
        print(f"  Simplifying CORZEN done. Took {toc - tic:0.4f} seconds")

        print(" Computing gravitational forces")
        tic = time.perf_counter()
        Grav = Js.transpose()*tauE
        toc = time.perf_counter()
        print(f"  Computing Grav done. Took {toc - tic:0.4f} seconds")
        tic = time.perf_counter()
        self.Grav = self.simpl.execute( Grav )
        toc = time.perf_counter()
        print(f"  Simplifying Grav done. Took {toc - tic:0.4f} seconds")

    
    def computeCoMJacobian(self):
        Jlist = []
        for i in range(self.Kinematics.countJoints()):
            TSi = self.Kinematics.CachedPartialTransforms[i]*trans(self.Bodies[i].CoM[0], self.Bodies[i].CoM[1], self.Bodies[i].CoM[2])
            Jlist.append(Jacobian(TSi, self.Kinematics.JointValueSymbols))
        return BlockMatrix([[J] for J in Jlist]).as_explicit()

    def computeMassMatrix(self):
        return Matrix.diag([Matrix.diag(b.Mass*Matrix.eye(3),b.Inertia, unpack=True) for b in self.Bodies], unpack=True)

    def computeAngularBodyVelocity(self):
        Olist = []
        for T in self.Kinematics.CachedPartialTransforms:
            # get the rotational submatrix
            R = T.extract([0, 1, 2], [0, 1, 2])
            # transpose to get the inverse
            R_inv = R.transpose()
            # get the skew symmetric angular velocity matrix
            O = diff(R,t)*R_inv
            Olist.append(O)

        return Matrix.diag([Matrix.diag(Matrix.zeros(3),O, unpack=True) for O in Olist], unpack=True)

    def computeGravitationalForce(self):
        return BlockMatrix([ [BlockMatrix([[b.Mass*self.Gravitation],[Matrix.zeros(3,1)]]).as_explicit()] for b in self.Bodies]).as_explicit()

