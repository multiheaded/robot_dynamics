import time, sys, getopt
from sympy import Symbol, latex, pi
from sympy.abc import g
from sympy.matrices import Matrix

from kinematics import SymbolicJointValues, Kinematics
from dynamics import Dynamics, Body
from simplifier import Simplifier
from transforms import DHSet

from sympy.interactive.printing import init_printing
init_printing(use_unicode=True, wrap_line=False, pretty_print=True, latex_mode=True)

def SymbolicConstantVector(name, N):
  return Matrix([Symbol('%s%d'%(name, i+1), real=True) for i in range(N)])

def printHelp():
    print('main.py -s <simplifier>')
    print(' simplifier is \"trigsimp\", \"simplify\", or \"None\"')

def main(argv):
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

    try:
        opts, args = getopt.getopt(argv,"hs:",["simplifier="])
    except getopt.GetoptError:
        printHelp()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('main.py -s <simplifier>')
            print(' simplifier is \"trigsimp\", \"simplify\", or \"None\"')
            sys.exit()
        elif opt in ("-s", "--simplifier"):
            simp = Simplifier(arg)

    print(f"Using \"%s\" simplifier"%simp.name)    

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
    SCARA = Kinematics(DHParameters, q, simp)
    toc = time.perf_counter()
    print(f"Kinematics done. Took {toc - tic:0.4f} seconds")

    Bodies = [ Body(CoM=Matrix([x[i], y[i], z[i]]), Mass=m[i]) for i in range(6)]
    Gravitation = Matrix([0,g,0])

    tic = time.perf_counter()
    DynamicSCARA = Dynamics(SCARA, Bodies, Gravitation, simp)
    toc = time.perf_counter()
    print(f"Dynamics done. Took {toc - tic:0.4f} seconds")

    print("Inverting Theta")
    tic = time.perf_counter()
    T_inv = simp.execute(DynamicSCARA.Theta.inv())
    toc = time.perf_counter()
    print(f"Forward dynamics done. Took {toc - tic:0.4f} seconds")

    tic = time.perf_counter()
    L = simp.execute(DynamicSCARA.Theta.cholesky(False))
    toc = time.perf_counter()
    print(f"Decomposition done. Took {toc - tic:0.4f} seconds")

    tic = time.perf_counter()
    testmat = simp.execute(L*L.transpose())
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

if __name__ == "__main__":
   main(sys.argv[1:])