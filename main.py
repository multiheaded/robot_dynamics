import sys, getopt

from sympy import latex
from simplifier import Simplifier
from URx import computeURxRobot
from SCARA import computeSCARA
from STANF import computeSTANF

from sympy.interactive.printing import init_printing
init_printing(use_unicode=True, wrap_line=False, pretty_print=True, latex_mode=True)

def printHelp():
    print('main.py -r <robot> -s <simplifier> [-k] [-j] [-d]')
    print(' robot is either URx, SCARA, or STANF')
    print(' simplifier is \"trigsimp\", \"simplify\", or \"None\"')
    print(' -k print the forward kinematics')
    print(' -j print the jacobian')
    print(' -d print the dynamic matrics')

def main(argv):
    simp = Simplifier("None")

    computeMechanism = computeSCARA
    printKinematics = False
    printJacobian = False
    printDynamics = False

    try:
        opts, args = getopt.getopt(argv,"hr:s:kjd",["robot=", "simplifier="])
    except getopt.GetoptError:
        printHelp()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            printHelp()
            sys.exit()
        elif opt in ("-s", "--simplifier"):
            simp = Simplifier(arg)
        elif opt in ("-r", "--robot"):
            if arg == "URx":
                computeMechanism = computeURxRobot
            elif arg == "SCARA":
                computeMechanism = computeSCARA
            elif arg == "STANF":
                computeMechanism = computeSTANF
        elif opt in ("-k"):
            printKinematics = True
        elif opt in ("-j"):
            printJacobian = True
        elif opt in ("-d"):
            printDynamics = True

    kin, dyn = computeMechanism(simp)
    
    if printKinematics:
        print(latex(kin.forwardKinematics()))

    if printJacobian:
        print(latex(kin.SymbolicJacobian))

    if printDynamics:
        print(latex(dyn.Theta))
        print("")
        print(latex(dyn.CORZEN))
        print("")
        print(latex(dyn.Grav))


if __name__ == "__main__":
   main(sys.argv[1:])