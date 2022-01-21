from sympy.simplify import simplify, trigsimp

class Simplifier:
    execute = None
    
    def __init__(self, name = "simplify"):
        if name == "simplify":
            self.execute = simplify
        elif name == "trigsimp":
            self.execute = trigsimp
        else:
            self.execute = lambda x : x
