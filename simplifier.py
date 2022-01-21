from unicodedata import name
from sympy.simplify import simplify, trigsimp

class Simplifier:
    execute = None
    name = None
    
    def __init__(self, name = "simplify"):
        self.name = name
        if name == "simplify":
            self.execute = simplify
        elif name == "trigsimp":
            self.execute = trigsimp
        else:
            name = "None"
            self.execute = lambda x : x
