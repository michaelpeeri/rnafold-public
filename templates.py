import re


_symbol = re.compile("<<(\S+)>>")

class Evaluator:
    def __init__(self, globals, locals):
        self._globals = globals
        self._locals = locals
    
    def __call__(self, match):
        expression = match.group(1)

        try:
            return str(eval(expression, self._globals, self._locals))
        except Exception as err:
            print("Error in expression <<{}>>".format(expression))
            raise err


def generate(filein, fileout, globals, locals):

    repl = Evaluator(globals, locals)
        
    with open(filein, "r") as fin, open(fileout, "w") as fout:
        
        for line in fin.readlines():
            line = _symbol.sub(repl, line)
            fout.write(line)



