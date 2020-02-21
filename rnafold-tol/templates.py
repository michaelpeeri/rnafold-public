# RNAFold - Analyze mRNA folding bias (Local Folding Energy) based on randomizations.
# Copyright (C) 2016-2020 Michael Peeri
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
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



