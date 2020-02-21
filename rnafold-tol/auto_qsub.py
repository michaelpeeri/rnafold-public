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
import sys
from templates import generate

count = 0

template = open(sys.argv[1],"r").read()


for index,element in enumerate(sys.stdin.readlines()):
    element = element.rstrip()
    
    generate(sys.argv[1], "%s_%03d.job.sh" % (sys.argv[1], count), {"element":element, "index":index}, {})
    
    count += 1

print("Created %d files."%(count,))
