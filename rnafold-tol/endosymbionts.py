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

knownEndosymbionts = frozenset((107806,203907,322098,331104,228908,1236703,1116230,218497,115713,353152,347515,5693,36329,99287,400667,83332,203267,1227812,272631,508771,224914,262768,272633,169963,227377,266834,264462,862908,283166,1321371,1208920,138677,331113,1165094))

def isEndosymbiont(taxId):
    return int(int(taxId) in knownEndosymbionts)
