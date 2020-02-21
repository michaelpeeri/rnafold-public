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
from ete3 import NCBITaxa


ncbiTaxa = NCBITaxa()

def getDomainForSpecies(taxId):
    lineage = frozenset(ncbiTaxa.get_lineage(taxId))
    if 2 in lineage:
        return "Bacteria"
    elif 2157 in lineage:
        return "Archaea"
    elif 2759 in lineage:
        return "Eukaryotes"
    else:
        raise Exception("Couldn't determine domain of species {}".format(taxId))

