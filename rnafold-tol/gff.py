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
# Access GFF data sources
import gffutils  # See: http://daler.github.io/gffutils/contents.html

# Get a GFF FeatureDB object for a given .gff file; Implement policies for this
def createGffDb( filename, variant ):

    db = gffutils.create_db(filename, dbfn="%s.db" % filename, force=True, keep_order=True,
                            sort_attribute_values=True,
                            merge_strategy='merge',
                            id_spec=['ID', 'Name'])
    return db



