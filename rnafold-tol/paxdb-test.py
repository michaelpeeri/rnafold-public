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
import csv
import requests
from time import sleep
from data_helpers import SpeciesCDSSource


f = sys.argv[1]
taxid = int(sys.argv[2])

#http://rest.ensemblgenomes.org/xrefs/id/CAB12275?content-type=application/json;all_levels=1

server = "http://rest.ensemblgenomes.org"
ext = "/xrefs/id/%s?content-type=application/json;all_levels=1"


def getRecord(protid):
        r = requests.get(
            server + ext % protid)
            #headers={ "Content-Type" : "application/json"})
 
        if not r.ok:
            r.raise_for_status()
            sys.exit()

        # Rate-limit requests
        sleep(0.5)
 
        decoded = r.json()
        return decoded

for protid in SpeciesCDSSource(taxid):
    record = getRecord(protid)
    print(record)
        

#with open(f, 'r') as csvfile:
#    for row in csv.reader(csvfile, delimiter='\t'):
#        #['3706992', '224308.Bsubs1_010100004063', '4.46']

#        paxId = row[1].split(".")[1]
#        pa = float(row[2])

        
        

