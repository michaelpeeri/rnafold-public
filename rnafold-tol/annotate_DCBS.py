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
import _distributed
from data_helpers import allSpeciesSource
from DCBS import annotateDCBS

def runDistributed():
    import _distributed
    import dask

    scheduler = _distributed.open()
    delayedCalls = []
    
    for taxId in allSpeciesSource():
        call = dask.delayed( annotateDCBS )(taxId)
        delayedCalls.append( call )

    print("Starting %d calls..." % len(delayedCalls))
    
    futures = scheduler.compute(delayedCalls) # submit all delayed calculations; obtain futures immediately

    try:
        _distributed.progress(futures) # wait for all calculations to complete
    except Exception as e:
        print(E)
    print("\n")

    print("Waiting for all tasks to complete...")
    _distributed.wait(futures)

    results = {}
    errorsCount = 0
    newValuesCount = 0
    oldValuesCount = 0
    for f in futures:
        try:
            (taxId, DCBS, isFreshValue) = scheduler.gather(f)
            results[taxId] = DCBS
            if isFreshValue:
                newValuesCount += 1
            else:
                oldValuesCount += 1
            
        except Exception as e:
            print(e)
            errorsCount += 1
            
    print("Finished %d species with %d errors" % (len(results), errorsCount))
    print("{} new values; {} old values".format(newValuesCount, oldValuesCount))
    return results
            
        
print(runDistributed())        
