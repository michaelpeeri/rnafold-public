import _distributed
from data_helpers import allSpeciesSource, getSpeciesProperty
from ENCprime import annotateENcPrime

def runDistributed():
    import _distributed
    import dask

    scheduler = _distributed.open()
    delayedCalls = []
    
    for taxId in allSpeciesSource():

        if not getSpeciesProperty( taxId, "ENc-prime" )[0] is None:
            continue

        print(taxId)
        
        call = dask.delayed( annotateENcPrime )(taxId)
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
            (taxId, ENc, ENc_prime, isFreshValue) = scheduler.gather(f)
            results[taxId] = (ENc, ENc_prime)
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
