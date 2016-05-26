# Extract data (from redis) to create RNA energy plot for native CDS fold energy and shuffled CDS fold energy
# Use plot_rnafold_energy_vs_cds_length.r to create plot.
import sys
import redis
import config

# Configuration
# Species (taxIds) for inclusion
taxIdsForProcessing = [3055, 556484]


r = redis.StrictRedis(host=config.host, port=config.port, db=config.db)



for taxIdForProcessing in taxIdsForProcessing:
    #print("#Procesing %d sequences for tax-id %d (%s)..."
    #    % (r.scard("species:taxid:%d:CDS" % taxIdForProcessing),
    #    taxIdForProcessing,
    #    r.get("species:taxid:%d:name" % taxIdForProcessing)))

    skipped = 0
    selected = 0
    for protId in r.sscan_iter("species:taxid:%d:CDS" % taxIdForProcessing):
        # Filters
        # Skip partial CDSs
        if(r.exists("CDS:taxid:%d:protid:%s:partial" % (taxIdForProcessing, protId))):
            skipped += 1
            continue

        # Skip entries missing the native folding energy
        if(not r.exists("CDS:taxid:%d:protid:%s:computed:rna-fold-0:energy" % (taxIdForProcessing, protId))):
            skipped += 1
            continue

        # Skip entries missing the shuffled folding energy
        if(not r.exists("CDS:taxid:%d:protid:%s:computed:rna-fold-for-shuffled-0:energy" % (taxIdForProcessing, protId))):
            skipped += 1
            continue

        selected += 1

        # Read field values
        cdsLength = int(r.strlen("CDS:taxid:%d:protid:%s:seq" % (taxIdForProcessing, protId)))
        foldEnergy = float(r.get("CDS:taxid:%d:protid:%s:computed:rna-fold-0:energy" % (taxIdForProcessing, protId)))
        foldEnergySuffled = float(r.get("CDS:taxid:%d:protid:%s:computed:rna-fold-for-shuffled-0:energy" % (taxIdForProcessing, protId)))

        # Write all data for this protein in CSV format
        print("%d,%s,%d,%f,%f" % (taxIdForProcessing, protId, cdsLength, foldEnergy, foldEnergySuffled))

        #print("#%d selected, %d skipped (%d total)" % (selected, skipped, selected+skipped))
