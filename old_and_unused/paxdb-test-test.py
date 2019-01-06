import csv
import sys

thres = float(sys.argv[1])

aboveCount = 0
belowCount = 0
total = 0

with open("4932-WHOLE_ORGANISM-integrated.txt.v2.converted.csv", "rb") as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    for row in reader:
        val3 = float(row[2])
        if( val3 < thres ):
            belowCount += 1
        else:
            aboveCount += 1
        total += 1

print("Below: %.2g%% Above: %.2g%%" % (float(belowCount)/total, float(aboveCount)/total))

        
