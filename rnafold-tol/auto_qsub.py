import sys
from templates import generate

count = 0

template = open(sys.argv[1],"r").read()


for index,element in enumerate(sys.stdin.readlines()):
    element = element.rstrip()
    
    generate(sys.argv[1], "%s_%03d.job.sh" % (sys.argv[1], count), {"element":element, "index":index}, {})
    
    count += 1

print("Created %d files."%(count,))
