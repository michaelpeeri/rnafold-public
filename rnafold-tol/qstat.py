import re
import subprocess
from cStringIO import StringIO


class qstat(object):

    def __init__(self):
        self.jobs = []
        self.newjob = None
    
    def handleNewJobLine(self, line):
        # Save the existing (previous) job record
        if self.newjob:
            self.jobs.append(self.newjob)

        # Construct the new job record
        assert(line.startswith("Job Id: "))
        newJobId = line[8:]
        isJobArray = False
        dotpos = newJobId.find(".")
        if newJobId[dotpos-2]=="[":  # job arrays
            dotpos -= 2
            isJobArray = True
            
        numericId = int(newJobId[:dotpos])
        self.newjob = {"id":newJobId, "numeric-id":numericId, "qstat.py:is_job_array":isJobArray}

    def handleVarLine(self, line):
        parts = line.split(" = ")
        assert(len(parts)==2)

        self.newjob[parts[0]] = parts[1]


    def handleContinuationLine(self, line):
        pass

    def processQstatLine(self, indentLevel, isContinuation, line):

        if isContinuation:
            return self.handleContinuationLine(line)

        elif indentLevel==0:
            return self.handleNewJobLine(line)

        elif indentLevel==4:
            return self.handleVarLine(line)

        else:
            print(line)
            print(indentLevel)
            assert(False)


    def __call__(self):
        for line in StringIO(subprocess.check_output(("qstat", "-f"), shell=False)):
            line = line.rstrip()

            if not line: continue

            if ord(line[0])==9:
                isContinuation = True
                line = line[1:]
            else:
                isContinuation = False

            indentLevel = len(line) - len(line.lstrip())


            self.processQstatLine(indentLevel, isContinuation, line.strip() )
            
        return self.jobs


if __name__=="__main__":
    import sys
    jobs = qstat()()
    print(len(jobs))
    print(jobs[0])

    numRunning = sum([1 for job in jobs if job['job_state'] == 'R'])
    numQueded  = sum([1 for job in jobs if job['job_state'] == 'Q'])
    print("Running: {} Queued: {}".format(numRunning, numQueded))
    sys.exit(0)
