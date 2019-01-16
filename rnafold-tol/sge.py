import subprocess
import os.path
import sys

def qsub_shell_script(script, **kw):
    
    if not "name" in kw:
        kw["name"] = os.path.basename(script)

    out = None
    try:
        out = subprocess.check_output((
            "qsub",
            "-q", "tamirs",
            "-N", kw["name"][:9],
            "-l", "mem=1800mb,pmem=1800mb,vmem=3000mb,pvmem=3000mb,cput=4:00:00",
            script), shell=False, stderr=subprocess.STDOUT # capture errors
        )
    except subprocess.CalledProcessError as e:
        print(out)
        print(e)
        raise

if __name__=="__main__":
    if len(sys.argv)>=2:
        command = sys.argv[1]

        if( command == "run" ):
            if len(sys.argv)<3:
                raise Exception("run: missing arguments")
            sys.exit(qsub_shell_script(sys.argv[2]))
        else:
            raise Exception("Unknown command: %s" % command)
