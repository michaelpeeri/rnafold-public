import timeit
import subprocess
from config_utils import getMysqlConnectionParams 


def mysql_select():

    # time ~/anaconda2N/bin/mysql -u rnafold --password=123456789 -h server -D rnafold -P 14404 --execute="select alphabet, source, count(*) from sequences2 group by alphabet, source;"
    sqlStatement = "select alphabet, source, count(*) from sequences2 where id between 5000000 and 25000000 group by alphabet, source;"

    (user, passw, host, port, dbname) = getMysqlConnectionParams()
    subprocess.call( ("/tamir1/mich1/anaconda2N/bin/mysql", "-u", user, "--password={}".format(passw), "-h", host, "-D", dbname, "-P", port, "-e", sqlStatement ), shell=False )




def measure(func):
    t = timeit.Timer(stmt=func)
    return t.timeit(number=1)

print(measure(mysql_select))
    
    
    
