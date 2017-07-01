import logging
# Site configuration file (for server details, etc.)

# redis
host = "power5"
port = 6379
db = 0
password = "rnafold"

# sqlite
sqlite_base_path = "/home/michael/rnafold/data"
def make_sqlite_host_connection(filename):
    return 'sqlite:///%s/%s.db' % (sqlite_base_path, filename)

# mysql
mysql_host_connection = 'mysql://rnafold:--password--@ec2-54-234-142-88.compute-1.amazonaws.com/rnafold'

# Pushover
# Application token
computation_monitor_app = '--pushover-key--'
# Group key
computation_monitor_group = '--pushover-key--'

# codonW
codonwBasePath = "/home/michael/src/codonW"

def process_id():
    import subprocess
    import os
    return "%s_%d" % (subprocess.check_output('hostname')[:-1], os.getpid())
    

# Logging
logging.basicConfig(filename='process_%s.log' % process_id(), level=logging.INFO)

