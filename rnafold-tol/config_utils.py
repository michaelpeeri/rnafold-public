import re
import config


reMysqlConn = re.compile("mysql://(\w+):([^@]+)@([^:]+):(\d+)/(\w+)")

def getMysqlConnectionParams():
    return reMysqlConn.match(config.mysql_host_connection).groups()


    
