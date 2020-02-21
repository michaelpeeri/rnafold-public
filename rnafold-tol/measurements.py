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
    
    
    
