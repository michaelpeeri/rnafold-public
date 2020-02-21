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
import logging
# Site configuration file (for server details, etc.)

# redis
host = "---redis-server---"
port = 6379
db = 0
password = "---password---"

# sqlite
sqlite_base_path = "/a/home/cc/students/enginer/mich1/rnafold/data/"
def make_sqlite_host_connection(filename):
    return 'sqlite:///%s/%s.db' % (sqlite_base_path, filename)

# mysql
mysql_host_connection = 'mysql://---username---:---password---@---host---:---port---/rnafold'
run_without_mysql_server = False

# data files
base_data_dir    = "/a/home/cc/students/enginer/mich1/rnafold/data/"
ensembl_data_dir = "/a/home/cc/students/enginer/mich1/rnafold/data/Ensembl/"


# Pushover
# Application token
computation_monitor_app = '---pushover-token---'
# Group key
computation_monitor_group = '---pushover-group---'

# Matlab
MatlabPath = "/usr/local.cc/bin/matlab"
# codonW
codonwBasePath = '/a/home/cc/students/enginer/mich1/src/codonW'
# ENCprime
ENCprimeBasePath = '/a/home/cc/students/enginer/mich1/src/ENCprime-master/bin/'

def process_id():
    import subprocess
    import os
    return "%s_%d" % (subprocess.check_output('hostname')[:-1], os.getpid())
 
# logging
def initLogging():
	logging.basicConfig(filename='process_%s.log' % process_id(), level=logging.INFO)


