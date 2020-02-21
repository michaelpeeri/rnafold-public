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
import sys
import json
import requests
import config


"""
Send simple notification via Push-over.
Note: Limited to 7500 messages/month.
TODO: Add rate limiting, better error handling.
"""
def notify(message):
    payload = { 'token': config.computation_monitor_app, 'user': config.computation_monitor_group, 'message': message }
    r = requests.post('https://api.pushover.net/1/messages.json', data=payload)
    if( r.status_code != 200 ):
        raise Exception("Pushover returned code %d, content %s" % r.status_code, r.json())
    return (r.status_code, r.json())


def defaultNotification():
    from socket import gethostname
    notify("Notification (%s)" % gethostname())


if __name__=="__main__":
    defaultNotification()
    sys.exit(0)
