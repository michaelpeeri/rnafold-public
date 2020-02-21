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
from timeit import default_timer # import the best platform- and version-specific timer

"""
Performance timer for profiling;
call start() and stop() around measured code;
stats() will return cummulative time spent in section
"""
class Timer(object):
    def __init__(self):
        self.reset()

    def reset(self):
        self._totalDuration = 0.0
        self._numIntervals = 0
        self._start = None

    def start(self):
        if( not self._start is None ):
            raise Exception("Timer started while already running (duration=%f secs)" % (default_timer()-self._start))
        self._start = default_timer()

    def abort(self):
        self._start = None

    def stop(self):
        duration = default_timer() - self._start
        assert(duration >= 0.0)
        self._start = None
        self._totalDuration += duration
        self._numIntervals += 1

    def stats(self):
        return (self._totalDuration, self._numIntervals)

