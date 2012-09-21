#!/usr/bin/python
# SWiM - a semi-Lagrangian, semi-implicit shallow water model in
# Cartesian coordiates
# Copyright (C) 2008-2012 Christian Lerrahn
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

# plottimes.py
from Scientific.IO.NetCDF import *
import sys

def usage():
	print 'Usage: swimplot.py <NetCDF file>'
	sys.exit(2)


if len(sys.argv) < 1:
	usage()

# get file name from cmd line
filename = sys.argv[1]

# load data from netCDF file
netcdf = NetCDFFile(filename,'r')

# read variable to be plotted
timevar = netcdf.variables['t'].getValue()

for t in range(len(timevar)):
	print "%0.0f" % (timevar[t])
