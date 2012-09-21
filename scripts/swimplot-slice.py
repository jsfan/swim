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

# swimplot-slice.py
import getopt
import MA
import sys
import pylab as p
from Scientific.IO.NetCDF import *
from numpy import *

def usage():
	print 'Usage: swimplot-slice.py [<options>] <NetCDF file> <time dependent variable> <dimension (x/y)>'
	print "\t--minval\tMinimum value for colour scale"
	print "\t--maxval\tMaximum value for colour scale"
	print "\t-p, --prefix\tprefix for output files"
	print "\t-c, --position=\tposition to slice through 2D field"
	sys.exit(2)

try:
        (opts,args) = getopt.getopt(sys.argv[1:],'p:c:',('minval=','maxval=','prefix=','position='))
except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
	usage()

if len(args) < 2:
	usage()

# get file name from cmd line
filename = args[0]

# load data from netCDF file
netcdf = NetCDFFile(filename,'r')

# read variable name
var = args[1]

# slice in x or y
dim = args[2]

# read coordinate data
lat = netcdf.variables['latitude'].getValue()

# read variable to be plotted
pvar = netcdf.variables[var].getValue()
timevar = netcdf.variables['t'].getValue()
maxt = len(pvar)

minval = 0
maxval = 0
pos = '0'
ipos = 0
prefix = ''

for opt, arg in opts:
        if opt == '--minval':
                minval = int(arg)
        elif opt == '--maxval':
                maxval = int(arg)
        elif opt in ('--prefix','-p'):
                prefix = arg + '-'
        elif opt in ('--position','-c'):
		pos = arg
                ipos = int(arg)

if dim != 'x' and dim != 'X':
	dim = 'y'
else:
	dim = 'x'

if ipos == 0:
	if dim == 'x':
		ipos = int(round(shape(pvar)[1]/2))
	else:
		ipos = int(round(shape(pvar)[2]/2))
	pos = repr(ipos)

if minval == maxval:
	if dim == 'x':
		minval = int(floor(MA.minimum(pvar[:,:,ipos])))
		maxval = int(ceil(MA.maximum(pvar[:,:,ipos])))
	else:
		minval = int(floor(MA.minimum(pvar[:,ipos,:])))
		maxval = int(ceil(MA.maximum(pvar[:,ipos,:])))
	minval = minval - abs(minval)*0.05
	maxval = maxval + abs(maxval)*0.05


sys.stdout.write("Plotting "+var+" (sliced at "+dim+" position "+pos+")\nPlotting time step ")

for t in range(maxt):
    
	sys.stdout.write(repr(t)+'...')
	sys.stdout.flush()
	
	if dim == 'x':
		p.plot(pvar[t,:,ipos])	
	else:
		p.plot(pvar[t,ipos,:])
   
	ax = p.gca()
	ax.set_ylim(minval,maxval)
		
	filename = prefix + "%s-sliced%04d%s.png" % (var,t,dim)
	ptime = "%0.3f" % (timevar[t]/3600)
	p.title(var + ' at t=' + ptime + 'h')
	p.savefig(filename)
	
#	p.show()
	p.clf()

sys.stdout.write("done\n")
