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

# epsilon-plot.py
import getopt
import sys
import matplotlib
matplotlib.use('Agg')
import pylab as p
from Scientific.IO.NetCDF import *
from numpy import *

def usage():
	print 'Usage: epsilon-plot.py [<options>] <NetCDF file>'
	print "\t--minval\tMinimum value for colour scale"
	print "\t--maxval\tMaximum value for colour scale"
	print "\t-d,--dynamic\tOff-centring is dynamic => plot maks for every time step"
	print "\t-t, --transpose\ttranspose value matrix before plotting"
	print "\t-p, --prefix\tprefix for output files"
        sys.exit(2)     

try:
        (opts,args) = getopt.getopt(sys.argv[1:],'tp:d',('minval=','maxval=','transpose','prefix=','dynamic'))
except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
	usage()

if len(args) < 1:
	usage()

# get file name from cmd line
filename = args[0]

# load data from netCDF file
netcdf = NetCDFFile(filename,'r')

# read coordinate data
lat = netcdf.variables['latitude'].getValue()

# read dimensions
x = netcdf.dimensions['x']
y = netcdf.dimensions['y']

# read variable to be plotted
timevar = netcdf.variables['t'].getValue()

minval = 0
maxval = 0

tp = 0
prefix = ''

dyn = 0

for opt, arg in opts:
        if opt == '--minval':
                minval = int(arg)
        elif opt == '--maxval':
                maxval = int(arg)
        elif opt in ('-t','--transpose'):
                tp = 1
        elif opt in ('-d','dynamic'):
                dyn = 1
        elif opt in ('--prefix','-p'):
                prefix = arg + '-'

if tp:
        print 'Value matrix will be transposed.'

print 'Boundaries for values are %0.2f and %0.2f' % (minval,maxval)

step = int(ceil((maxval-minval)/1024.0))

for i in range(1,4):
	# read variable to be plotted
	pvar = netcdf.variables['epsilon'+repr(i)].getValue()
	timevar = netcdf.variables['t'].getValue()
	maxt = len(pvar)

	# colourbar range not set on cmdline
	if (minval == maxval):
	        minval = int(floor(1.01*array(pvar).min()))
	        maxval = int(ceil(1.01*array(pvar).max()))

	print 'Plotting epsilon'+repr(i)+'...'

	if dyn == 1:
		for t in range(maxt):
			
			print 'Plotting time step '+repr(t)
			
			
			p.imshow(pvar[t,:,:].transpose(),origin='lower',interpolation='nearest',vmin=minval,vmax=maxval)
			p.colorbar(orientation="horizontal")

			filename = "epsilon%04d.png" % (t)
			ptime = "%0.3f" % (timevar[t]/3600)
			p.title('epsilon'+repr(i)+' at t=' + ptime + 'h')
			p.savefig(filename)

			#p.show()
			p.clf()
	else:
		
		print 'Plotting first epsilon'+repr(i)+' as global.'
		
		p.imshow(pvar[0,:,:].transpose(),origin='lower',interpolation='nearest')
		p.colorbar(orientation="horizontal")
		
		filename = "epsilon%01d.png" % (i)
		p.title('epsilon'+repr(i))
		p.savefig(filename)
		
		#p.show()
		p.clf()
