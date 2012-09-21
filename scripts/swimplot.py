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

# swimplot.py
import getopt
import sys
# Scientific NetCDF API
from Scientific.IO.NetCDF import *
# Nio NetCDF API
#from Nio import *
from numpy import *

def usage():
	print 'Usage: swimplot.py [<options>] <NetCDF file> <time dependent variable>'
	print "\t--colour-min=<minimum value>\tMinimum of colour scale (white if under), implies --minval"
	print "\t--colour-max=<maximum value>\tMaximum of colour scale (black if over), implies --maxval"
	print "\t-d=<time dependent variable>, --diff=<time dependent variable>\tPlot difference <main variable>-<this variable>"
	print "\t-i, --interactive\tPlot interactively (requires X)"
	print "\t--minval=<minimum value>\tMinimum value for colour scale"
	print "\t--maxval=<maximum value>\tMaximum value for colour scale"
	print "\t-n, --dry-run\tDo not plot anything."
	print "\t-t, --transpose\tTranspose value matrix before plotting"
	print "\t-p, --prefix=<prefix>\tPrefix for output files"
	print "\t--quiver=<scale factor>\tOverlay with quiver pot of velocity (scale with factor)"
	sys.exit(2)


try:
        (opts,args) = getopt.getopt(sys.argv[1:],'tpind:',('minval=','maxval=','colour-min=','colour-max=','transpose','prefix=','interactive','diff=','dry-run','quiver='))
except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
	usage()

if len(args) < 2:
	usage()

# get file name from cmd line
filename = args[0]

# load data from netCDF file
# Scientific
netcdf = NetCDFFile(filename,'r')
# Nio
#netcdf = open_file(filename,'r')

# read variable name
var = args[1]

# read coordinate data
# Scientific
lat = netcdf.variables['latitude'].getValue()
# Nio
#lat = netcdf.variables['latitude'].get_value()

# read velocity data
# Scientific
#u = netcdf.variables['u'].getValue()
#v = netcdf.variables['v'].getValue()
# Nio
#u = netcdf.variables['u'].get_value()
#v = netcdf.variables['v'].get_value()

# read dimensions
x = netcdf.dimensions['x']
y = netcdf.dimensions['y']

# read variable to be plotted
# Scientific
pvar = netcdf.variables[var].getValue()
timevar = netcdf.variables['t'].getValue()
# Nio
#pvar = netcdf.variables[var].get_value()
#timevar = netcdf.variables['t'].get_value()

maxt = len(pvar)

minx = 0
maxx = x
miny = 0
maxy = y

minval = 0
maxval = 0

mincol = minval
maxcol = maxval
mncf = False
mxcf = False

dry = False
tp = False
interactive = False
prefix = ''

qscale = None

for opt, arg in opts:
        if opt == '--minval':
                minval = int(arg)
        elif opt == '--maxval':
                maxval = int(arg)
        elif opt in ('-t','--transpose'):
                tp = True
	elif opt in ('--colour-min'):
                mincol = int(arg)
		mncf = True
        elif opt in ('--colour-max'):
                maxcol = int(arg)
		mxcf = True
        elif opt in ('--interactive','-i'):
                interactive = True
	elif opt in ('--prefix','-p'):
                prefix = arg + '-'
        elif opt in ('--diff','-d'):
		# Scientific
                pvar = pvar-netcdf.variables[arg].getValue()
		# Nio
		#pvar = pvar-netcdf.variables[arg].get_value()
		var = arg + '-' + var
        elif opt == '-x':
                (mins,maxs) = arg.split('-')
                minx = int(mins)
                maxx = int(maxs)
        elif opt == '--minx':
                minx = int(arg)
        elif opt == '--maxx':
                maxx = int(arg)
        elif opt == '-y':
                (mins,maxs) = arg.split('-')
                miny = int(mins)
                maxy = int(maxs)
        elif opt == '--miny':
                miny = int(arg)
        elif opt == '--maxy':
                maxy = int(arg)
	elif opt in ('--dry-run', '-n'):
		dry = True
	elif opt == '--quiver':
                qscale = int(arg)

import matplotlib
if not interactive:
	matplotlib.use('Agg')
import pylab as p

if tp:
        print 'Value matrix will be transposed.'

extend='neither'
if mncf:
	print 'Colourbar will have lower boundary of %0d' % (mincol)
	# colour minimum implies value minimum
	minval = mincol
	extend='min'
if mxcf:
	print 'Colourbar will have upper boundary of %0d' % (maxcol)
	# colour maximum implies value minimum
	maxval = maxcol
	if mncf:
		extend='both'
	else:
		extend='max'

# create colour map
colmap = matplotlib.cm.Paired
if mncf:
	colmap.set_under(color='w')
if mxcf:
	colmap.set_over(color='k')

cnorm = matplotlib.colors.Normalize()
if mncf and mxcf:
	cnorm = matplotlib.colors.Normalize(vmin=mincol,vmax=maxcol,clip=False)
elif mncf:
	cnorm = matplotlib.colors.Normalize(vmin=mincol,clip=False)
elif mxcf:
	cnorm = matplotlib.colors.Normalize(vmax=maxcol,clip=False)

# value range not set on cmdline
if (minval == maxval):
        minval = int(floor(1.01*array(pvar).min()))
        maxval = int(ceil(1.01*array(pvar).max()))

print 'Boundaries for values are %0.2f and %0.2f' % (minval,maxval)

if not dry:
	p.imshow(pvar[0,:,:],origin='lower',interpolation='nearest',vmin=minval,vmax=maxval,cmap=colmap,norm=cnorm)
	p.colorbar(orientation='horizontal',cmap=colmap,norm=cnorm,extend=extend);
	filename = prefix + "%s-colorbar.png" % (var)
	p.title('')
	p.show()
	p.savefig(filename)
	p.clf()

for t in range(maxt):


	print 'Plotting time step ' + repr(t)

	if tp:
		pvart = pvar[t,minx:maxx,miny:maxy].transpose()
	else:
		pvart = pvar[t,minx:maxx,miny:maxy]

	if not dry:
		p.imshow(pvart,origin='lower',interpolation='nearest',vmin=minval,vmax=maxval,cmap=colmap,norm=cnorm)

		if qscale != None:
			X,Y = p.meshgrid(linspace(0,int(x-1),int(x/2)),linspace(0,int(x-1),int(x/2)))
			if tp:
				p.quiver(X,Y,v[t,::10,::10],u[t,::10,::10],scale=qscale)
			else:
				p.quiver(X,Y,u[t,::10,::10],v[t,::10,::10],scale=qscale)
		filename = prefix + "%s%04d.png" % (var,t)
		ptime = "%0.3f" % (timevar[t]/3600)
		p.title(var + ' at t=' + ptime + 'h')
		p.savefig(filename)
		
		p.show()
		p.clf()
