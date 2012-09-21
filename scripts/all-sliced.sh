#!/bin/bash
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

# all-sliced.sh
export SCRIPTPATH=~/uni/swim/scripts

if [ -z $1 ]
then
	echo "No file specified."
else
	export nc=$1
	`/usr/bin/ncdump $nc | perl -ni -e'if ($_ =~ /^[\sA-Za-z0-9]+\s=.*$/) { s/^\s+([A-Za-z0-9]+)\s=.*$/$1/g; $a = $_; chomp $a; $low = 0; $high = 0 } else { @b = split /,/; foreach $value (@b) { if ($value < $low) { $low = $value * 1.0 } elsif ($value > $high) { $high = $value * 1.0 } } if ( $_ =~ /[\s0-9,];\s*$/) { print "export $a"."_low=$low $a"."_high=$high\n"; }}'`
	/usr/bin/ncdump $nc | perl -ni -e'if ($_ =~ /^[\sA-Za-z0-9]+\s=.*$/) { s/^\s+([A-Za-z0-9]+)\s=.*$/$1/g; $a = $_; chomp $a; $low = 0; $high = 0 } else { @b = split /,/; foreach $value (@b) { if ($value < $low) { $low = $value * 1.0 } elsif ($value > $high) { $high = $value * 1.0 } } if ( $_ =~ /[\s0-9,];\s*$/) { print "export $a"."_low=$low $a"."_high=$high\n"; }}'
	echo "Processing file $nc"
	for varname in phi u v
	do
		export low=`printenv ${varname}_low`
		export high=`printenv ${varname}_high`
		lowlim=`echo "$low*1.2" | bc -l`
		highlim=`echo "$high*1.2" | bc -l`
		echo "${varname} limits are ${low} and ${high}"
		echo "Plotting with limits ${lowlim} and ${highlim}"
		${SCRIPTPATH}/swimplot-slice.py $nc $varname $2 y $lowlim $highlim
		${SCRIPTPATH}/swimplot-slice.py $nc $varname $3 x $lowlim $highlim
#		if [ ${varname} = "phi" ]
#		then
#			~/swim/plot/swimplot-slice.py $nc $varname $2 y $4 $5
#			~/swim/plot/swimplot-slice.py $nc $varname $3 x $4 $5
#		else
#			~/swim/plot/swimplot-slice.py $nc $varname $2 y $6 $7
#			~/swim/plot/swimplot-slice.py $nc $varname $3 x $6 $7
#		fi
		montage -verbose -geometry 1024x768+5+5 `for j in $(seq 0 5 70); do echo ${varname}-sliced$(printf %0.4d $j)x.png; done` ${varname}-slicedx.jpg
		montage -verbose -geometry 1024x768+5+5 `for j in $(seq 0 5 70); do echo ${varname}-sliced$(printf %0.4d $j)y.png; done` ${varname}-slicedy.jpg
	done
fi
