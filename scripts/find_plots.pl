#!/usr/bin/perl
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

# find_plots.pl
use strict;

my @fulltimes = ( 4, 8, 11, 15, 72, 168 );
my @fdiff = ( 10., 10., 10., 10., 10., 10. );
my @ftidx = ( 0, 0, 0, 0, 0, 0 );
my @comptimes = ( 4, 25, 168 );
my @cdiff = ( 10., 10., 10. );
my @ctidx = ( 0, 0, 0 );

my $ptime;
my $i = 0;
while (defined($ptime = <STDIN>)) {

	chomp $ptime;

	my $time;
	my $j = 0;
	foreach $time (@fulltimes) {
		my $crdiff = abs($ptime - $time*3600.)/3600.;
		if ($crdiff < $fdiff[$j]) {
			$fdiff[$j] = $crdiff;
			$ftidx[$j] = $i;
		}
		#print "$ptime $time $crdiff $fdiff[$j] $ftidx[$j] $j\n";
		$j++;
	}

	$j = 0;
	foreach $time (@comptimes) {
		my $crdiff = abs($ptime - $time*3600.)/3600.;
		if ($crdiff < $cdiff[$j]) {
			$cdiff[$j] = $crdiff;	
			$ctidx[$j] = $i;
		}
		$j++;
	}
	$i++;
}
for ($i=0;$i<(scalar @ftidx);$i++) {
	print "$ftidx[$i]\n";
}
for ($i=0;$i<(scalar @ctidx);$i++) {
	print STDERR "$ctidx[$i]\n";
}
