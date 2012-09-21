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

# plotcombine-horizontal.sh
if [ "$1" == "" ]
then
	echo "Usage: plotcombine.sh <series name> [comparison plot]"
	exit
fi
# read in files to plot
count=0
list=""
while read plotnumber
  do
  let count=${count}+1
  plotno=`printf "%04d" ${plotnumber}`
  list=`echo "${list} "``echo ${plotno}`
done
name=`echo horizontal-$1`
comp=$2

for i in phi u v
  do
  #for j in ${i}*[0-9].png; do convert -crop 560x600+119+0 ${j} `basename ${j} .png`-cropped.png; done
  #convert -crop 640x80+90+450 ${i}-colorbar.png ${i}-colorbar-cropped.png
  if [ "${comp}" != "" ]
      then
      cpath=`echo ${comp} | sed -e"s:%%:${i}:"`
      montage -geometry 280x300 -tile ${count}x1 `for j in ${list}; do echo ${i}${j}-cropped.png; done` ${i}-${name}-short.png
      montage -geometry `echo "${count}*280" | bc`x300 -tile 1x2 ${cpath} ${i}-${name}-short.png ${i}-${name}-comparison.png
      montage -geometry `echo "${count}*280" | bc`x600\>+0+0 -gravity south -tile 1x2 ${i}-colorbar-cropped.png ${i}-${name}-comparison.png ${i}-${name}.tmp.png
      convert -crop `echo "${count}*280" | bc`x680+0+520 ${i}-${name}.tmp.png ${i}-${name}-comparison.png
  else
      if [ `echo "${count}%2" | bc` == "1" ]
	  then
	  let count=count+1
      fi
      ytiles=`echo "${count}/2" | bc`
      montage -geometry 280x300 -tile ${ytiles}x2 `for j in ${list}; do echo ${i}${j}-cropped.png; done` ${i}-${name}.png
      montage -geometry `echo "${ytiles}*280" | bc`x600\>+0+0 -gravity south -tile 1x2 ${i}-colorbar-cropped.png ${i}-${name}.png ${i}-${name}.tmp.png
      convert -crop `echo "${ytiles}*280" | bc`x680+0+520 ${i}-${name}.tmp.png ${i}-${name}.png
  fi
  rm ${i}-${name}.tmp.png
done
