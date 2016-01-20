#!/bin/bash
if [ "$1" = "" ] ||  [ "$2" = "" ]
then
echo Usage: $0 filename element-name
file=out_2d_00001
element=1
	
else
file=$1
element=$2
fi

cat << EOF  > gnu.plt
set size square
set pm3d map
#set cntrparam levels 100
#set samples 100
#set contour base
splot '$file' using $element
pause -1
EOF
gnuplot gnu.plt
