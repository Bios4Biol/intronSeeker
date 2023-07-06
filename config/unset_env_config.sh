#/bin/bash

#<IntronSeeker searches introns by splice-realigning reads on contigs.>
#Copyright (C) <2019> INRAE
#<Sarah Maman, Philippe Bardou, Emilien Lasguignes, Faustine Oudin, Floréal Cabanettes, Christophe Klopp>
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <https://www.gnu.org/licenses/>.


if [ -e $CONDA_PREFIX/.perl5lib ]
then
	export PERL5LIB=`cat $CONDA_PREFIX/.perl5lib`
	rm $CONDA_PREFIX/.perl5lib
fi
if [ -e $CONDA_PREFIX/.per5lib ]
then
	export PER5LIB=`cat $CONDA_PREFIX/.per5lib`
	rm $CONDA_PREFIX/.per5lib
fi
if [ -e $CONDA_PREFIX/.perl_local_lib_root ]
then
	export PERL_LOCAL_LIB_ROOT=`cat $CONDA_PREFIX/.perl_local_lib_root`
	rm $CONDA_PREFIX/.perl_local_lib_root
fi
if [ -e $CONDA_PREFIX/.perl_mb_opt ]
then
	export PERL_MB_OPT=`cat $CONDA_PREFIX/.perl_mb_opt`
	rm $CONDA_PREFIX/.perl_mb_opt
fi
if [ -e $CONDA_PREFIX/.perl_mm_opt ]
then
	export PERL_MM_OPT=`cat $CONDA_PREFIX/.perl_mm_opt`
	rm $CONDA_PREFIX/.perl_mm_opt
fi

if [ -e $CONDA_PREFIX/.path ]
then
	export PATH=`cat $CONDA_PREFIX/.path`
	rm $CONDA_PREFIX/.path
fi
