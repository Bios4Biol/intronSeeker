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


if [ -n $PERL5LIB ]
then
	echo $PERL5LIB > $CONDA_PREFIX/.perl5lib
	export PERL5LIB=""
fi
if [ -n $PER5LIB ]
then
	echo $PER5LIB > $CONDA_PREFIX/.per5lib
	export PER5LIB=""
fi
if [ -n $PERL_LOCAL_LIB_ROOT ]
then
	echo $PERL_LOCAL_LIB_ROOT > $CONDA_PREFIX/.perl_local_lib_root
	export PERL_LOCAL_LIB_ROOT=""
fi
if [[ -n $PERL_MB_OPT ]]
then
	echo $PERL_MB_OPT > $CONDA_PREFIX/.perl_mb_opt
	export PERL_MB_OPT='--install_base ""'
fi
if [ -n $PERL_MM_OPT ]
then
	echo $PERL_MM_OPT > $CONDA_PREFIX/.perl_mm_opt
	export PERL_MM_OPT=INSTALL_BASE=
fi


echo $PATH > $CONDA_PREFIX/.path
