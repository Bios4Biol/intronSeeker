#!/bin/bash

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


echo -e "\n### Begin installing intronSeeker...\n"

# Add the launcher and the scripts to the PATH by .bashrc
HERE=$(dirname $(realpath $0))

# Export the conda functions in the subshell
source activate
source $CONDA_PREFIX/etc/profile.d/conda.sh

# Conda environment installing
echo -e "### Conda environment installing... \n"
if [ $(conda env list | grep 'ISeeker_environment' | wc -l) = 0 ]
then
    mamba env create --file $HERE/config/environment.yml
    
    if [ $? = 0 ] 
    then 
        echo -e "### Conda environment installed. \n"
    
        # Conda environment refinement
        echo -e "### Conda environment configuration...\n"
    
        conda activate ISeeker_environment
    else
        code=$?
        echo "FAIL : Conda environment installation FAILED."
        exit $code
    fi
    
else
    echo -e "Conda environment already installed. \n"
    echo -e "### Conda environment configuration...\n"
    
    conda activate ISeeker_environment
fi




if [ $? = 0 ]
then
    ln -s $CONDA_PREFIX/lib/libcrypto.so.1.1 $CONDA_PREFIX/lib/libcrypto.so.1.0.0
    
    mkdir -p $CONDA_PREFIX/etc/conda/activate.d/
    cp config/set_env_config.sh $CONDA_PREFIX/etc/conda/activate.d/
    echo -e "export PATH=$HERE:$HERE/scripts/:\$PATH" >> $CONDA_PREFIX/etc/conda/activate.d/set_env_config.sh
    
    mkdir -p $CONDA_PREFIX/etc/conda/deactivate.d/
    cp config/unset_env_config.sh $CONDA_PREFIX/etc/conda/deactivate.d/
    
    cp -a $HERE/bin/assemblathon_stats.pl $CONDA_PREFIX/bin/
    cp -a $HERE/bin/FAlite.pm $CONDA_PREFIX/bin/
    
    sed -i 's#/tmp/build/80754af9/perl_1527832170752/_build_env/#'"$CONDA_PREFIX"'/#g' $CONDA_PREFIX/lib/5.26.2/x86_64-linux-thread-multi/Config.pm
    sed -i 's#/tmp/build/80754af9/perl_1527832170752/_build_env/#'"$CONDA_PREFIX"'/#g' $CONDA_PREFIX/lib/5.26.2/x86_64-linux-thread-multi/Config_heavy.pl
    sed -i 's#/tmp/build/80754af9/perl_1527832170752/_build_env/#'"$CONDA_PREFIX"'/#g' $CONDA_PREFIX/lib/5.26.2/x86_64-linux-thread-multi/CORE/config.h
    
    conda deactivate
else
    code=$?
    echo "FAIL : Conda environment activation FAILED."
    exit $code
fi

if [ $? = 0 ]
then
    echo -e "### Conda environment configured.\n"
    
    echo -e "### Grinder Installation...\n"
    
    conda activate ISeeker_environment
else
    code=$?
    echo "FAIL : Conda environment configuration FAILED."
    exit $code
fi

if [ $? = 0 ]
then
    mkdir $CONDA_PREFIX/tmp
    cd $CONDA_PREFIX/tmp/
    
    echo -e "\n###\tGetOpt::Euclid Installation...\n"
    
    wget https://cpan.metacpan.org/authors/id/F/FA/FANGLY/Getopt-Euclid-0.4.5.tar.gz
else
    code=$?
    echo "FAIL : Conda environment activation FAILED."
    exit $code
fi

if [ $? = 0 ]
then
    tar -xzf Getopt-Euclid-0.4.5.tar.gz
    cd Getopt-Euclid-0.4.5
    perl Makefile.PL
    make
    make install
else
    code=$?
    echo "FAIL : Perl lib (GetOpt::Euclid) downloading FAILED."
    exit $code
fi

if [ $? = 0 ]
then
    cd ..
    echo -e "\n###\tGetOpt::Euclid Installed\n"
    echo -e "\n###\tTest::Number::Delta Installation...\n"
    wget https://cpan.metacpan.org/authors/id/D/DA/DAGOLDEN/Test-Number-Delta-1.06.tar.gz
else 
    code=$?
    echo "FAIL : Perl lib (GetOpt::Euclid) installation FAILED."
    exit $code
fi

if [ $? = 0 ]
then
    tar -xzf Test-Number-Delta-1.06.tar.gz
    cd Test-Number-Delta-1.06
    perl Makefile.PL
    make 
    make install
else 
    code=$?
    echo "FAIL : Perl lib (Test::Number::Delta) download FAILED."
    exit $code
fi
 

if [ $? = 0 ]
then
    cd ..
    
    echo -e "\n###\tTest::Number::Delta Installed\n"
    echo -e "\n###\tMath::Random::MT Installation...\n"
    
    wget https://cpan.metacpan.org/authors/id/F/FA/FANGLY/Math-Random-MT-1.17.tar.gz
else
    code=$?
    echo "FAIL : Perl lib (Test::Number::Delta) installation FAILED."
    exit $code
fi

if [ $? = 0 ]
then
    tar -xzf Math-Random-MT-1.17.tar.gz
    cd Math-Random-MT-1.17
    perl Makefile.PL
    make
    make install
else
    code=$?
    echo "FAIL : Perl lib (Math::Random::MT) download FAILED."
    exit $code
fi

if [ $? = 0 ]
then
    cd ..
    
    echo -e "\n###\tMath::Random::MT Installed\n"
    
    wget https://sourceforge.net/projects/biogrinder/files/latest/download
else
    code=$?
    echo "FAIL : Perl lib (Math::Random::MT) installation FAILED."
    exit $code
fi

if [ $? = 0 ]
then
    tar -xzf download
    cd Grinder-0.5.4
    perl Makefile.PL
    make 
    make install
else
    code=$?
    echo "FAIL : Grinder download FAILED."
    exit $code
fi

cd ..
cd $HERE

chmod -R 777 $CONDA_PREFIX/tmp/
rm -r $CONDA_PREFIX/tmp/

echo -e "\n### Grinder Installed\n"

echo -e "### Installation done.\n"
