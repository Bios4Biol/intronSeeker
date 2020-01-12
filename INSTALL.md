How to Install ?
================

Hereafter, you will find all the details about the intronSeeker installation.

Requirements
------------

IntronSeeker requires Python version 3.6 or above.

It also needs Python packages and external softwares to work 
correctly.
So, to make install easier, a conda enviromnent (grinder excluded)
has been created  : it is the [environment.yml](https://forgemia.inra.fr/emilien.lasguignes/intronSeeker/blob/master/config/environment.yml) 
file available in intronSeeker/config/ directory. 

For an easy install, a setup.sh script was develloped to install conda environment, to configure it and to install Grinder 
but a conda installation is needed . 
To install Miniconda, follow [this tutorial](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) or,
to install Ananconda, follow [this one](https://docs.anaconda.com/anaconda/install/).

Installation procedure
----------------------

### intronSeeker installation and environment configuration.

To install intronSeeker, download all files and directories found at  
[URL](https://forgemia.inra.fr/faustine.oudin/Script_unigene/tree/modification_emilien),
then, open a bash session if you are not already in bash (bash command) and run the setup.sh script :

 ```diff
 + ./setup.sh
 ```

This script will  install the Conda environment, configure the newly installed
environment and install Grinder and its dependencies (Perl modules) wich are not yet
available in conda. Grinder installation is performed in the Conda environment, so, despite
the missing conda package, its installation will only affect your
conda environment (i.e. you will be able to run Grinder only when the environment is activated).
You should get these messages in the standard output :

```
Begin installing intronSeeker...

Directory added to the PATH

Conda environment installing...

Solving environment: done

Downloading and Extracting Packages
libcurl-7.63.0       | 550 KB    | ###################################################### | 100% 
[...]
perl-uri-1.74        | 54 KB     | ###################################################### | 100%
Preparing transaction: done
Verifying transaction: done
Executing transaction: done
#
# To activate this environment, use
#
#     $ conda activate ISeeker_environment
#
# To deactivate an active environment, use
#
#     $ conda deactivate

Conda environment installed.

Conda environment configuration...

Conda environment configured.

Installation done.

```

### Test installation.

NB. Installation can be quite long (15 minutes) due to the numerous grinder perl dependencies. 

Before each use of intronSeeker, activate the conda environment with :

```diff
+ source activate ISeeker_environment
```

Now, your command prompt should be like this :

```
(ISeeker_environment) elasguignes@node027 ~ $
```

Finally, to test the installation, run the command :

```diff
+ intronSeeker checkInstall
```

You should get this message in standard ouput : 

```
gffread testing...OK ! 

grinder testing...OK ! 

hisat2 testing...OK ! 

star testing...OK ! 

samtools testing...OK ! 

transdecoder testing...OK ! 

diamond testing...OK ! 

All the dependencies are correctly installed
```
