How to Install ?
================

Hereafter, you will find all the details about the intronSeeker installation.

Requirements
------------

IntronSeeker requires Python version 3.6 or above. You don't need to change your Python version because this is managed by your Conda environment.

It also needs Python packages and external softwares to work 
correctly (all these dependancies and their versions are detailed in the file 
requirements.txt). So, to make install easier, a conda enviromnent (grinder excluded)
has been created  : Conda environment 
(it is the [environment.yml](https://forgemia.inra.fr/emilien.lasguignes/intronSeeker/-/blob/master/config/environment.yml) file)

For an easy install, conda environments is recommended. 
Follow [this tutorial](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) to install Miniconda.


## intronSeeker installation and environment configuration on genobioinfo
--------------------------------------------------------------------------

First, open a bash session if you are not already in bash (bash command).

### Clone intronSeeker code from Git.

Download all files and directories found at [URL](https://forgemia.inra.fr/emilien.lasguignes/intronSeeker) 

```diff
$ git clone https://forgemia.inra.fr/emilien.lasguignes/intronSeeker.git
Cloning into 'intronSeeker'...
remote: Enumerating objects: 2364, done.
remote: Counting objects: 100% (2334/2334), done.
remote: Compressing objects: 100% (798/798), done.
remote: Total 2364 (delta 1538), reused 2313 (delta 1523), pack-reused 30
Receiving objects: 100% (2364/2364), 539.80 MiB | 54.32 MiB/s, done.
Resolving deltas: 100% (1556/1556), done.
Checking out files: 100% (153/153), done.
```

### For a faster installation, we advise you to use mamba, available in the latest version of conda

```diff
$ module load devel/Miniconda/Miniconda3 
$ conda update -n base conda
$ conda install -n base conda-libmamba-solver
```

### Set up intronSeeker

```diff
$ cd intronSeeker/
$ ./setup.sh    
...
### Grinder Installed

### Installation done.
```

### Activate ISeeker_environment 

```diff
$ conda activate ISeeker_environment 
```

### If conda activate not found, initialize your shell

Message:
```diff
CommandNotFoundError: Your shell has not been properly configured to use 'conda activate'.
```

To initialize your shell, run

    $ conda init <SHELL_NAME>

Currently supported shells are:
  - bash
  - fish
  - tcsh
  - xonsh
  - zsh
  - powershell

==> For changes to take effect, close and re-open your current shell. <==

```diff
$ ssh -XY user@genobioinfo.toulouse.inrae.fr
(base) [user@genobioinfo1 intronSeeker]$ conda activate ISeeker_environment
```

### intronSeeker set up

Run the setup.sh script :

 ```diff
 $ /bin/bash setup.sh
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
$ source activate ISeeker_environment
```

Now, your command prompt should be like this :

```
(ISeeker_environment) [smaman@genobioinfo1 intronSeeker]$ 
```

Finally, to test the installation, run the command :

```diff
(ISeeker_environment) $ intronSeeker checkInstall
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

To see intronSeeker help:

```diff
(ISeeker_environment) $ intronSeeker -h

Program: intronSeeker
Version: 1.0

```
