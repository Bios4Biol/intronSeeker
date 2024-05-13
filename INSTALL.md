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

For an easy install, creating a conda environment is recommended. 
Follow [this tutorial](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) to install Miniconda.
For a faster installation, we advise you to use mamba, already available in the latest version of conda


## intronSeeker installation and environment configuration on genobioinfo
--------------------------------------------------------------------------

First, open a bash session if you are not already in bash (bash command).

### Clone intronSeeker code from Git.

Download all files and directories found at [URL](https://forgemia.inra.fr/emilien.lasguignes/intronSeeker) 

```diff
$ git clone https://github.com/Bios4Biol/intronSeeker.git
Cloning into 'intronSeeker'...
remote: Enumerating objects: 2364, done.
remote: Counting objects: 100% (2334/2334), done.
remote: Compressing objects: 100% (798/798), done.
remote: Total 2364 (delta 1538), reused 2313 (delta 1523), pack-reused 30
Receiving objects: 100% (2364/2364), 539.80 MiB | 54.32 MiB/s, done.
Resolving deltas: 100% (1556/1556), done.
Checking out files: 100% (153/153), done.
```

### Install intronSeeker environment.


```diff
$ cd intronSeeker
$ conda env create -f config/environment.yml
```

You can speedup intronSeeker installation by using [`mamba`](https://mamba.readthedocs.io/en/latest/) or [`libmamba`](https://conda.github.io/conda-libmamba-solver/)


### Activate ISeeker_environment and test installation.

Before each use of intronSeeker, activate the conda environment with :

```diff
$ conda activate ISeeker_environment
```

Now, your command prompt should be like this :

```
(ISeeker_environment) [...]$ 
```

Finally, to test the installation, run the command :

```diff
(ISeeker_environment) $ ./intronSeeker checkInstall
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
(ISeeker_environment) $ ./intronSeeker -h

Program: intronSeeker
Version: 1.0

```

## Add intronSeeker into you conda env PATH (once)

intronSeeker can be run from anywhere. 
In order to avoid to type the full PATH to intronSeeker program, you can create a symbolic link into the `ISeeker_environment` environment.
You can use the following command when `ISeeker_environment` environment is loaded to do this:

```diff
ln -s "${PWD}"/intronSeeker "${CONDA_PREFIX}/bin/"
```

Then you can run it from anywhere:


```diff
(ISeeker_environment) $ intronSeeker -h

Program: intronSeeker
Version: 1.0

```