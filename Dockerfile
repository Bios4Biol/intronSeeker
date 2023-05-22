FROM continuumio/miniconda3:latest

#Replace conda by mamba
#To install and set the new solver, run the following commands:

RUN conda update -n base conda
RUN conda install -n base conda-libmamba-solver
RUN conda config --set solver libmamba  

# Copy intronSeeker's code in the futur Docker image:
COPY . /intronSeeker/

# intronSeeker install:
WORKDIR /intronSeeker/
RUN /bin/bash /intronSeeker/setup.sh

# source activate not functional, user SHELL tu run command
# https://pythonspeed.com/articles/activate-conda-dockerfile/
RUN ./entrypoint.sh intronSeeker checkInstall
RUN echo "End install and check installation..."

ENTRYPOINT ["./entrypoint.sh"]