FROM continuumio/miniconda3:latest

#Replace conda by mamba
#To install and set the new solver, run the following commands:

RUN conda update -n base conda
RUN conda install -n base conda-libmamba-solver
RUN conda config --set solver libmamba  

# Copy intronSeeker's code in the futur Docker image:
COPY config/environment.yml  .
RUN ls /

# intronSeeker install:
RUN /bin/bash /intronSeeker/setup.sh
RUN source activate ISeeker_environment 
RUN intronSeeker checkInstall
RUN echo "End install and check installation..." 
