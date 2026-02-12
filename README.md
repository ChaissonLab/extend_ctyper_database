This is for the one-time installation of the ctyper database.


##1. Preliminaries; cyclecloud setup.

The Azure CycleCloud Workspace for Slurm should be created as a resource.  This will create the /shared NFS directory visible to all nodes.

The image should be Ubuntu 22.04, slurm 25.05.2 (current defaults with Cyclecloud)

The resources that we have used are:

Scheduler VM Type
Standard_D8as_v4
Login node VM Type
Standard_D8as_v4

HTC VM Type
Standard_E64-16s_v6, Standard_D8s_v6
The maximum number of nodes should scale with what you are comfortable running

The pipeline does not use HPC/GPU nodes, the maximum number of nodes for these should be set to 0.
HPC VM Type
GPU VM Type

The Extend_ori directory should be copied to /shared

##2. Environment setup.
###2.1 g++
Most of the environment is set up using conda, however some binaries are built from source,
which requires g++ to be installed.


```
sudo apt update
sudo apt install -y build-essential git
```

###2.2 conda

You will want to install conda in a location that is accessible on NFS.

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x ./Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh
```
Follow the instructions, and ensure that the root directory is under `/share`
