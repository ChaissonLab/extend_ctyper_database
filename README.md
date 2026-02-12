This is for the one-time installation of the ctyper database on Azure Cyclecloud clusters.


##1. Preliminaries; cyclecloud setup.

The Azure CycleCloud Workspace for Slurm should be created as a resource.  This will
create the /shared NFS directory visible to all nodes.

The image should be Ubuntu 22.04, slurm 25.05.2 (current defaults with Cyclecloud)

The resources that we have used are:

Scheduler VM Type
Standard_D8as_v4
Login node VM Type
Standard_D8as_v4

HTC VM Type
Standard_E64-16s_v6, Standard_D8s_v6
The maximum number of nodes should scale with what you are comfortable running

The pipeline does not use HPC/GPU nodes, the maximum number of nodes for these 
should be set to 0.


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

###2.3 Installing the database construction scripts
First, change to the Extend_ori directory:
```cd /share/Extend_ori```

Next, clone this repository, and run the setup script:
```git clone https://github.com/ChaissonLab/extend_ctyper_database
cd extend_ctyper_database
source setup_env.sh
cp -r scripts ../
cp Snakefile ..```

Your environment to run database construction is complete.  Within the
 directory /share/Extend_ori, you should have a folder scripts, and a file Snakefile.

##3. Configuring database construction build.
The instructions on configuring files for building the database are available here:
https://github.com/Walfred-MA/PATs/

The only difference is that the `config.json` file that needs to be created
does not need to specify a path for scripts. The following can be used for
a `config.json` file to appropriately run Snakemake.

```
{
  "slurm": " ",
  "QueryPath": "query_pathes_withrefs.txt",
  "TargetFolder": "groups",
  "TempFolder": "snaketemp",
  "genelist": "./genes.bed",
  "ReferencePrefix": "NC_0609",
  "NumPartitions": 1
}
```
##4. Set up files.
Create a directory `/share/Extend_ori/Assemblies`, and copy all of the repeat masked assemblies that will be added to the database to this folder.
 The assemblies should have one FASTA file per haplotype as produced by hifiasm.  Assemblies produced with older assembly algorithms that produce
 primary and secondary contigs are *not* compatible with ctyper.


You will create a file  `query_pathes_withrefs.txt` (spelling needs to be exact, and to match what is in the .json file for QueryPath) that has two columns, one the name of the 

##5. Run database construction 
snakemake -k  --cluster "sbatch --partition htc  --time=500:00:00 {resources.slurm_extra}" --default-resources "mem_mb=3000" --jobs 500  --rerun-incomplete  --notemp --latency-wait 100 --resources mem_gb=1000

##6. Merge database files.
Once the snakemake run is done, the results need to be merged. You can follow step F here:
https://github.com/Walfred-MA/PATs/?tab=readme-ov-file#f-find-merge-and-index-compiled-k-mer-matrices
