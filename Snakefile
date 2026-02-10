from pathlib import Path
import os
import subprocess
import hashlib
import json
import shutil
import sys
import math

# --- PATCH Snakemake bug: Rule.__eq__ vs strings ---
# --- PATCH Snakemake bug: Rule.__eq__ vs strings ---
from snakemake.rules import Rule as _Rule

def _safe_rule_eq(self, other):
	if not isinstance(other, _Rule):
		return False
	return self.name == other.name and self.output == other.output

_Rule.__eq__ = _safe_rule_eq
# --- end patch ---

# Snakemake and working directories
SD = os.path.dirname(workflow.snakefile)


# Load config
with open("config.json") as f:
	config = json.load(f)
	
slurm        = config["slurm"]
QueryPath    = config["QueryPath"]
ScriptFolder = SD + "/scripts"
TargetFolder = config["TargetFolder"]
TempFolder   = "snaketemp"
blockmergesize = 80000

numcpu = 8

############################
# Pre-check for tools
############################
def precheck():
	required_tools = [
		"blastn", "makeblastdb",
		"winnowmap",
		"bedtools",
		"samtools",
		"minimap2"
	]
	missing = [tool for tool in required_tools if shutil.which(tool) is None]
	
	if missing:
		sys.stderr.write(
			f"\n[ERROR] Missing required tools: {', '.join(missing)}\n"
		)
		sys.stderr.write("Install them and ensure they are in PATH.\n")
		sys.exit(1)
	else:
		print("[CHECK] All required tools found in PATH.")
		
precheck()

def hotspot_done_for_prefix(wc):
        part = prefix_to_part[wc.prefix]
        return f"{TempFolder}/{part}.hotspots.done"


############################
# Make BLAST DBs
############################
def prepare(allfolders):
	for folder in allfolders:
		thetarget = f"{TargetFolder}/{folder}/{folder}_samples.fasta_fixed.fa_loci.txt.fasta_graph.FA"
		dbfile = thetarget + "_db.ndb"
		if not os.path.isfile(dbfile):
			os.system(
				f"bash {ScriptFolder}/runmakeblastdb -in {thetarget} -out {thetarget}_db"
			)
			
############################
# Read queries
############################

nlines = 0
with open(QueryPath, "r") as f:
	for i,line in enumerate(f):
		nlines += 1
		
numpart = math.ceil(nlines/numcpu)
allparts = [[] for x in range(numpart)]
queries = {}
prefix_to_part = {}
with open(QueryPath, "r") as f:
	for i,line in enumerate(f):
		line = line.split()
		qname, qfile = line[0], line[1]
		queries[qname] = qfile
		os.makedirs(f"{TempFolder}/{qname}", exist_ok=True)
		
		allparts[i % numpart].append(qname)
		prefix_to_part[qname] = (i % numpart)
		
for part in range(numpart):
	partfile = Path(f"{TempFolder}/{part}_hpart")
	
	if partfile.exists():
		continue
	
	with open(partfile, mode='w') as w: 
		for i,qname in enumerate(allparts[part]):
			w.write(f"{qname}\t{queries[qname]}\n")
	
	with open(str(partfile)+".output", mode='w') as w: 
		for i,qname in enumerate(allparts[part]):
			w.write(f"{TempFolder}/{qname}/{qname}\n")
			
############################
# Write all kmer file list
############################

Folders = os.path.join(TempFolder, "allfolders.list")
if not os.path.exists(Folders) or os.path.getsize(Folders) < 3:
	folders_l = []
	with os.scandir(TargetFolder) as it:
		for entry in it:
			if entry.is_dir() == False:
				continue
			name = entry.name
			kmers_path  = os.path.join(TargetFolder, name, f"{name}_samples.fasta_fixed.fa_qc.fa_kmers.txt")
			graph_path  = os.path.join(TargetFolder, name, f"{name}_samples.fasta_fixed.fa_loci.txt.fasta_graph.FA")
			if os.path.isfile(kmers_path) and os.path.isfile(graph_path):
				folders_l.append(entry.name)
	folders_l.sort()
	os.makedirs(os.path.dirname(Folders), exist_ok=True)
	with open(Folders, "w") as f:
		f.write("\n".join(folders_l) + ("\n" if folders_l else ""))
	print(f"Wrote {len(folders_l)} entries to {Folders}")
with open(Folders) as f:
	allfolders = [line.strip() for line in f if line.strip()]

Kmers = os.path.join(TempFolder, "allkmers.list")

if not os.path.exists(Kmers) or os.path.getsize(Kmers) < 3:
	allkmers = [
		f"{TargetFolder}/{name}/{name}_samples.fasta_fixed.fa_qc.fa_kmers.txt"
		for name in allfolders
	]
	with open(Kmers, "w") as f:
		for item in allkmers:
			f.write(item + "\n")
	print(f"Wrote {len(allkmers)} entries to {Kmers}")
	
############################
# Write all graph file list
############################
Targets = os.path.join(TempFolder, "allgraphs.list")

if not os.path.exists(Targets) or os.path.getsize(Targets) < 3:
	allgraphs = [
		f"{TargetFolder}/{name}/{name}_samples.fasta_fixed.fa_loci.txt.fasta_graph.FA"
		for name in allfolders
	]
	with open(Targets, "w") as f:
		for item in allgraphs:
			f.write(item + "\n")
	print(f"Wrote {len(allgraphs)} entries to {Targets}")

prepare(allfolders)


KmerBin = Kmers + ".bin"

############################
# RULE ALL  (fixed commas)
############################
rule all:
	input:
		KmerBin,
		hotspots=[f"{TempFolder}/{part}.hotspots.done" for part in range(numpart)],
		fastas=[f"{TempFolder}/{q}/{q}_fasta.txt" for q in queries]
		
		
############################
# KMER BIN
############################
rule kmer_bin:
	input:
		kmers=Kmers
	output:
		bin=KmerBin
	params:
		script=ScriptFolder
	resources:
		mem_mb=30000,
		slurm_extra="--mem=30G -c 1"
	threads: 1
	shell:
		"""
		touch temp.out
		{params.script}/kmer_searcher -T {input.kmers} -i temp.out -o temp2.out
		"""
		
		
############################
# HOTSPOT FINDING
############################
rule hotspots:
	input:
		part =TempFolder +"/{part}_hpart",
		bin = KmerBin
	output:
		marker = TempFolder + "/{part}.hotspots.done"
	params:
		script = ScriptFolder,
		kmers = Kmers
	resources:
		mem_mb = 30000,
		slurm_extra = f"--mem=30G -c {numcpu}"
	threads: numcpu
	run:
		print("running hotspots on", wildcards.part)
		shell("{params.script}/kmer_searcher -T {params.kmers} -I {input.part} -O {input.part}.output -c 100 -n {threads} -p \"\" ")
		shell("touch {output.marker} || true")
############################
# FASTA EXTRACTION
############################
rule fastas:
	input:
		done  = hotspot_done_for_prefix,
	output:
		text = TempFolder + "/{prefix}/{prefix}_fasta.txt"
	params:
		script=ScriptFolder,
		targets=Targets,
		temp=TempFolder
	resources:
		mem_mb=64000,
		slurm_extra=f"--mem=64G -c {numcpu}"
	threads: numcpu
	run:
		print("running fasta on", wildcards.prefix)
		queryfile = queries[wildcards.prefix]
		inputfile = params.temp+f"/{wildcards.prefix}/{wildcards.prefix}_hotspot.txt"
		shell(
			"python {params.script}/AllGraphAlign.py "
			"-i {inputfile} -g {params.targets} "
			"-s {wildcards.prefix} -q {queryfile} -o {output.text} "
			"-t {threads}"
		)
