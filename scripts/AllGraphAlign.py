#!/usr/bin/env python3

from __future__ import annotations

import argparse
import collections as cl
import pathlib
import subprocess
import os
from typing import Dict, List, Tuple
import multiprocessing as mp
import time
import signal

from LocalGraphAlign import GetRegions

BLAST_IDENTITY = 90
BLAST_THREADS  = 1
OVERLAP_GAP    = 30000
EXTENDSIZE = 30000

# globals for workers
QUERYSEQ = None
OUTPUT_PATH = None
OUT_LOCK = None
LOG_LOCK = None

def reverse_complement(seq: str) -> str:
	trans_table = str.maketrans("ATCGatcg", "TAGCtagc")
	return seq.translate(trans_table)[::-1]


def sh(cmd: str, *, check: bool = False) -> None:
	subprocess.run(cmd, shell=True, check=check)
	
	
def read_fasta(path: pathlib.Path) -> Dict[str, str]:
	seqs: Dict[str, str] = {}
	with path.open(mode = 'r') as fh:
		head, buf = "", []
		for ln in fh:
			ln = ln.rstrip()
			if ln.startswith(">"):
				if head:
					seqs[head] = "".join(buf)
				head, buf = ln[1:].split()[0], []
			else:
				buf.append(ln)
		if head:
			seqs[head] = "".join(buf)
	return seqs

def merge_intervals(spans: List[Tuple[int, int]], gap: int = OVERLAP_GAP) -> List[Tuple[int, int]]:
	"""Greedy merge of (start,end) pairs allowing ≤ *gap* bp separation."""
	if not spans:
		return []
	
	spans.sort()
	merged = [list(spans[0])]
	for s, e in spans[1:]:
		last = merged[-1]
		if s - last[1] <= gap:
			last[1] = max(last[1], e)
		else:
			merged.append([s, e])
			
	return [(max(0, s - gap // 2), e + gap // 2) for s, e in merged]


def output_regions(results: List[List], align_lines: List[str], output: str):
	# use global lock
	with OUT_LOCK:
		with open(output, mode="a") as w:
			for result in results:
				gindex, prefix, chrom, strd, start, end = result
				w.write(f"{gindex}\t{prefix}\t{chrom}\t{strd}\t{start}\t{end}\n")
		if align_lines:
			align_out = output + "_align.out.txt"
			with open(align_out, mode="a") as w2:
				for result, meta in zip(results, align_lines):
					gindex, prefix, chrom, strd, start, end = result
					# keep original columns, then add combined first-string
					w2.write(
						f"{gindex}\t{prefix}\t{chrom}\t{strd}\t{start}\t{end}\t{meta}\n"
					)
					
					
					
def process_block(
	gindex: int,
	hotspot: Dict[str, List[List[int]]],
	sample: str,
	target: str,
	nthread: int = 2
	):
	# use globals in worker
	query = QUERYSEQ
	output = OUTPUT_PATH
	
	prefix = target.split("/")[-1].split("_")[0]
	
	folder = os.path.dirname(target)
	os.makedirs(os.path.join(folder, "samples"), exist_ok=True)
	hotspot_fa = f"{folder}/samples/{sample}_hotspot.txt.fa"
	
	regions = cl.defaultdict(list)
	for chrom, region in hotspot.items():
		regions[chrom] = merge_intervals(region)
		
	allresults = []
	index = 1
	
	alignresults = []
	t0 = time.perf_counter()
	for contig, regionlist in regions.items():
		seq_ref = query[contig]
		for beg, end in regionlist:
			
			seq = seq_ref[beg:end]
			# write query fasta for this hotspot
			with open(hotspot_fa, "w") as out:
				out.write(f">{sample}_{index}\t{contig}:{beg}-{end}+\n")
				out.write(seq + "\n")
			index += 1
			lowercase = sum(1 for c in seq if c.islower())
			outfile = hotspot_fa + "_align.out"
			print(f"working on file:  {hotspot_fa}, {contig}:{beg}-{end} on {target}")
			# NOTE: adjust args to GetRegions to match your actual signature
			results = GetRegions(hotspot_fa, target, outfile, lowercase, nthread)
			
			newresults = []
			for result in results:
				if result[0] == '+':
					newresults.append([result[0], beg + result[1], beg + result[2]])
				else:
					newresults.append([result[0], end - 1 - result[2], end - 1 - result[1]])
					
			index_start = len(allresults)
			if len(newresults):
				index_start = len(allresults)
				allresults.extend([[gindex,prefix, contig] + x for x in newresults])
				index_end = len(allresults)
				# NEW: read first lines from *_align.out.txt and hotspot_fa
				align_txt_path = hotspot_fa + "_align.out.txt"
				firststring = ""
				secondstring = ""
				
				with open(align_txt_path, "r") as f1:
					secondstring = f1.readline().rstrip("\n")
					
				with open(hotspot_fa, "r") as f2:
					firststring = f2.readline().rstrip("\n")
					
				meta = firststring + "\t" + secondstring +"\t" + f"{index_start}_{index_end}\n"
				# extend results and alignresults in parallel
				alignresults.append(meta)
				
				
			#Path(hotspot_fa).unlink(missing_ok=True)
			#Path(hotspot_fa + "_align.out").unlink(missing_ok=True)
				
	t1 = time.perf_counter()
	
	with LOG_LOCK:
		print(f"[TimeLog:] {prefix} processing took {t1 - t0:.2f} seconds")
	if allresults:
		output_regions(allresults, alignresults, output)
		

def _process_block_entry(q: mp.Queue, args_tuple):
	"""
	Run process_block(*args_tuple) in a dedicated process group so the parent
	can kill the whole group on timeout.
	"""
	try:
		# Put this worker in its own process group (Linux/Unix).
		# Parent can then kill the group via os.killpg(pid, ...).
		os.setsid()
	except Exception:
		# If setsid is not available/allowed, we still proceed.
		pass
		
	try:
		process_block(*args_tuple)
		q.put(("ok", None))
	except Exception as e:
		q.put(("err", repr(e)))
		
def run_tasks_with_restarts(tasks, nthreads, timeout_sec=3600, max_retries=3):
	"""
	Run tasks with:
		- bounded parallelism
		- retries on exceptions
		- NO retry on timeout; timeout => kill + skip

	Each task runs in its own mp.Process so we can hard-kill it.
	"""
	processes = max(1, int(nthreads // 2))
	
	# retry counter per task index
	retries = {i: 0 for i in range(len(tasks))}
	remaining = cl.deque(range(len(tasks)))  # task indices to schedule
	running = {}  # i -> dict(proc=..., q=..., start=..., args=...)
	
	timed_out = set()
	
	def _kill_task(proc: mp.Process):
		"""Best-effort kill of a task process and its process group."""
		try:
			# If child did os.setsid(), its pid is a process-group leader.
			os.killpg(proc.pid, signal.SIGKILL)
		except Exception:
			# Fallback: kill just the process
			try:
				proc.terminate()
			except Exception:
				pass
				
		try:
			proc.join(timeout=2)
		except Exception:
			pass
			
		if proc.is_alive():
			try:
				proc.kill()
			except Exception:
				pass
			try:
				proc.join(timeout=2)
			except Exception:
				pass
				
	with LOG_LOCK:
		print(f"[Info] Running {len(tasks)} tasks with up to {processes} parallel workers; "
				f"timeout={timeout_sec}s; max_retries={max_retries}")
		
	while remaining or running:
		# launch new tasks if slots available
		while remaining and len(running) < processes:
			i = remaining.popleft()
			if i in timed_out:
				continue
			
			# Optionally adjust per-task threads as you did before
			gindex, hotspot, sample, target, _ = tasks[i]
			threads = max(1, int(nthreads // 2))
			args_tuple = (gindex, hotspot, sample, target, threads)
			
			q = mp.Queue(maxsize=1)
			proc = mp.Process(target=_process_block_entry, args=(q, args_tuple))
			proc.daemon = False
			proc.start()
			
			running[i] = {"proc": proc, "q": q, "start": time.time(), "args": args_tuple}
			
		# monitor running tasks
		now = time.time()
		finished = []
		
		for i, info in list(running.items()):
			proc: mp.Process = info["proc"]
			q: mp.Queue = info["q"]
			t0 = info["start"]
			
			# timeout check
			if proc.is_alive() and (now - t0) > timeout_sec:
				gindex, hotspot, sample, target, threads = info["args"]
				with LOG_LOCK:
					print(f"[Timeout] graph={target} (gindex={gindex}) exceeded {timeout_sec}s, killing and skipping (no retry)")
				_kill_task(proc)
				retries[i] += 1
				if retries[i] < max_retries:
					with LOG_LOCK:
						print(f"[Timeout-Retry] graph={target} (gindex={gindex}), retry {retries[i]}/{max_retries}")
					finished.append(i)
					remaining.append(i)
				else:
					with LOG_LOCK:
						print(f"[Timeout-GiveUp] graph={target} (gindex={gindex}), timed out {retries[i]} times; giving up")
					timed_out.add(i)
					finished.append(i)
				continue
			
			# completion check
			
			if not proc.is_alive():
				proc.join(timeout=0.1)
				
				# Default to error unless we explicitly receive "ok"
				status = ("err", "no status returned")
				try:
					status = q.get_nowait()
				except Exception:
					status = ("err", "no status returned")
					
				if status[0] == "ok":
					finished.append(i)
				else:
					retries[i] += 1
					if retries[i] < max_retries:
						with LOG_LOCK:
							print(f"[Warning] task {i} failed with {status[1]}, "
									f"retry {retries[i]}/{max_retries}")
						remaining.append(i)
					else:
						with LOG_LOCK:
							print(f"[Error] task {i} failed with {status[1]} "
									f"after {max_retries} attempts; giving up")
					finished.append(i)
					
		for i in finished:
			# cleanup
			info = running.pop(i, None)
			if info is not None:
				try:
					info["q"].close()
				except Exception:
					pass
					
		if not finished:
			time.sleep(0.05)
			
	with LOG_LOCK:
		ok_cnt = len(tasks) - len(timed_out)
		print(f"[Info] All tasks finished. Timed out tasks (skipped): {len(timed_out)}. "
				f"Non-timeout completed/failed: {ok_cnt}.")

		
def cli() -> argparse.Namespace:
	p = argparse.ArgumentParser(description="Hotspot → local graph → regions")
	p.add_argument("-i", "--input",   required=True, help="hotspot txt (per locus)")
	p.add_argument("-q", "--query",   required=True, help="reference FASTA")
	p.add_argument("-g", "--graphs",  required=True, help="file listing graph targets, one per line")
	p.add_argument("-s", "--samplename",  required=True, help="sample name")
	p.add_argument("-o", "--output",  required=True, help="output BED")
	p.add_argument("-t", "--threads",   type=int, default=4, help="number of threads")
	p.add_argument("--timeout", type=int, default=3600, help="per-task timeout in seconds (default: 3600)")
	return p.parse_args()

def main() -> None:
	
	# <<< CHANGED: Force 'fork' so QUERYSEQ is inherited and shared (Linux only)
	try:
		mp.set_start_method("fork")
	except RuntimeError:
		# Start method already set (e.g., if this is run from another context)
		pass
		
	global QUERYSEQ, OUTPUT_PATH, OUT_LOCK, LOG_LOCK
	
	args = cli()
	
	# 1) targets
	targets = []
	with open(args.graphs, mode = 'r') as f:
		for line in f:
			line = line.strip()
			if line:
				targets.append(line)
				
	# 2) reference (load ONCE in parent)
	QUERYSEQ = read_fasta(pathlib.Path(args.query))
	
	# 3) hotspots input -> build structure:
	#    hotspots[gindex][contig] = [[beg, end], ...]
	hotspots = cl.defaultdict(lambda: cl.defaultdict(list))
	with open(args.input, mode = 'r') as f:
		for line in f:
			parts = line.strip().split()
			contig = parts[0]
			beg0 = int(parts[-2])
			end0 = int(parts[-1])
			
			beg = max(0, beg0 - EXTENDSIZE)
			end = min(len(QUERYSEQ[contig]), end0 + EXTENDSIZE)
			
			gindex = int(parts[1])  # you did this in your code
			hotspots[gindex][contig].append([beg, end])
			
	# 4) prepare output
	OUTPUT_PATH = args.output
	open(OUTPUT_PATH, "w").close()
	
	# 5) shared lock
	OUT_LOCK = mp.Lock()
	LOG_LOCK = mp.Lock()
	# 6) run in parallel
	tasks = []
	for gindex, hotspot in hotspots.items():
		# gindex are 1-based, target file is 0-based in your code
		target = targets[gindex - 1]
		tasks.append((gindex, hotspot, args.samplename, target, 2))
		
	run_tasks_with_restarts(tasks, args.threads, timeout_sec=args.timeout)


	
def _init_worker(output_path, lock, lock2):
	"""
	Worker initializer.

	NOTE: QUERYSEQ is NOT passed here; with 'fork' start method
	workers inherit the already-loaded QUERYSEQ from the parent,
	so the large dict is shared via copy-on-write and not copied.
	"""
	global OUTPUT_PATH, OUT_LOCK, LOG_LOCK
	OUTPUT_PATH = output_path
	OUT_LOCK = lock
	LOG_LOCK = lock2
	
	
if __name__ == "__main__":
	
	main()
