#!/usr/bin/env python3
import os
import subprocess
import argparse

# Adjust to your actual scripts path
script_folder = os.path.dirname(os.path.abspath(__file__))
def checkmask(fasta_file: str) -> int:
    bash_cmd = f'''
    grep -v '^>' "{fasta_file}" | tr -d '\\n' | tr -cd 'a-z' | wc -c
    '''
    result = subprocess.run(["bash", "-c", bash_cmd],
                            capture_output=True, text=True, check=True)
    return int(result.stdout.strip())

def callblastn(query: str, graphfile: str, output: str,
               nthreads: int, simi: float = 0.90, ifrepeat: int = 1) -> None:
    if ifrepeat:
        repeatopts = "-dust yes -lcase_masking"
    else:
        repeatopts = ""

    cmd = (
        "bash {}/runblastn -task megablast -query {} -db {} "
        "-gapopen 10 -gapextend 2 -word_size 30 -perc_identity {} {} "
        "-evalue 1e-200 -outfmt 17 -out {} -num_threads {} -max_target_seqs 100"
    ).format(
        script_folder,
        query,
        graphfile + "_db",
        int(simi * 100),
        repeatopts,
        output,
        nthreads,
    )
    os.system(cmd)

def callalign(query: str, graphfile: str, output: str,
              nthreads: int, simi: float = 0.90,
              ifhighqual: int = 1, ifrepeat: int = 1) -> None:
    if ifhighqual:
        callblastn(query, graphfile, output, nthreads, simi, ifrepeat)

        if (not os.path.exists(output) or os.path.getsize(output) <= 10):
            print(f"WARNING: alignment fail :{output}, realigning\n")
            callblastn(query, graphfile, output, nthreads, simi, ifrepeat)

        if ifrepeat:
            cmd = "{}/runwinnowmap.sh {} {} {} {} 300 > {}".format(
                script_folder, query, graphfile, nthreads, simi, output + "_wm"
            )
            os.system(cmd)
            if (not os.path.exists(output + "_wm") or os.path.getsize(output + "_wm") <= 10):
                os.system(cmd)
            os.system(f"( cat {output}_wm >> {output} && rm {output}_wm ) || true")
    else:
        cmd = "{}/runwinnowmap.sh {} {} {} {} 150 masking > {}".format(
            script_folder, query, graphfile, nthreads, simi, output
        )
        os.system(cmd)

def run_align(graphfile: str, thefile: str, simi: float,
              ifhighqual: int = 1) -> None:
    lowcasecount = checkmask(thefile)

    ifrepeat = 0
    if lowcasecount > 5000:
        ifrepeat = 1

    output = thefile + "_align.out"
    callalign(thefile, graphfile, output, 1, simi, ifhighqual, ifrepeat)

def main():
    parser = argparse.ArgumentParser(
        description="Run alignment for a single FASTA file"
    )
    parser.add_argument("fasta", help="Path to query fasta (.fa) file")
    parser.add_argument(
        "--graphfile",
        default="",
        help="Graph FASTA file (without _db suffix)",
    )
    parser.add_argument(
        "--simi", type=float, default=0.9,
        help="Similarity threshold (default: 0.9)",
    )
    parser.add_argument(
        "--ifhighqual", type=int, default=1,
        help="ifhighqual flag (0/1), default 0 for this batch",
    )
    args = parser.parse_args()

    print(f"Running alignment for {args.fasta}")
    run_align(args.graphfile, args.fasta, args.simi, ifhighqual=args.ifhighqual)

if __name__ == "__main__":
    main()

