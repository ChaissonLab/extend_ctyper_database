import argparse
import os
import shutil
from pathlib import Path

def Summary(input_file, output_dir):
    with open(input_file, 'r') as infile:
        for line in infile:
            line = line.strip()
            if not line:
                continue
            
            fa_path = Path(line).resolve()
            prefix = fa_path.stem.split('_')[0]
            target_folder = Path(output_dir) / prefix
            target_folder.mkdir(parents=True, exist_ok=True)

            # Copy original .FA file
            shutil.copy2(fa_path, target_folder / fa_path.name)

            # Copy .cache.json file
            cache_file = fa_path.with_name(fa_path.name.replace('.FA', 'cache.json'))
            if cache_file.exists():
                shutil.copy2(cache_file, target_folder / cache_file.name)
            else:
                print(f"[WARNING] Missing cache file: {cache_file}")

            # Copy ._qc.fa_kmers.txt file
            qc_file = fa_path.with_name(fa_path.name.replace("_loci.txt.fasta_graph.FA", "_qc.fa_kmers.txt"))
            if qc_file.exists():
                shutil.copy2(qc_file, target_folder / qc_file.name)
            else:
                print(f"[WARNING] Missing qc kmer file: {qc_file}")

            # Handle partitions
            original_folder = fa_path.parent
            partitions_dir = original_folder / 'partitions'
            if partitions_dir.exists() and partitions_dir.is_dir():
                for subdir in partitions_dir.iterdir():
                    if subdir.is_dir():
                        new_subdir = target_folder / 'partitions' / subdir.name
                        new_subdir.mkdir(parents=True, exist_ok=True)
                        kmer_file = subdir / f"{prefix}{subdir.name}_samples.fasta_kmer.list"
                        if kmer_file.exists():
                            shutil.copy2(kmer_file, new_subdir / kmer_file.name)
                        else:
                            print(f"[WARNING] Missing kmer list: {kmer_file}")
            else:
                print(f"[WARNING] No partitions directory in {original_folder}")

def main(args):
    Summary(args.input, args.output)

def run():
    """
    Parse arguments and run
    """
    parser = argparse.ArgumentParser(description="Copy graph and cache files into organized folders.")
    parser.add_argument("-i", "--input", help="path to input data file", dest="input", type=str, required=True)
    parser.add_argument("-o", "--output", help="path to output folder", dest="output", type=str, required=True)

    parser.set_defaults(func=main)
    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    run()

