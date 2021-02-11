import sys
import argparse
from pathlib import Path

from Bio import SeqIO

from strainge_benchmarks import open_compressed


def main():
    parser = argparse.ArgumentParser(
        description='Create a big map which links contigs to reference genomes'
    )

    parser.add_argument(
        'ref', nargs='+',
        help="Reference genomes."
    )

    parser.add_argument(
        '-o', '--output', default=sys.stdout, type=argparse.FileType('w'),
        help="Output file."
    )

    args = parser.parse_args()

    print("ref", "contig", sep='\t', file=args.output)
    for ref in args.ref:
        ref_path = Path(ref)
        ref_name = Path(ref_path.name.replace(".gz", "")).with_suffix("").name

        with open_compressed(ref) as f:
            for record in SeqIO.parse(f, "fasta"):
                print(ref_name, record.id, sep='\t', file=args.output)


if __name__ == '__main__':
    main()
