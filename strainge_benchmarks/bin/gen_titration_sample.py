#!/usr/bin/env python3
"""
Mix reads from one sample randomly with reads from another sample, at a given
sequencing depth and relative abundance.
"""

import sys
import random
import logging
import argparse
import tempfile
import subprocess
from pathlib import Path

from Bio import SeqIO
from sh import seqtk

from strainge_benchmarks import open_compressed

logging.basicConfig()
logger = logging.getLogger()
logger.setLevel(logging.INFO)


SHUFFLE_SHELL = """
paste  <(cat {r1} "{b1}") \\
       <(cat {r2} "{b2}") \\
    | paste - - - - \\
    | shuf \\
    | awk -F'\\t' '{{OFS="\\n"; print $1,$3,$5,$7 > "{out_prefix}.1.fq";\
        print $2,$4,$6,$8 > "{out_prefix}.2.fq"}}'

gzip "{out_prefix}.1.fq"
gzip "{out_prefix}.2.fq"
"""


def main():
    parser = argparse.ArgumentParser(
        description="Generate artificial titration (spike-in) samples."
    )

    parser.add_argument(
        '-R1', type=Path,
        help="FastQ reads file 1"
    )

    parser.add_argument(
        '-R2', type=Path,
        help="FastQ reads file 2"
    )

    parser.add_argument(
        '-B1', type=Path,
        help="Background sample FASTQ reads file 1"
    )

    parser.add_argument(
        '-B2', type=Path,
        help="Background sample FASTQ reads file 2"
    )

    parser.add_argument(
        '-d', '--depth', type=int,
        help="Depth of sequencing"
    )

    parser.add_argument(
        '-s', '--seed', type=int, default=None, required=False,
        help="Random number generator seed."
    )

    parser.add_argument(
        '-o', '--output-prefix',
        help="Output filename prefix"
    )

    args = parser.parse_args()

    seed = args.seed if args.seed else random.randrange(1, 2**32 - 1)
    logger.info("Random seed: %d", seed)

    # Determine read length for both samples
    in_sample_read_len = None
    with open_compressed(args.R1) as f:
        read = next(SeqIO.parse(f, "fastq"))
        in_sample_read_len = len(read.seq)

    base_sample_read_len = None
    with open_compressed(args.B1) as f:
        read = next(SeqIO.parse(f, "fastq"))
        base_sample_read_len = len(read.seq)

    if in_sample_read_len != base_sample_read_len:
        logger.error("Read length of input sample and base sample is not the "
                     "same (%d != %d)", in_sample_read_len,
                     base_sample_read_len)

        return 1

    # Determine number of reads to sample
    total_read_pairs = int(args.depth / (in_sample_read_len * 2))

    logger.info("Counting number of reads in input sample...")
    with open_compressed(args.R1) as f:
        p = subprocess.run(f'grep "^@" | wc -l',
                           stdin=f, stdout=subprocess.PIPE,
                           shell=True, executable='/bin/bash')

    num_reads_in = int(p.stdout)
    reads_from_bg = total_read_pairs - num_reads_in

    logger.info("Total number of read pairs required: %d", total_read_pairs)
    logger.info("Number of read pairs in input sample: %d", num_reads_in)
    logger.info("Number of read pairs from background sample: %d",
                reads_from_bg)

    logger.info("Subsampling reads from background sample...")
    b1 = tempfile.NamedTemporaryFile(mode='w')
    b2 = tempfile.NamedTemporaryFile(mode='w')
    seqtk.sample("-s", seed, args.B1, reads_from_bg,
                 _out=b1)
    seqtk.sample("-s", seed, args.B2, reads_from_bg,
                 _out=b2)

    if args.R1.name.endswith('.gz'):
        r1 = f'<(zcat "{args.R1}")'
        r2 = f'<(zcat "{args.R2}")'
    else:
        r1 = f'"{args.R1}"'
        r2 = f'"{args.R2}"'

    logger.info("Shuffling reads and gzipping output files...")
    subprocess.run(
        SHUFFLE_SHELL.format(b1=b1.name, b2=b2.name, r1=r1, r2=r2,
                             out_prefix=args.output_prefix),
        shell=True, executable='/bin/bash'
    )

    b1.close()
    b2.close()


if __name__ == '__main__':
    r = main()
    sys.exit(r if r else 0)
