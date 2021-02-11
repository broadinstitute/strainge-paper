"""
Mask regions with high density of variant calls. This is to try to give MIDAS
a bit of slack in our evaluation where lots and lots of variants where called
in certain regions (likely due to its crappy E. coli assembly).
"""

import sys
import argparse

import numpy
import pysam
import pandas

from strainge.io.variants import boolean_array_to_bedfile


def main():
    parser = argparse.ArgumentParser(
        description="Mask MIDAS high density SNP regions"
    )
    parser.add_argument('midas_vcf', help="Midas SNP calls as VCF")
    parser.add_argument(
        '-o', '--output', type=argparse.FileType('w'), default=sys.stdout,
    )

    args = parser.parse_args()

    vcf = pysam.VariantFile(args.midas_vcf)
    variants_per_contig = {}
    for contig in vcf.header.contigs:
        length = vcf.header.contigs[contig].length
        variants_per_contig[contig] = numpy.zeros((length,))

    for record in vcf.fetch():
        variants_per_contig[record.contig][record.pos-1] = 1

    for contig, variants in variants_per_contig.items():
        rolling = pandas.Series(variants).rolling(int(1e3),
                                                  min_periods=1).sum()
        high_density = rolling.values > 50

        boolean_array_to_bedfile(~high_density, args.output, contig, 25)


if __name__ == '__main__':
    main()
