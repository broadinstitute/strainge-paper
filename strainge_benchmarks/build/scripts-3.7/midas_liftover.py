import re
import sys
import logging
import argparse
from pathlib import Path
from datetime import datetime

import pysam
import pandas

from strainge_benchmarks import open_compressed

logging.basicConfig()
logger = logging.getLogger()
logger.setLevel(logging.INFO)

VCF_TEMPLATE = """\
##fileformat=VCFv4.1
##fileDate={date}
##source=midas_liftover.py
"""

cs_ops = re.compile(r"(:[0-9]+|\*[a-z][a-z]|[=\+\-][A-Za-z]+)")


def get_aligned_pairs_cs(segment):
    cs = segment.get_tag('cs')

    if not segment.query_sequence:
        return

    ref_pos = segment.reference_start
    query_pos = segment.query_alignment_start
    logger.debug("Ref start: %d, query start: %d", ref_pos, query_pos)
    for match in cs_ops.finditer(cs):
        op = match.group(1)

        if op[0] == ":":
            # Query and ref match
            num = int(op[1:])
            query_seq = segment.query_sequence[query_pos:query_pos+num]

            yield from zip(
                range(ref_pos, ref_pos+num),
                range(query_pos, query_pos+num),
                query_seq,
                query_seq
            )

            ref_pos += num
            query_pos += num
        elif op[0] == "*":
            # Substitution
            yield (ref_pos, query_pos, op[1].upper(), op[2].upper())

            ref_pos += 1
            query_pos += 1
        elif op[0] == "-":
            deletion = op[1:]
            del_len = len(deletion)

            yield from zip(
                range(ref_pos, ref_pos+del_len),
                [None] * del_len,
                deletion,
                [None] * del_len
            )

            ref_pos += del_len
        elif op[0] == "+":
            insertion = op[1:]
            ins_len = len(insertion)

            yield from zip(
                [None] * ins_len,
                range(query_pos, query_pos+ins_len),
                [None] * ins_len,
                insertion
            )

            query_pos += ins_len
        else:
            logger.warning("Unknown operation: %s", op)

    logger.debug("Ref end: %d, query end: %d", ref_pos, query_pos)


def main():
    parser = argparse.ArgumentParser(
        description="Liftover MIDAS SNP calls to a VCF based on a different "
                    "reference genome."
    )

    parser.add_argument(
        'ref_to_ref_aln', type=Path,
        help="The reference-to-reference whole genome alignment SAM/BAM. The "
             "target reference for the VCF should be the same as the reference"
             " used in this alignment file."
    )

    parser.add_argument(
        'midas_snps', type=Path,
        help="The output file from `run_midas.py snps` (for a single species)."
    )

    parser.add_argument(
        '-o', '--output', type=argparse.FileType('w'), default=sys.stdout,
        help="Output file, defaults to stdout."
    )

    args = parser.parse_args()

    logger.info("Loading SNPs called by MIDAS...")
    with open_compressed(args.midas_snps) as f:
        midas_snps = pandas.read_csv(f, index_col=['ref_id'],
                                     sep='\t')

        midas_snps = midas_snps[midas_snps['depth'] > 0]
        midas_snps['ref_pos'] -= 1
        midas_snps = midas_snps.reset_index().set_index(['ref_id', 'ref_pos'])

        logger.info("Determine major allele...")
        argmax = (midas_snps[['count_a', 'count_c', 'count_g', 'count_t']]
                  .idxmax(axis=1))
        midas_snps['major_allele'] = argmax.map({
            'count_a': 'A', 'count_c': 'C', 'count_g': 'G', 'count_t': 'T'})

        logger.info("%d positions with >0 reads", midas_snps.shape[0])

    bam_file = pysam.AlignmentFile(args.ref_to_ref_aln)

    vcf_tpl = VCF_TEMPLATE.format(
        date=datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    )
    vcf_header = pysam.VariantHeader()
    for line in vcf_tpl.split('\n'):
        line = line.strip()
        if not line:
            continue

        vcf_header.add_line(line)

    for scaffold in bam_file.references:
        vcf_header.contigs.add(scaffold,
                               length=bam_file.get_reference_length(scaffold))

    vcf_out = pysam.VariantFile(args.output, "w", header=vcf_header)

    untranslated = 0
    same_as_target = 0
    snp = 0

    for segment in bam_file.fetch():
        if segment.query_name not in midas_snps.index.levels[0]:
            logger.info("Skipping segment %s (no alleles present)",
                        segment.query_name)
            continue

        qstart = segment.query_alignment_start
        qend = segment.query_alignment_end

        logger.info("Alignment of %s (%d:%d) (supplementary: %s)",
                    segment.query_name, qstart, qend, segment.is_supplementary)

        segment_snps = midas_snps.loc[segment.query_name]

        pos_map = {
            qpos: (rpos, rbase, qbase)
            for rpos, qpos, rbase, qbase in get_aligned_pairs_cs(segment)
            if qpos is not None and rpos is not None
        }

        pos_set = set(p for p in segment_snps.index
                      if p >= qstart and p < qend)
        pos_translated = set()

        for pos in segment_snps.index:
            if pos < qstart or pos >= qend:
                continue

            if pos not in pos_map:
                continue

            ref_pos, ref_base, query_base = pos_map[pos]

            pos_translated.add(pos)

            if segment_snps.loc[pos, 'major_allele'] == ref_base:
                same_as_target += 1
                continue

            r = vcf_out.new_record()
            r.chrom = segment.reference_name
            r.pos = ref_pos + 1
            r.id = "."
            r.ref = ref_base
            r.alts = [segment_snps.loc[pos, 'major_allele']]
            r.qual = 50
            r.filter.add("PASS")

            vcf_out.write(r)
            snp += 1

        pos_untranslated = pos_set - pos_translated

        if pos_untranslated:
            logger.info("Could not translate %d positions on %s "
                        "due to missing alignment.", len(pos_untranslated),
                        segment.query_name)
            untranslated += len(pos_untranslated)

    logger.info("Total MIDAS positions with >0 reads: %d", midas_snps.shape[0])
    logger.info("Total SNPs on new reference %d/%d", snp, midas_snps.shape[0])
    logger.info("Total positions which match the target reference allele: "
                "%d/%d", same_as_target, midas_snps.shape[0])
    logger.info("Total positions on MIDAS reference that could not be "
                "translated to new ref: %d/%d", untranslated,
                midas_snps.shape[0])


if __name__ == '__main__':
    r = main()
    sys.exit(int(r) if r is not None else 0)
