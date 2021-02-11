"""
Artificially introduce small variants and bigger deletions into a known
reference genome.
"""

import random
import logging
import argparse
from pathlib import Path
from datetime import datetime

import pysam
import numpy
import pandas
from scipy.stats import lognorm
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

from strainge_benchmarks import open_compressed

logging.basicConfig()
logger = logging.getLogger()
logger.setLevel(logging.INFO)

VCF_TEMPLATE = """\
##fileformat=VCFv4.1
##fileDate={date}
##source=mutate_ref.py
##reference={reference}
"""


def open_fasta_file(filename: Path):
    """Open an optionally compressed FASTA file"""

    with open_compressed(filename) as f:
        yield from SeqIO.parse(f, "fasta")


def parse_gff_atrr(attr_str):
    attr = {}

    for definition in attr_str.split(';'):
        k, v = definition.split('=', maxsplit=1)
        k = k.strip()
        v = v.strip()

        attr[k] = v

    return attr


def read_gene_annotations(gff: Path):
    genes = []
    with open_compressed(gff) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            parts = line.split('\t')

            if parts[2] != 'gene':
                continue

            genes.append({
                'chrom': parts[0],
                'start': int(parts[3]) - 1,
                'end': int(parts[4]),
                'strand': parts[6],
                'orig_str': line.strip(),
                **parse_gff_atrr(parts[8])
            })

    return pandas.DataFrame(genes).set_index('chrom')


def mutate_ref(ref: Path, genes: Path, similarity, snp_weight, ins_weight,
               del_weight, gene_deletions, seed=None):
    logger.info("Target ANI: %g", similarity)

    weights_sum = snp_weight + ins_weight + del_weight
    snp_prob = snp_weight / weights_sum
    ins_prob = ins_weight / weights_sum
    del_prob = del_weight / weights_sum
    prob_vector = [snp_prob, ins_prob, del_prob]

    gene_annots = None
    if gene_deletions and genes:
        gene_annots = read_gene_annotations(genes)
        logger.info("Read gene annotations, total number of genes: %d",
                    gene_annots.shape[0])
    elif not genes:
        logger.warning("No gene annotations file given, skipping gene "
                       "deletions.")

    if seed:
        numpy.random.seed(seed)

    BASES = set("ACTG")
    for scaffold in open_fasta_file(ref):
        length = len(scaffold.seq)
        sequence = str(scaffold.seq)
        pos_adjust = numpy.zeros((length,), dtype=numpy.int32)

        num_mutations = int(round(
            length * ((100.0 - similarity) / 100)))
        mean_dist = int(round(length / num_mutations))
        logger.info("Processing scaffold %s... length: %d", scaffold.id,
                    length)
        logger.info("Generating %d mutations, average distance between"
                    "mutations: %d.", num_mutations, mean_dist)

        mutation_deltas = lognorm.rvs(scale=mean_dist, s=0.3,
                                      size=num_mutations, random_state=seed)
        mutation_deltas = numpy.round(mutation_deltas)

        old_ref_pos = 0
        cur_ref_pos = 0
        new_ref = []
        mutations = []
        for delta in mutation_deltas:
            old_ref_pos = cur_ref_pos
            cur_ref_pos += int(delta)

            if cur_ref_pos >= length:
                break

            ref_base = sequence[cur_ref_pos-1]

            mut_type = numpy.random.choice(list(range(3)), p=prob_vector)
            if mut_type == 0:
                # SNP
                bases = tuple(BASES - {ref_base})
                snp_base = random.choice(bases)

                new_ref.append(sequence[old_ref_pos:cur_ref_pos-1])
                new_ref.append(snp_base)

                mutations.append({
                    'chrom': scaffold.id,
                    'pos': cur_ref_pos - 1,
                    'ref': ref_base,
                    'alt': snp_base
                })
            elif mut_type == 1:
                # Insertion
                new_base = random.choice(tuple(BASES))
                new_ref.append(sequence[old_ref_pos:cur_ref_pos])
                new_ref.append(new_base)

                mutations.append({
                    'chrom': scaffold.id,
                    'pos': cur_ref_pos - 1,
                    'ref': ref_base,
                    'alt': ref_base + new_base
                })

                pos_adjust[cur_ref_pos:] += 1
            else:
                # Deletion
                new_ref.append(sequence[old_ref_pos:cur_ref_pos-1])

                mutations.append({
                    'chrom': scaffold.id,
                    'pos': cur_ref_pos - 2,
                    'ref': sequence[cur_ref_pos-2] + ref_base,
                    'alt': sequence[cur_ref_pos-2]
                })

                pos_adjust[cur_ref_pos:] -= 1

        new_ref = "".join(new_ref)
        mutations = pandas.DataFrame(mutations)
        genes_to_delete = None

        logger.info("%d mutations actually generated", mutations.shape[0])

        if gene_annots is not None:
            scaffold_genes = (gene_annots.loc[scaffold.id].copy()
                              .reset_index()
                              .set_index(['start', 'end']))

            genes_to_delete = scaffold_genes.sample(
                frac=gene_deletions / 100)
            logger.info("Deleting %g%% genes (total: %d)...", gene_deletions,
                        genes_to_delete.shape[0])

            for start, end in genes_to_delete.index:
                adj_start = start + pos_adjust[start]
                adj_end = end + pos_adjust[end]

                # Cut it out
                new_ref = new_ref[:adj_start] + new_ref[adj_end:]
                pos_adjust[end:] -= end - start

                ix = (mutations['pos'] >= start) & (mutations['pos'] < end)
                if mutations[ix].shape[0] > 0:
                    logger.info("%d mutations fall in a gene to be deleted.",
                                mutations[ix].shape[0])

                mutations = mutations[~ix].copy()

        scaffold.seq = Seq(new_ref, IUPAC.extended_dna)
        scaffold.description += " (mutated)"

        yield scaffold, mutations, genes_to_delete


def main():
    parser = argparse.ArgumentParser(
        description="Introduce SNPs, indels and large deletions in a given "
                    "reference genome."
    )

    parser.add_argument(
        'reference', type=Path,
        help="Path to reference genome"
    )

    parser.add_argument(
        'genes', type=Path, nargs='?', default=None,
        help="Optional path to a GFF file with genes for this reference."
    )

    parser.add_argument(
        '--seed', type=int, required=False, default=None,
        help="Random number generator seed."
    )

    parser.add_argument(
        '-o', '--output-prefix', required=True,
        help="Output filename prefix, this script generates both a FASTA and "
             "a VCF file."
    )

    mut_group = parser.add_argument_group(
        'Mutation generation settings',
        description="Configure how many and what kind of mutations are "
                    "generated. The mutation type is determined at random "
                    "proportional to the weights given. Defaults to 80% "
                    "SNPs, 10% insertions, and 10% deletions."
    )

    mut_group.add_argument(
        '-s', '--similarity', type=float, required=True,
        help="Target average nucleotide identity to the original reference."
    )

    mut_group.add_argument(
        '--snp-weight', type=int, default=80, required=False,
        help="Weight for generating a SNP. Defaults to %(default)d."
    )

    mut_group.add_argument(
        '--ins-weight', type=int, default=10, required=False,
        help="Weight for generating an insertion. Defaults to %(default)d."
    )
    mut_group.add_argument(
        '--del-weight', type=int, default=10, required=False,
        help="Weight for generating a deletion. Defaults to %(default)d."
    )

    mut_group.add_argument(
        '-d', '--gene-deletions', type=float, default=0.0, required=False,
        help="Percentage of genes to delete, requires a gene annotation GFF"
             "file."
    )

    args = parser.parse_args()

    Path(args.output_prefix).parent.mkdir(parents=True, exist_ok=True)

    out_fasta = open(args.output_prefix + ".fa", "w")
    out_vcf = args.output_prefix + ".mutations.vcf"
    out_gff = open(args.output_prefix + ".genes_deleted.gff", "w")

    ref = args.reference.resolve()
    vcf_tpl = VCF_TEMPLATE.format(
        date=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        reference=f"file://{ref}"
    )
    vcf_header = pysam.VariantHeader()
    for line in vcf_tpl.split('\n'):
        line = line.strip()
        if not line:
            continue

        vcf_header.add_line(line)

    for scaffold in open_fasta_file(ref):
        vcf_header.contigs.add(scaffold.id, length=len(scaffold.seq))

    vcf_out = pysam.VariantFile(out_vcf, "w", header=vcf_header)

    for new_scaffold, mutations, gene_deletions in mutate_ref(
        args.reference, args.genes, args.similarity, args.snp_weight,
        args.ins_weight, args.del_weight, args.gene_deletions, args.seed
    ):
        SeqIO.write(new_scaffold, out_fasta, "fasta")

        for ix in mutations.index:
            r = vcf_out.new_record()
            r.chrom = mutations.loc[ix, 'chrom']
            r.pos = mutations.loc[ix, 'pos'] + 1
            r.id = "."
            r.ref = mutations.loc[ix, 'ref']
            r.alts = [mutations.loc[ix, 'alt']]
            r.qual = 50
            r.filter.add("PASS")

            vcf_out.write(r)

        if gene_deletions is not None:
            for gene in gene_deletions['orig_str']:
                out_gff.write(f"{gene}\n")

    out_fasta.close()
    vcf_out.close()
    out_gff.close()


if __name__ == '__main__':
    main()
