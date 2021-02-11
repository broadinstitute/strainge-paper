"""
Script to generate configuration for benchmark data generation
==============================================================

"""

import sys
import random
import logging
import argparse
from pathlib import Path
from itertools import chain

import yaml
import h5py
import pandas
from snakemake.utils import update_config

import strainge_benchmarks as sb


STRAINGST_SYNTHETIC = {
    'num_mixes': 50,  # Per changed variable
    'single_strain_cov': [0.1, 1.0, 10.0, 100.0],
    'mixes_equal_cov': [0.1, 1.0, 10.0],
    'unequal_abundance_cov': [
        (100, 10),
        (100, 1.0),
        (100, 0.1),
        (10, 1.0),
        (10, 0.1),
        (1.0, 0.1)
    ]
}


STRAINGR_SYNTHETIC = {
    'num_mixes': 50,
    'single_strain_cov': [0.1, 1.0, 10.0],
    'mixes_equal_cov': [0.1, 1.0, 10.0],
    'unequal_abundance_cov': [
        (10, 1.0, 0.1),
        (10, 1.0),
        (1.0, 0.1)
    ],
    'ani': [99.99, 99.9, 99],
    'deletion_sizes': [1e3, 10e3, 100e3]
}


logging.basicConfig()
logger = logging.getLogger()
logger.setLevel(logging.INFO)


def single_strain_samples(num_mixes, coverages, references, db_refs):
    refs_not_in_db = references - db_refs

    num_in_db = int(num_mixes / 2)

    if refs_not_in_db:
        num_not_in_db = num_mixes - num_in_db
    else:
        num_not_in_db = 0

    random_in_db = random.sample(db_refs, num_in_db)

    if refs_not_in_db:
        random_not_in_db = random.sample(refs_not_in_db, num_not_in_db)
    else:
        random_not_in_db = []

    for cov in coverages:
        logger.info("Generating single strain samples at coverage %gx...", cov)
        for ref in random_in_db:
            sample_id = f"simulated/single/{ref}_{cov:g}x"

            # Each sample gets its own unique seed for reproducible ART
            # read generation
            seed = random.randrange(sys.maxsize)

            yield {
                'name': sample_id,
                'seed': seed,
                'reference1': ref,
                'coverage1': cov
            }

        for ref in random_not_in_db:
            sample_id = f"simulated/single/{ref}_{cov:g}x"

            # Each sample gets its own unique seed for reproducible ART
            # read generation
            seed = random.randrange(sys.maxsize)

            yield {
                'name': sample_id,
                'seed': seed,
                'reference1': ref,
                'coverage1': cov
            }


def equal_abundance_mixes(num_mixes, coverages, references, db_refs):
    for num in [2, 3, 4]:
        mixes = [random.sample(references, num) for i in range(num_mixes)]

        for cov in coverages:
            logger.info("Generating equal abundance mixes with %d strains at "
                        "coverage %gx", num, cov)
            for i, refs in enumerate(mixes):
                sample_id = f"simulated/mixes/equal{num}_{cov:g}x_{i}"

                reads_prefix = Path(
                    "benchmark_data_gen/" + sample_id).resolve()
                reads_prefix = str(reads_prefix / reads_prefix.parts[-1])

                # Each sample gets its own unique seed for reproducible ART
                # read generation
                seed = random.randrange(sys.maxsize)

                meta = {
                    'name': sample_id,
                    'seed': seed
                }

                for j, ref in enumerate(refs):
                    meta[f'reference{j+1}'] = ref
                    meta[f'coverage{j+1}'] = cov

                yield meta


def unequal_abundance_mixes(num_mixes, coverages, references, db_refs):
    # Per number of strains per mix, keep references constant so that we only
    # vary the coverage levels across mixes
    # So precompute which random references will be used
    mixes = {}
    for num_strains in set(len(cov) for cov in coverages):
        mixes[num_strains] = [random.sample(references, num_strains)
                              for i in range(num_mixes)]

    for cov in coverages:
        mix_list = mixes[len(cov)]

        logger.info("Generating mixes with unequal abundance, coverage: %s",
                    ":".join(f"{c:g}" for c in cov))

        for i, refs in enumerate(mix_list):
            cov_str = "_".join(f"{c:g}x" for c in cov)
            sample_id = f"simulated/mixes/unequal_{cov_str}_{i}"

            # Each sample gets its own unique seed for reproducible ART
            # read generation
            seed = random.randrange(sys.maxsize)

            meta = {
                'name': sample_id,
                'seed': seed
            }

            for j, ref in enumerate(refs):
                meta[f'reference{j+1}'] = ref
                meta[f'coverage{j+1}'] = cov[j]

            yield meta


def generate_titration_samples(sample_list, config):
    logger.info("Creating entries for titration samples...")

    for sample in sample_list:
        if not sample['name'].startswith("simulated"):
            continue

        meta = dict(sample)
        meta['name'] = meta['name'].replace(
            "simulated/", f"titration/")
        meta['seed'] = random.randrange(sys.maxsize)

        yield meta


def generate_straingst_samples(output, references_meta, strainge_db,
                               config=None, seed=None, *args, **kwargs):
    base_config = {}
    update_config(base_config, STRAINGR_SYNTHETIC)

    if config:
        with open(config) as f:
            config = yaml.load(f, Loader=yaml.CLoader)

        config = update_config(base_config, config)
    else:
        config = base_config

    references = sb.read_references_meta(references_meta)
    references = set(references.index)

    with h5py.File(strainge_db) as h5:
        db_refs = set()
        for key in h5:
            if isinstance(h5[key], h5py.Group):
                db_refs.add(key)

    if not seed:
        seed = random.randrange(sys.maxsize)

    random.seed(seed)
    logger.info("Initializing RNG with seed %d", seed)

    num_mixes = config['num_mixes']
    new_samples = list(chain(
        single_strain_samples(num_mixes, config['single_strain_cov'],
                              references, db_refs),
        equal_abundance_mixes(num_mixes, config['mixes_equal_cov'],
                              references, db_refs),
        unequal_abundance_mixes(num_mixes, config['unequal_abundance_cov'],
                                references, db_refs)
    ))

    titration_samples = generate_titration_samples(new_samples, config)
    new_samples.extend(titration_samples)

    new_df = pandas.DataFrame(new_samples).set_index('name',
                                                     verify_integrity=True)

    output.write(f"# Synthetic samples generated with seed {seed}\n")
    new_df.to_csv(output, sep='\t')


def generate_straingr_samples(output, references_meta, strainge_db,
                              config=None, seed=None, *args, **kwargs):
    base_config = {}
    update_config(base_config, STRAINGR_SYNTHETIC)

    if config:
        with open(config) as f:
            config = yaml.load(f, Loader=yaml.CLoader)

        config = update_config(base_config, config)
    else:
        config = base_config

    references = sb.read_references_meta(references_meta)
    references = set(references.index)

    with h5py.File(strainge_db) as h5:
        db_refs = set()
        for key in h5:
            if isinstance(h5[key], h5py.Group):
                db_refs.add(key)

    # Only keep references in the database
    references &= db_refs

    if not seed:
        seed = random.randrange(sys.maxsize)

    random.seed(seed)
    logger.info("Initializing RNG with seed %d", seed)

    num_mixes = config['num_mixes']
    new_samples = []
    for sample in list(chain(
        single_strain_samples(num_mixes, config['single_strain_cov'],
                              references, db_refs),
        equal_abundance_mixes(num_mixes, config['mixes_equal_cov'],
                              references, db_refs),
        unequal_abundance_mixes(num_mixes, config['unequal_abundance_cov'],
                                references, db_refs)
    )):
        # Mutate metadata to include target ANI
        for ani in config['ani']:
            meta = sample.copy()
            meta['name'] = sample['name'] + f"_ani_{ani:g}"
            meta['seed'] = random.randrange(sys.maxsize)
            meta['ani'] = ani

            new_samples.append(meta)

    titration_samples = generate_titration_samples(new_samples, config)
    new_samples.extend(titration_samples)

    new_df = pandas.DataFrame(new_samples).set_index('name',
                                                     verify_integrity=True)

    output.write(f"# Synthetic samples generated with seed {seed}\n")
    new_df.to_csv(output, sep='\t')


def main():
    parser = argparse.ArgumentParser(
        description="Manage lists of samples"
    )

    parser.add_argument(
        '-o', '--output', type=argparse.FileType('w'), default=sys.stdout,
        help="Output file"
    )

    parser.set_defaults(func=None)
    subparsers = parser.add_subparsers()

    gen_straingst_parser = subparsers.add_parser(
        'gen-straingst',
        help="Generate database entries for straingst benchmark samples."
    )
    gen_straingst_parser.set_defaults(func=generate_straingst_samples)

    gen_straingst_parser.add_argument(
        'references_meta',
        help="TSV with metadata on reference genomes (created by "
             "prepare_strainge_db.py)"
    )

    gen_straingst_parser.add_argument(
        'strainge_db',
        help="StrainGE database HDF5 file"
    )

    gen_straingst_parser.add_argument(
        '-c', '--config', type=Path, required=False, default=None,
        help="Optional YAML configuration file specifying parameters for "
             "straingst samples."
    )

    gen_straingst_parser.add_argument(
        '-s', '--seed', type=int, default=None, required=False,
        help="Seed for random number generator (optional)."
    )

    gen_straingr_parser = subparsers.add_parser(
        'gen-straingr',
        help="Generate database entries for straingr benchmark samples."
    )
    gen_straingr_parser.set_defaults(func=generate_straingr_samples)

    gen_straingr_parser.add_argument(
        'references_meta',
        help="TSV with metadata on reference genomes (created by "
             "prepare_strainge_db.py)"
    )

    gen_straingr_parser.add_argument(
        'strainge_db',
        help="StrainGE database HDF5 file"
    )

    gen_straingr_parser.add_argument(
        '-c', '--config', type=Path, required=False, default=None,
        help="Optional YAML configuration file specifying parameters for "
             "straingr samples."
    )

    gen_straingr_parser.add_argument(
        '-s', '--seed', type=int, default=None, required=False,
        help="Seed for random number generator (optional)."
    )

    args = parser.parse_args()

    if not args.func:
        print("No subcommand given, please see", sys.argv[0], "--help for "
              "more information")
        return 1

    status = args.func(**vars(args))

    return status if status else 0


if __name__ == '__main__':
    status = main()

    sys.exit(status)
