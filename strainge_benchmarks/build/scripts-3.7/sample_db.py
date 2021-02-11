"""
Sample database management for the StrainGE benchmarking pipeline
=================================================================

This is a helper script to manage lists of samples.
"""

import os
import sys
import logging
import argparse
from pathlib import Path

import pandas

import strainge_benchmarks as sb

EXTENSIONS = ["fq", "fastq", "fq.gz", "fastq.gz", "bam"]

logging.basicConfig()
logger = logging.getLogger()
logger.setLevel(logging.INFO)


def remove_sample_extension(name):
    return ".".join(
        p for p in name.split('.') if p not in {'bam', 'fastq', 'fq', 'gz'}
    )


def get_samples(filenames, remove_prefixes=None):
    for path in filenames:
        extension_checks = (path.name.endswith(ext) for ext in EXTENSIONS)
        if not any(extension_checks):
            logger.warning("Ignoring %s because it has an invalid extension."
                           " Allowed extensions: %s", path,
                           ",".join(EXTENSIONS))
            continue

        f1 = remove_sample_extension(path.name)

        # Check for read pairs
        if f1.endswith(".1"):
            # Split from the end of the string at ".1.", then add ".2.".
            # There's no `rreplace` in python.
            pair_path = path.parent / ".2.".join(path.name.rsplit(".1.", 1))

            if not pair_path.is_file():
                logger.warning("Ignoring %s, because the file with its mate "
                               "pairs doesn't exist. Checked for file: %s.",
                               path, pair_path)
                continue

            reads = [path, pair_path]
            sample_name = f1[:-2]
        elif f1.endswith("_R1"):
            # Split from the end of the string at "_R1", then add "_R2". Splits
            # only once at the end, so likely doesn't mess up anything else in
            # the filename.
            # There's no `rreplace` in python.
            pair_path = path.parent / "_R2".join(path.name.rsplit("_R1", 1))

            if not pair_path.is_file():
                logger.warning("Ignoring %s, because the file with its mate "
                               "pairs doesn't exist. Checked for file: %s.",
                               path, pair_path)
                continue

            reads = [path, pair_path]
            sample_name = f1[:-3]
        elif f1.endswith(".2") or f1.endswith("_R2"):
            # Ignore these files because we'll check them when we check the
            # corresponding .1 or _R1 file.
            continue
        else:
            reads = [path]
            sample_name = f1

        sample_name = str(path.parent / sample_name)
        if remove_prefixes:
            for prefix in remove_prefixes:
                if sample_name.startswith(prefix):
                    sample_name = sample_name[len(prefix):]
                    break

        parts = sample_name.split('/')
        if parts[-1] == parts[-2]:
            # Sometimes sample filename are contained in their own directory
            # with the same name, remove this redundancy in the sample name.
            sample_name = "/".join(parts[:-1])

        if len(reads) == 1:
            yield {
                'name': sample_name,
                'reads1': reads[0]
            }
        else:
            yield {
                'name': sample_name,
                'reads1': reads[0],
                'reads2': reads[1]
            }


def add_samples(database, files, fofn, remove_common_path, create, force=False,
                *args, **kwargs):
    if not force and create and database.is_file():
        logger.error("create: database file %s already exists!", database)
        return 1
    elif not create and not database.is_file():
        logger.error("append: database %s doesn't exist!", database)
        return 1

    if fofn:
        logger.info("Reading FOFN...")
        for f in fofn:
            files.append(Path(f.strip()))

    remove_prefixes = None
    if remove_common_path:
        common_path = os.path.commonpath(files) + '/'
        remove_prefixes = [common_path]

        logger.info("Common path to samples: %s, removing from ID.",
                    common_path)

    logger.info("Collecting samples...")
    samples = list(get_samples(files, remove_prefixes))
    new_samples_df = pandas.DataFrame(samples).set_index('name',
                                                         verify_integrity=True)

    if create:
        new_samples_df.to_csv(database, sep='\t')
    else:
        samples = sb.SampleList(database)
        samples.append(new_samples_df)

        samples.samples.to_csv(database, sep='\t')

    logger.info("Done. saved to %s.", database)


def main():
    parser = argparse.ArgumentParser(
        description="Manage lists of samples"
    )

    parser.add_argument(
        '-d', '--database', type=Path, required=True,
        help="Path to sample database CSV"
    )
    parser.set_defaults(func=None)

    subparsers = parser.add_subparsers()

    create_parser = subparsers.add_parser(
        'create', help="Create a new database based on a list of samples")
    create_parser.set_defaults(func=add_samples, create=True)

    create_parser.add_argument(
        'files', type=Path, nargs='*', metavar='FILE',
        help="List of sample files. This script will try to match files"
             "in case of paired-end reads."
    )

    create_parser.add_argument(
        '-f', '--fofn', type=argparse.FileType('r'), default=None,
        required=False,
        help="It's also possible to give a file of filenames where each line"
             "contains a path to a sample file. Use '-' to denote standard "
             "input."
    )

    create_parser.add_argument(
        '--force', action="store_true", required=False,
        help="Force creation of a new database even if the file exists."
    )

    create_parser.add_argument(
        '-c', '--remove-common-path', action="store_true", required=False,
        help="Each sample gets an ID based on its path. Use this option to "
             "remove a common path for all samples, before generating the ID. "
             "Example: if you have samples 'samples/t00.bam' and "
             "'samples/t01.bam', with this option enabled the sample ID's "
             "would be 't00' and 't01' respectively, because 'samples/' is "
             "a common path prefix."
    )

    append_parser = subparsers.add_parser(
        'append', help="Append new samples to an existing database.")
    append_parser.set_defaults(func=add_samples, create=False)

    append_parser.add_argument(
        'files', type=Path, nargs='*', metavar='FILE',
        help="List of sample files. This script will try to match files"
             "in case of paired-end reads."
    )

    append_parser.add_argument(
        '-f', '--fofn', type=argparse.FileType('r'), default=None,
        required=False,
        help="It's also possible to give a file of filenames where each line"
             "contains a path to a sample file. Use '-' to denote standard "
             "input."
    )

    append_parser.add_argument(
        '-c', '--remove-common-path', action="store_true", required=False,
        help="Each sample gets an ID based on its path. Use this option to "
             "remove a common path for all samples, before generating the ID. "
             "Example: if you have samples 'samples/t00.bam' and "
             "'samples/t01.bam', with this option enabled the sample ID's "
             "would be 't00' and 't01' respectively, because 'samples/' is "
             "a common path prefix."
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
