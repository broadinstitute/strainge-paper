"""
This module contains all sorts of helper functions to aid in the origanisation
of the snakemake benchmark workflow for StrainGE.
"""

import bz2
import gzip
import shlex
import functools
from pathlib import Path
from contextlib import contextmanager

import h5py
import pandas

config = None


def init(snakemake_config):
    global config
    config = snakemake_config


@contextmanager
def open_compressed(filename):
    if not isinstance(filename, Path):
        filename = Path(filename)

    if filename.suffix == ".gz":
        f = gzip.open(filename, "rt")
    elif filename.suffix == ".bz2":
        f = bz2.open(filename, "rt")
    else:
        f = open(filename)

    yield f

    f.close()


def get_config_key(*path, default=None):
    """Traverses the `config` dict hierarchy and tries to find
    the given key. If not found returns the given default value."""

    lookup_dict = config
    for key in path[:-1]:
        lookup_dict = lookup_dict.get(key, {})

    return lookup_dict.get(path[-1], default)


def get_resource_cfg(*path, default=None):
    return get_config_key('resources', *path, default=default)


def read_references_meta(tsv):
    return pandas.read_csv(tsv, sep='\t',
                           names=['ref_id', 'path', 'accession'],
                           index_col='ref_id')


@functools.lru_cache(maxsize=8)
def get_strainge_db_k(db_path):
    with h5py.File(db_path) as h5:
        return h5.attrs['k']


def bt2_index_files(path):
    ext = [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]
    return [f"{path}{e}" for e in ext]


class SampleList:
    """
    Helper class to manage all samples, where the reads are stored, and
    obtaining the right input files for a snakemake rule.
    """

    def __init__(self, sample_tsv):
        self.samples = pandas.read_csv(sample_tsv, sep='\t', comment='#',
                                       index_col='name')

    def append(self, new_samples_df):
        self.samples = pandas.concat([self.samples, new_samples_df],
                                     verify_integrity=True)

    def get_sample_input(self, sample):
        reads = [self.samples.loc[sample, 'reads1']]
        reads2 = self.samples.loc[sample, 'reads2']

        if reads2 and isinstance(reads2, str):
            reads.append(reads2)

        return reads

    def get_sample_refs(self, sample):
        sample_meta = self.samples.loc[sample]

        ref_coverages = []
        for i in range(1, 5):
            ref_key = f'reference{i}'
            cov_key = f'coverage{i}'

            if ref_key in sample_meta and cov_key in sample_meta:
                ref = sample_meta[ref_key]
                cov = sample_meta[cov_key]

                if isinstance(ref, str):
                    ref_coverages.append((ref, cov))

        return ref_coverages

    def get_sample_rng_seed(self, sample):
        return self.samples.loc[sample, 'seed']

    def get_target_ani(self, sample):
        if 'ani' in self.samples:
            ani = self.samples.loc[sample, 'ani']
            return ani if ani else 100.0

        return 100.0


def zcat_if_compressed(path):
    """Some programs don't accept compressed files. This function uses `zcat`
    and bash process substitution to return a filename to the uncompressed
    version."""

    if str(path).endswith(".gz"):
        return "<(zcat {})".format(shlex.quote(path))

    return path
