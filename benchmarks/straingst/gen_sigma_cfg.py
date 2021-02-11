"""
Make a configuration file for sigma
"""

import sys
import argparse
from pathlib import Path
from configparser import ConfigParser

SIGMA_BASE_CFG = Path("sigma/bin/sigma_config.cfg").resolve()
SIGMA_REF_DB = Path("sigma/ecoli").resolve()


def main():
    parser = argparse.ArgumentParser(
        description="Generate sample configuration to run Sigma"
    )

    parser.add_argument(
        '-r', '--reads', nargs=2, type=Path,
        help="Specify FASTQ files with reads"
    )

    parser.add_argument(
        '-o', '--output', type=argparse.FileType('w'), default=sys.stdout,
        help="Output file. Defaults to stdout."
    )

    args = parser.parse_args()

    sample_cfg = ConfigParser()
    sample_cfg.optionxform = str
    sample_cfg.read(SIGMA_BASE_CFG)

    sample_cfg['Data_Info']['Reference_Genome_Directory'] = str(SIGMA_REF_DB)

    name1 = args.reads[0].name.replace(".1.fq.gz", ".trimmed.1.fq")
    name2 = args.reads[1].name.replace(".2.fq.gz", ".trimmed.2.fq")
    sample_cfg['Data_Info']['Paired_End_Reads_1'] = str(
        args.reads[0].resolve().with_name(name1))
    sample_cfg['Data_Info']['Paired_End_Reads_2'] = str(
        args.reads[1].resolve().with_name(name2))

    sample_cfg.write(args.output)


if __name__ == '__main__':
    main()
