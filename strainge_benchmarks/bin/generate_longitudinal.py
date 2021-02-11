"""
Generate synthetic longitudinal dataset
"""

import argparse
from pathlib import Path

from sh import art_illumina


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-r', '--references-meta', type=argparse.FileType('r'),
        help="Path to reference genome metadata sheet."
    )

    parser.add_argument(
        'specification', type=argparse.FileType('r'),
        help="Dataset specification CSV"
    )

    parser.add_argument(
        '-o', '--output-dir', type=Path,
        help="Output directory"
    )

    args = parser.parse_args()

    if not args.output_dir.is_dir():
        args.output_dir.mkdir(parents=True, exist_ok=True)



def




if __name__ == '__main__':
    main()
