#!/usr/bin/env python3
"""Repair zero-coded allele 1 values in a PLINK .bim file."""

import argparse
from pathlib import Path


def repair_bim(input_path: Path, output_path: Path) -> int:
    """Write a repaired BIM file and return the number of repaired sites."""
    repaired_count = 0

    with input_path.open("r", encoding="utf-8") as infile, output_path.open(
        "w", encoding="utf-8"
    ) as outfile:
        for line_number, line in enumerate(infile, start=1):
            fields = line.split()

            if not fields:
                continue
            if len(fields) != 6:
                raise ValueError(
                    f"{input_path}:{line_number}: expected 6 columns, found {len(fields)}"
                )

            allele1, allele2 = fields[4], fields[5]

            if allele1 == "0" and allele2 != "0":
                fields[4] = allele2
                repaired_count += 1

            outfile.write("\t".join(fields) + "\n")

    return repaired_count


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Replace allele1=0 with allele2 when allele2 is nonzero."
    )
    parser.add_argument(
        "--input",
        required=True,
        type=Path,
        help="Input PLINK .bim file.",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output path for the repaired PLINK .bim file.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    repaired_count = repair_bim(args.input, args.output)
    print(f"Replaced allele1=0 at {repaired_count} site(s).")


if __name__ == "__main__":
    main()
