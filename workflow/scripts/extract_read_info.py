#!/usr/bin/env python
"""extract_read_info.py: Extract read quality, read length, and haplotype from BAM."""
__author__ = "Adriana Sedeno"

import sys
import pysam


def process_bam(bam_file, output_file):
    bam = pysam.AlignmentFile(bam_file, "rb")

    with open(output_file, "w") as out:
        out.write("ReadName\tReadQuality\tReadLength\tHaplotype\n")

        for read in bam:
            rq = read.get_tag("rq") if read.has_tag("rq") else "NA"
            read_length = len(read.query_sequence)
            hp = read.get_tag("HP") if read.has_tag("HP") else "0"
            out.write(
                f"{read.query_name}\t{rq}\t{read_length}\t{hp}\n"
            )

    bam.close()
    print(f"Output written to {output_file}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <input.bam> <output.txt>")
        sys.exit(1)

    bam_file = sys.argv[1]
    output_file = sys.argv[2]
    process_bam(bam_file, output_file)
