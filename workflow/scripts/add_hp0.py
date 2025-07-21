#!/usr/bin/env python3

import sys
import pysam


def add_hp_tag(input_bam, output_bam):
    with pysam.AlignmentFile(input_bam, "rb") as in_bam:
        with pysam.AlignmentFile(output_bam, "wb", header=in_bam.header) as out_bam:
            for read in in_bam:
                if not read.has_tag("HP"):
                    read.set_tag("HP", 0, value_type="i")
                out_bam.write(read)
    pysam.index(output_bam)
    print("Successfully wrote {} with HP:i:0 tag and created index.".format(output_bam))


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: add_hp_tag.py <input.bam> <output.bam>")
        sys.exit(1)
    input_bam = sys.argv[1]
    output_bam = sys.argv[2]
    add_hp_tag(input_bam, output_bam)
