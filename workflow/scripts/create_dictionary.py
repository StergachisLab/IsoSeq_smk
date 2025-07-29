#!/usr/bin/env python3

import argparse
import sqlite3
import gzip
import os
import re


def parse_args():
    parser = argparse.ArgumentParser(description="Efficiently create dictionary from large files using SQLite indexing for Snakemake.")
    parser.add_argument("collapsed_file", help="Path to collapsed file")
    parser.add_argument("haplotag_db", help="Path to SQLite DB containing haplotags")
    parser.add_argument("classification_db", help="Path to SQLite DB containing classification data")
    parser.add_argument("output_file", help="Path for output gzipped TSV")
    return parser.parse_args()


def process_files(collapsed_file, haplotag_db, classification_db, output_file):
    conn_class = sqlite3.connect(classification_db)
    cur_class = conn_class.cursor()

    conn_hap = sqlite3.connect(haplotag_db)
    cur_hap = conn_hap.cursor()

    pattern_sample = re.compile(r'(PS\d+)')
    pattern_condition = re.compile(r'(treated|untreated)')

    with open(collapsed_file) as f, gzip.open(output_file, 'wt') as gz_out:
        gz_out.write('\t'.join([
            'read_id', 'isoform', 'structural_category', 'associated_gene',
            'associated_transcript', 'subcategory', 'haplotype', 'sample', 'condition'
        ]) + '\n')

        for line in f:
            parts = line.strip().split('\t')
            read_id, isoform = parts[0], parts[1]
            sample_match = pattern_sample.search(read_id)
            condition_match = pattern_condition.search(read_id)

            sample = sample_match.group(1) if sample_match else "NA"
            condition = condition_match.group(1) if condition_match else "NA"

            cur_class.execute('SELECT structural_category, associated_gene, associated_transcript, subcategory FROM classification WHERE isoform=?', (isoform,))
            class_row = cur_class.fetchone()
            if class_row is None:
                class_row = ('unassigned', 'NA', 'NA', 'unassigned')

            cur_hap.execute('SELECT haplotype FROM haplotags WHERE read_id=?', (read_id,))
            hap_row = cur_hap.fetchone()
            haplotype = hap_row[0] if hap_row else 'NA'

            gz_out.write('\t'.join([
                read_id,
                isoform,
                class_row[0],
                class_row[1],
                class_row[2],
                class_row[3],
                haplotype,
                sample,
                condition
            ]) + '\n')

    conn_class.close()
    conn_hap.close()


def main():
    args = parse_args()
    process_files(args.collapsed_file, args.haplotag_db, args.classification_db, args.output_file)


if __name__ == '__main__':
    main()
