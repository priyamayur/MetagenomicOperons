#!/usr/bin/env python3
"""Normalize assembled nucleotide FASTA IDs, optionally writing an old-to-new map."""
import argparse
import csv
import os
import tempfile
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('-af', '--assembledfile', required=True)
    p.add_argument('-o', '--output', required=True, help='Output FASTA file')
    p.add_argument('-map', '--mapping', help='Optional original-to-normalized TSV')
    args = p.parse_args()
    target = Path(args.output)
    target.parent.mkdir(parents=True, exist_ok=True)
    if args.mapping and Path(args.mapping).resolve() == target.resolve():
        p.error('Mapping and FASTA output paths must differ')
    seen = set()
    mapping = []
    fd, temp = tempfile.mkstemp(prefix='.poem_contigs_', dir=target.parent)
    try:
        with os.fdopen(fd, 'w') as out:
            for number, record in enumerate(SeqIO.parse(args.assembledfile, 'fasta'), 1):
                if record.id in seen:
                    raise ValueError(f'Duplicate contig ID: {record.id}')
                seen.add(record.id)
                normalized = f'contig_{number}'
                mapping.append((record.id, normalized))
                SeqIO.write(SeqRecord(record.seq, id=normalized, description=''), out, 'fasta')
        if not mapping:
            raise ValueError('No contigs found in the input FASTA')
        if args.mapping:
            mapping_target = Path(args.mapping)
            mapping_target.parent.mkdir(parents=True, exist_ok=True)
            with mapping_target.open('w', newline='') as fh:
                writer = csv.writer(fh, delimiter='\t')
                writer.writerow(['original_contig_id', 'poem_contig_id'])
                writer.writerows(mapping)
        os.replace(temp, target)
    finally:
        if os.path.exists(temp):
            os.remove(temp)
    print(f'Normalized {len(mapping)} contig FASTA headers: {target}')


if __name__ == '__main__':
    main()
