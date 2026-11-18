#!/usr/bin/env python3
"""Convert CDS GFF to POEM coordinates; optionally synchronize external protein IDs.

Standard Prodigal mode: -i gene.coords.gff -o gene_prediction_output
External mode: add -pf proteins.faa -af normalized_contigs.fa -map mapping.tsv
(External GFF seqids must match the ORIGINAL contig FASTA IDs.)
Supports single-interval CDS records; rejects ambiguous protein matching.
"""
import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path
from urllib.parse import quote, unquote
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def parse_attrs(value):
    attrs = {}
    for piece in value.strip().split(';'):
        if '=' in piece:
            key, val = piece.split('=', 1)
            attrs[key.strip()] = unquote(val.strip())
    return attrs


def read_mapping(path):
    if not path:
        return {}
    with open(path, newline='') as stream:
        rows = list(csv.DictReader(stream, delimiter='\t'))
    if not rows or not {'original_contig_id', 'poem_contig_id'}.issubset(rows[0]):
        raise ValueError('Invalid/empty contig ID mapping TSV')
    mapping = {row['original_contig_id']: row['poem_contig_id'] for row in rows}
    if len(mapping) != len(rows):
        raise ValueError('Duplicate contig IDs in mapping TSV')
    return mapping


def read_features(path, mapping, lengths):
    counts = defaultdict(int)
    features = []
    locations = set()
    with open(path) as stream:
        for line_no, line in enumerate(stream, 1):
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.rstrip('\n\r').split('\t')
            if len(fields) != 9:
                raise ValueError(f'GFF line {line_no}: expected 9 tab-separated fields')
            old, source, kind, start, end, score, strand, phase, attrs = fields
            if kind != 'CDS':
                continue
            if strand not in ('+', '-'):
                raise ValueError(f'GFF line {line_no}: unknown CDS strand {strand}')
            s, e = int(start), int(end)
            if s < 1 or e < s:
                raise ValueError(f'GFF line {line_no}: invalid CDS coordinates')
            new = mapping.get(old, old)
            if mapping and old not in mapping:
                raise ValueError(f'GFF line {line_no}: contig {old!r} not found in FASTA mapping')
            if lengths and (new not in lengths or e > lengths[new]):
                raise ValueError(f'GFF line {line_no}: CDS outside contig {new!r}')
            key = (new, s, e, strand)
            if key in locations:
                raise ValueError(f'GFF line {line_no}: duplicate CDS location {key}; multipart CDS not supported')
            locations.add(key)
            counts[new] += 1
            gene = f'{new}_{counts[new]}'
            attributes = parse_attrs(attrs)
            candidates = {attributes.get(k, '') for k in ('ID', 'protein_id', 'locus_tag', 'Name')}
            candidates |= {f'{old}_{counts[new]}', gene}
            features.append(dict(old=old, new=new, gene=gene, start=s, end=e,
                                 strand=strand, phase=phase, source=source,
                                 attrs=attributes, candidates=candidates - {''}))
    if not features:
        raise ValueError('No CDS features found in GFF')
    return features


def match_protein(feature, proteins, used):
    direct = [pid for pid in feature['candidates'] if pid in proteins and pid not in used]
    if len(direct) == 1:
        return direct[0]
    if len(direct) > 1:
        raise ValueError(f'Ambiguous protein IDs for {feature["gene"]}: {direct}')
    # Prodigal protein FASTA headers contain start/end/strand even if ID differs.
    matches = []
    for pid, record in proteins.items():
        if pid in used or not pid.startswith(feature['old'] + '_'):
            continue
        found = re.search(r'#\s*(\d+)\s*#\s*(\d+)\s*#\s*([+-]?1)\s*#', record.description)
        if found and (int(found.group(1)), int(found.group(2)), found.group(3)) == (
                feature['start'], feature['end'], '1' if feature['strand'] == '+' else '-1'):
            matches.append(pid)
    if len(matches) != 1:
        raise ValueError(f'Cannot uniquely match GFF CDS {feature["gene"]} with a protein. Matches: {matches}')
    return matches[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--input', required=True, help='CDS GFF file')
    parser.add_argument('-o', '--output', required=True, help='Gene prediction output directory')
    parser.add_argument('-pf', '--proteins', help='External protein FASTA for renaming')
    parser.add_argument('-af', '--contigs', help='Normalized nucleotide contig FASTA')
    parser.add_argument('-map', '--mapping', help='Original-to-normalized contig ID mapping')
    args = parser.parse_args()
    if args.proteins and (not args.contigs or not args.mapping):
        parser.error('External protein mode requires -pf, -af and -map')
    if not args.proteins and (args.contigs or args.mapping):
        parser.error('Use -af and -map only with -pf')

    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    mapping = read_mapping(args.mapping)
    contigs = SeqIO.to_dict(SeqIO.parse(args.contigs, 'fasta')) if args.contigs else {}
    features = read_features(args.input, mapping, {k: len(rec.seq) for k, rec in contigs.items()})
    proteins = {}
    used = set()
    if args.proteins:
        for record in SeqIO.parse(args.proteins, 'fasta'):
            if record.id in proteins:
                raise ValueError(f'Duplicate protein ID: {record.id}')
            proteins[record.id] = record
        if not proteins:
            raise ValueError('Empty protein FASTA')

    # Prevalidate all external matches before writing outputs.
    matches = []
    if proteins:
        for feat in features:
            pid = match_protein(feat, proteins, used)
            used.add(pid)
            matches.append(pid)
        if set(proteins) != used:
            raise ValueError(f'{len(set(proteins) - used)} proteins lack matching GFF CDS features')

    with (out / 'gene_coordinates.csv').open('w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['', 'contig_name', 'gene_name', 'start', 'end', 'strand'])
        for i, feat in enumerate(features):
            writer.writerow([i, feat['new'], feat['gene'], feat['start'], feat['end'], feat['strand']])

    if args.proteins:
        with (out / 'protein.translations.faa').open('w') as protein_out, \
             (out / 'dna.translations.faa').open('w') as dna_out, \
             (out / 'gene.coords.gff').open('w') as gff_out, \
             (out / 'gene_id_mapping.tsv').open('w', newline='') as id_map:
            writer = csv.writer(id_map, delimiter='\t')
            writer.writerow(['original_protein_id', 'poem_gene_id'])
            gff_out.write('##gff-version 3\n')
            for feat, pid in zip(features, matches):
                writer.writerow([pid, feat['gene']])
                SeqIO.write(SeqRecord(proteins[pid].seq, id=feat['gene'], description=''), protein_out, 'fasta')
                dna = contigs[feat['new']].seq[feat['start'] - 1:feat['end']]
                if feat['strand'] == '-':
                    dna = dna.reverse_complement()
                SeqIO.write(SeqRecord(dna, id=feat['gene'], description=''), dna_out, 'fasta')
                gff_out.write(f'{feat["new"]}\tPOEM-import\tCDS\t{feat["start"]}\t{feat["end"]}\t.\t{feat["strand"]}\t{feat["phase"]}\tID={quote(feat["gene"])}\n')
    print(f'Processed {len(features)} CDS features into {out / "gene_coordinates.csv"}')


if __name__ == '__main__':
    main()
