#!/usr/bin/env python3
"""Prepare multi-hit, reference-anchored FASTA inputs for tree-based orthology pruning."""
from pathlib import Path
from collections import defaultdict
import csv
import re
import sys

ROOT = Path('/mnt/sdb/Revised.ele')
OUT = ROOT / '04_multihit_orthology_screen'
PEP_DIR = ROOT / '01_inputs' / 'phylogenomics_69_pep'
CANDIDATES = OUT / '01_candidate_prefilter' / 'all_candidate_domain_coverage.tsv'
LOCUS_LIST = OUT / 'manifest' / 'SOG_5of5_1copy_478.txt'
SEED_DIR = OUT / 'seed_references'
RAW_DIR = OUT / '02_tree_candidate_fasta'

MIN_HMM_COVERAGE = 0.50
MAX_FULL_EVALUE = 1e-10
MIN_RELATIVE_BITSCORE = 0.50

# The original five strict single-copy seeds are preserved as fixed orthology anchors.
SEED_TO_OTU = {
    'Charabraunii': 'Charabraunii_cdhit',
    'Chlorokybusatmophyticus': 'Chlorokybusatmophyticus_cdhit',
    'ERR3487426': 'ERR3487426_cdhit',
    'Mesostigmaviride': 'Mesostigmaviride_cdhit',
    'SRR5949008': 'SRR5949008_cdhit',
}

def read_fasta(path):
    records = []
    header = None
    chunks = []
    with open(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    records.append((header, ''.join(chunks).replace('*', '')))
                header = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if header is not None:
            records.append((header, ''.join(chunks).replace('*', '')))
    return records

def safe_identifier(text):
    return re.sub(r'[^A-Za-z0-9_.-]+', '_', text).strip('_')

def as_float(value):
    return 0.0 if value == '0' else float(value)

loci = [line.strip() for line in LOCUS_LIST.read_text().splitlines() if line.strip()]
if len(loci) != 478 or len(set(loci)) != 478:
    raise SystemExit(f'Expected 478 unique seed loci, found {len(loci)}.')

# Read all reportable HMMER candidates and calculate the within-pair reference score.
rows = []
best_score = defaultdict(float)
with CANDIDATES.open() as handle:
    for row in csv.DictReader(handle, delimiter='\t'):
        score = float(row['full_sequence_bitscore'])
        key = (row['sample'], row['locus'])
        best_score[key] = max(best_score[key], score)
        row['_score'] = score
        row['_coverage'] = float(row['union_hmm_coverage'])
        row['_evalue'] = as_float(row['full_sequence_evalue'])
        rows.append(row)

selected = []
for row in rows:
    key = (row['sample'], row['locus'])
    if row['_coverage'] < MIN_HMM_COVERAGE:
        continue
    if row['_evalue'] > MAX_FULL_EVALUE:
        continue
    if row['_score'] < MIN_RELATIVE_BITSCORE * best_score[key]:
        continue
    selected.append(row)

# Group retained candidates by sample to fetch each PEP only once.
by_sample = defaultdict(list)
for row in selected:
    by_sample[row['sample']].append(row)

locus_records = defaultdict(list)
audit_rows = []
for sample, sample_rows in sorted(by_sample.items()):
    pep = PEP_DIR / f'{sample}.pep'
    if not pep.is_file():
        raise SystemExit(f'Missing active PEP file: {pep}')
    sequences = dict(read_fasta(pep))
    for row in sample_rows:
        target = row['target']
        if target not in sequences:
            raise SystemExit(f'Candidate {target} from {sample} was not found in {pep.name}')
        locus = row['locus']
        header = f"{sample}|CAND_{safe_identifier(target)}"
        locus_records[locus].append((header, sequences[target], row))

# Add the five original seed sequences to every locus under their active-taxonomy labels.
for locus in loci:
    seed_path = SEED_DIR / f'{locus}.fa'
    if not seed_path.is_file():
        raise SystemExit(f'Missing seed reference: {seed_path}')
    seed_records = read_fasta(seed_path)
    if len(seed_records) != 5:
        raise SystemExit(f'{seed_path.name} has {len(seed_records)} seed records, expected 5')
    for raw_header, sequence in seed_records:
        seed_name = raw_header.split('@', 1)[0]
        if seed_name not in SEED_TO_OTU:
            raise SystemExit(f'Unmapped seed label {seed_name} in {seed_path.name}')
        otu = SEED_TO_OTU[seed_name]
        identifier = raw_header.split('@', 1)[1] if '@' in raw_header else raw_header
        header = f"{otu}|SEED_{safe_identifier(identifier)}"
        locus_records[locus].append((header, sequence, {
            'sample': otu, 'locus': locus, 'target': raw_header,
            'target_length': len(sequence), 'hmm_length': '',
            'full_sequence_evalue': '', 'full_sequence_bitscore': '',
            'union_hmm_coverage': '', 'reported_domains': '',
            'source': 'fixed_seed_reference'
        }))

RAW_DIR.mkdir(parents=True, exist_ok=False)
per_locus = []
seen_audit = set()
for locus in loci:
    records = locus_records[locus]
    seen_headers = set()
    fasta_path = RAW_DIR / f'{locus}.fas'
    with fasta_path.open('w') as handle:
        for header, sequence, metadata in records:
            if not sequence:
                continue
            if header in seen_headers:
                # The seed record may also occur as a HMM candidate; retain both only
                # if their headers differ, otherwise avoid an identical duplicate header.
                continue
            seen_headers.add(header)
            handle.write(f'>{header}\n{sequence}\n')
            key = (locus, header)
            if key in seen_audit:
                continue
            seen_audit.add(key)
            source = metadata.get('source', 'hmm_candidate')
            audit_rows.append([
                locus, header, source, metadata['sample'], metadata['target'],
                metadata['target_length'], metadata['hmm_length'],
                metadata['full_sequence_evalue'], metadata['full_sequence_bitscore'],
                metadata['union_hmm_coverage'], metadata['reported_domains']
            ])
    total = len(seen_headers)
    candidate_headers = [h for h in seen_headers if '|CAND_' in h]
    candidate_otus = {h.split('|', 1)[0] for h in candidate_headers}
    per_locus.append([locus, total, len(candidate_headers), len(candidate_otus), 5])

with (RAW_DIR / 'candidate_sequence_audit.tsv').open('w', newline='') as handle:
    writer = csv.writer(handle, delimiter='\t')
    writer.writerow(['locus', 'tree_header', 'source', 'sample', 'original_target', 'target_length', 'hmm_length', 'full_sequence_evalue', 'full_sequence_bitscore', 'union_hmm_coverage', 'reported_domains'])
    writer.writerows(sorted(audit_rows))
with (RAW_DIR / 'per_locus_candidate_input_stats.tsv').open('w', newline='') as handle:
    writer = csv.writer(handle, delimiter='\t')
    writer.writerow(['locus', 'total_sequences_including_seed_references', 'hmm_candidates', 'active_taxa_with_hmm_candidate', 'fixed_seed_references'])
    writer.writerows(per_locus)
with (RAW_DIR / 'prefilter_parameters.tsv').open('w') as handle:
    handle.write('parameter\tvalue\n')
    handle.write(f'minimum_union_hmm_coverage\t{MIN_HMM_COVERAGE}\n')
    handle.write(f'maximum_full_sequence_evalue\t{MAX_FULL_EVALUE}\n')
    handle.write(f'minimum_relative_bitscore_within_sample_locus\t{MIN_RELATIVE_BITSCORE}\n')
    handle.write('selection_rule\tretain all candidates satisfying all three preliminary criteria; no single best hit is selected\n')
    handle.write('tree_pruning_plan\tuse reference-anchored, tree-based paralogy pruning after alignment and gene-tree inference\n')

print(f'loci_written={len(loci)}')
print(f'hmm_candidates_retained={len(selected)}')
print(f'total_tree_input_records={sum(x[1] for x in per_locus)}')
print(f'max_records_one_locus={max(x[1] for x in per_locus)}')
