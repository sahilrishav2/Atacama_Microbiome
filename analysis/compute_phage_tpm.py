#!/usr/bin/env python3
"""Bacteriophage abundance (TPM) from DIAMOND blastx hits of sequencing reads against a viral protein database.

Method
  1. Each read pair is one fragment. A fragment is assigned to the phage of its best-scoring DIAMOND hit
     (best bit score over both mates); hits below --min-identity are ignored.
  2. A hit belongs to a phage when the organism name in square brackets at the end of the subject title
     (NCBI RefSeq style: ">YP_123.1 capsid protein [Escherichia phage T4]") matches --name-filter.
  3. For every phage, rate = fragments / (total length of all its proteins in the database, in kilo-amino-acids).
  4. TPM = rate / (sum of the rates of all phages in that sample) * 1e6, so each sample sums to 1e6 over all phages.
  5. Phages with fewer than --min-fragments fragments in total (all samples together) are left out of the output
     tables. TPM is computed BEFORE this filter, so the filter does not change any value.

Input files in --hits-dir, per sample: <sample>_R1.tsv and <sample>_R2.tsv, made with
  diamond blastx --outfmt 6 qseqid sseqid pident length evalue bitscore stitle --max-target-seqs 1
Output files in --outdir: bacteriophage_abundance_tpm.csv, bacteriophage_fragment_counts.csv,
  bacteriophage_sample_summary.tsv, bacteriophage_reference_lengths.tsv
Only the Python standard library is used.
"""
import argparse
import csv
import gzip
import os
import re
import sys
from collections import defaultdict

ORGANISM = re.compile(r"\[([^\[\]]+)\]\s*$")


def organism_of(title):
    """Organism name from the trailing [brackets] of a FASTA/DIAMOND title, or None."""
    match = ORGANISM.search(title.strip())
    return match.group(1).strip() if match else None


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_reference_lengths(faa, name_re):
    """organism -> [number of proteins, total amino acids], for organisms whose name matches name_re."""
    totals = {}
    current, aa = None, 0

    def commit():
        if current is not None and aa > 0:
            entry = totals.setdefault(current, [0, 0])
            entry[0] += 1
            entry[1] += aa

    with open_text(faa) as handle:
        for line in handle:
            if line.startswith(">"):
                commit()
                organism = organism_of(line[1:])
                current = organism if (organism and name_re.search(organism)) else None
                aa = 0
            else:
                aa += len(line.strip())
    commit()
    return totals


def count_fragments(hit_files, totals, min_identity):
    """Returns (fragments per phage, number of fragments with a viral hit) for one sample."""
    best = {}  # fragment id -> (bit score, organism, percent identity)
    for path in hit_files:
        with open_text(path) as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t", 6)
                if len(fields) < 7:
                    continue
                query, _subject, identity, _length, _evalue, bits, title = fields
                fragment = re.sub(r"/[12]$", "", query)
                score = float(bits)
                if fragment in best and best[fragment][0] >= score:
                    continue  # on ties the first mate (R1) is kept
                best[fragment] = (score, organism_of(title), float(identity))
    counts = defaultdict(int)
    n_viral = 0
    for _score, organism, identity in best.values():
        if identity < min_identity:
            continue
        n_viral += 1
        if organism in totals:
            counts[organism] += 1
    return counts, n_viral


def to_tpm(counts, totals):
    rates = {org: n / (totals[org][1] / 1000.0) for org, n in counts.items()}
    scale = sum(rates.values())
    return {org: rate / scale * 1e6 for org, rate in rates.items()} if scale > 0 else {}


def write_matrix(path, rows, samples, data):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, quoting=csv.QUOTE_NONNUMERIC, lineterminator="\n")
        writer.writerow([""] + samples)
        for row in rows:
            writer.writerow([row] + [data[s].get(row, 0) for s in samples])


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--hits-dir", required=True, help="folder with <sample>_R1.tsv / <sample>_R2.tsv")
    parser.add_argument("--reference-faa", required=True, help="protein FASTA the DIAMOND database was built from")
    parser.add_argument("--outdir", default="results/14_phage_abundance")
    parser.add_argument("--name-filter", default="phage",
                        help="case-insensitive regex an organism name must match to count as a phage [phage]")
    parser.add_argument("--min-identity", type=float, default=30.0, help="minimum percent identity of the best hit [30]")
    parser.add_argument("--min-fragments", type=int, default=10,
                        help="keep phages with at least this many fragments over all samples [10]")
    parser.add_argument("--samples", nargs="+", required=True, help="sample names, in the order of the output columns")
    args = parser.parse_args()

    name_re = re.compile(args.name_filter, re.IGNORECASE)
    print(f"Reading protein lengths from {args.reference_faa} (organisms matching '{args.name_filter}')", file=sys.stderr)
    totals = read_reference_lengths(args.reference_faa, name_re)
    if not totals:
        sys.exit(f"ERROR: no organism matching '{args.name_filter}' with a [organism] tag found in {args.reference_faa}")
    print(f"  {len(totals)} phage(s) in the reference", file=sys.stderr)

    samples, counts, tpms, summary = [], {}, {}, []
    for sample in args.samples:
        files = [p for p in (os.path.join(args.hits_dir, f"{sample}_R{m}.tsv") for m in (1, 2)) if os.path.exists(p)]
        if not files:
            print(f"WARNING: no DIAMOND output for {sample}, sample left out of the tables", file=sys.stderr)
            continue
        samples.append(sample)
        counts[sample], n_viral = count_fragments(files, totals, args.min_identity)
        tpms[sample] = to_tpm(counts[sample], totals)
        summary.append((sample, n_viral, sum(counts[sample].values()), len(counts[sample])))
    if not samples:
        sys.exit(f"ERROR: no DIAMOND output found in {args.hits_dir}")

    total_counts = defaultdict(int)
    for sample in samples:
        for organism, n in counts[sample].items():
            total_counts[organism] += n
    kept = [o for o, n in total_counts.items() if n >= args.min_fragments]
    if not kept:
        sys.exit(f"ERROR: no phage has at least {args.min_fragments} fragments in total "
                 f"({len(total_counts)} phage(s) detected). Lower --min-fragments or check the DIAMOND output.")
    mean_tpm = {o: sum(tpms[s].get(o, 0.0) for s in samples) / len(samples) for o in kept}
    kept.sort(key=lambda o: (-mean_tpm[o], o))

    os.makedirs(args.outdir, exist_ok=True)
    rounded = {s: {o: round(v, 6) for o, v in tpms[s].items()} for s in samples}
    write_matrix(os.path.join(args.outdir, "bacteriophage_abundance_tpm.csv"), kept, samples, rounded)
    write_matrix(os.path.join(args.outdir, "bacteriophage_fragment_counts.csv"), kept, samples, counts)
    with open(os.path.join(args.outdir, "bacteriophage_sample_summary.tsv"), "w") as handle:
        handle.write("sample\tfragments_with_viral_hit\tfragments_assigned_to_phages\tphages_detected\n")
        for row in summary:
            handle.write("\t".join(map(str, row)) + "\n")
    with open(os.path.join(args.outdir, "bacteriophage_reference_lengths.tsv"), "w") as handle:
        handle.write("phage\tproteins\ttotal_aa\n")
        for organism in kept:
            handle.write(f"{organism}\t{totals[organism][0]}\t{totals[organism][1]}\n")

    print(f"{len(samples)} sample(s), {len(kept)} phage(s) kept "
          f"(of {len(total_counts)} detected); tables written to {args.outdir}", file=sys.stderr)


if __name__ == "__main__":
    main()
