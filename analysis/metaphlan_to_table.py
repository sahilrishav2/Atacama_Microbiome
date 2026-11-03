#!/usr/bin/env python3
"""Converts the merged MetaPhlAn table of step 3 into a taxa x samples CSV for the downstream R Markdown documents
(Microbial_community_analysis.Rmd, Co-occurrence_network.Rmd).

  * one rank is taken from the MetaPhlAn table (--rank, default species); rows are named after the taxon
    without the rank prefix (s__Cutibacterium_acnes -> Cutibacterium_acnes)
  * the values are MetaPhlAn relative abundances (percent), NOT read counts
  * sample names are cleaned (--sample-suffix is removed from the end of the names)
  * --keep-samples-from METADATA.csv keeps only the samples listed in that file and puts the columns in the same
    order, which the documents require for their metadata files (column --sample-column, default: first column)

Output: CSV with the taxa in the first column and one column per sample, in the same quoting style as R's write.csv.
Only the Python standard library is used.
"""
import argparse
import csv
import re
import sys

PREFIX = {"kingdom": "k__", "phylum": "p__", "class": "c__", "order": "o__", "family": "f__",
          "genus": "g__", "species": "s__"}
NON_SAMPLE = {"clade_name", "NCBI_tax_id", "clade_taxid", "additional_species"}


def read_metaphlan(path):
    with open(path) as handle:
        lines = handle.read().splitlines()
    header = next((i for i, line in enumerate(lines) if re.match(r"^#?clade_name", line)), None)
    if header is None:
        sys.exit(f"ERROR: no 'clade_name' header line in {path}; is this a merged MetaPhlAn table?")
    columns = lines[header].lstrip("#").split("\t")
    keep = [i for i, name in enumerate(columns) if name not in NON_SAMPLE]
    rows = []
    for line in lines[header + 1:]:
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < len(columns):
            continue
        rows.append((fields[0], [fields[i] for i in keep]))
    return [columns[i] for i in keep], rows


def sample_order(path, column):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle) if column else csv.reader(handle)
        if column:
            if column not in reader.fieldnames:
                sys.exit(f"ERROR: column '{column}' not found in {path} (columns: {', '.join(reader.fieldnames)})")
            return [row[column].strip() for row in reader if row[column].strip()]
        next(reader)  # header
        return [row[0].strip() for row in reader if row and row[0].strip()]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--metaphlan", default="results/03_metaphlan/merged_abundance_table.txt")
    parser.add_argument("--out", required=True)
    parser.add_argument("--rank", default="species", choices=sorted(PREFIX), help="taxonomic rank [species]")
    parser.add_argument("--sample-suffix", default="_metagenome$", help="regex removed from the end of sample names")
    parser.add_argument("--keep-samples-from", help="metadata CSV: keep only these samples, in this order")
    parser.add_argument("--sample-column", help="column of that CSV with the sample IDs [first column]")
    args = parser.parse_args()

    samples, rows = read_metaphlan(args.metaphlan)
    samples = [re.sub(args.sample_suffix, "", s) for s in samples]
    prefix = PREFIX[args.rank]

    table, seen = [], {}
    for clade, values in rows:
        label = clade.split("|")[-1]
        if not label.startswith(prefix):
            continue
        name = label[len(prefix):]
        if name in seen:
            seen[name] += 1
            print(f"WARNING: taxon name '{name}' occurs more than once at rank {args.rank}; renamed", file=sys.stderr)
            name = f"{name}_{seen[name]}"
        else:
            seen[name] = 1
        table.append((name, values))
    if not table:
        sys.exit(f"ERROR: no {args.rank}-level rows ('{prefix}...') found in {args.metaphlan}")

    order = list(range(len(samples)))
    if args.keep_samples_from:
        wanted = sample_order(args.keep_samples_from, args.sample_column)
        missing = [s for s in wanted if s not in samples]
        if missing:
            sys.exit(f"ERROR: {len(missing)} sample(s) of {args.keep_samples_from} are not in the MetaPhlAn table: "
                     f"{', '.join(missing[:5])}{' ...' if len(missing) > 5 else ''}\n"
                     f"       MetaPhlAn samples are named e.g. {', '.join(samples[:3])}")
        order = [samples.index(s) for s in wanted]
        print(f"Keeping {len(order)} of {len(samples)} samples, in the order of {args.keep_samples_from}", file=sys.stderr)

    with open(args.out, "w", newline="") as handle:
        writer = csv.writer(handle, quoting=csv.QUOTE_NONNUMERIC, lineterminator="\n")
        writer.writerow([""] + [samples[i] for i in order])
        for name, values in table:
            writer.writerow([name] + [float(values[i]) for i in order])
    print(f"Wrote {args.out}: {len(table)} {args.rank}-level taxa x {len(order)} samples "
          f"(MetaPhlAn relative abundance in percent, not counts)", file=sys.stderr)


if __name__ == "__main__":
    main()
