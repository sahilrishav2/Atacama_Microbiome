#!/usr/bin/env python3
"""Keeps the samples listed in a metadata CSV, in the order of that file, from a taxa x samples CSV table.

The downstream documents require the columns of every table to be in the same order as the samples of the metadata.
Input and output: CSV with the features in the first column and one column per sample. Values are passed on unchanged.
  --keep-samples-from METADATA.csv   [--sample-column NAME]   (default: the first column holds the sample IDs)
Only the Python standard library is used.
"""
import argparse
import csv
import sys


def number(text):
    for kind in (int, float):
        try:
            return kind(text)
        except ValueError:
            pass
    return text


def sample_ids(path, column):
    with open(path, newline="") as handle:
        if column:
            reader = csv.DictReader(handle)
            if column not in (reader.fieldnames or []):
                sys.exit(f"ERROR: column '{column}' not found in {path} (columns: {', '.join(reader.fieldnames or [])})")
            return [row[column].strip() for row in reader if row[column].strip()]
        reader = csv.reader(handle)
        next(reader, None)
        return [row[0].strip() for row in reader if row and row[0].strip()]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--table", required=True)
    parser.add_argument("--keep-samples-from", required=True)
    parser.add_argument("--sample-column")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    with open(args.table, newline="") as handle:
        rows = list(csv.reader(handle))
    if len(rows) < 2:
        sys.exit(f"ERROR: {args.table} has no data rows")
    header = rows[0][1:]
    wanted = sample_ids(args.keep_samples_from, args.sample_column)
    missing = [s for s in wanted if s not in header]
    if missing:
        sys.exit(f"ERROR: {len(missing)} sample(s) of {args.keep_samples_from} are not columns of {args.table}: "
                 f"{', '.join(missing[:5])}{' ...' if len(missing) > 5 else ''}\n"
                 f"       the table's samples are named e.g. {', '.join(header[:3])}")
    order = [header.index(s) + 1 for s in wanted]
    with open(args.out, "w", newline="") as handle:
        writer = csv.writer(handle, quoting=csv.QUOTE_NONNUMERIC, lineterminator="\n")
        writer.writerow([rows[0][0]] + wanted)
        for row in rows[1:]:
            writer.writerow([row[0]] + [number(row[i]) for i in order])
    print(f"Wrote {args.out}: {len(rows) - 1} rows x {len(wanted)} of {len(header)} samples "
          f"(order of {args.keep_samples_from})", file=sys.stderr)


if __name__ == "__main__":
    main()
