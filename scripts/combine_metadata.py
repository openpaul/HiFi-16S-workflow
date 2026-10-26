#!/usr/bin/env python3
"""Collapse a QIIME 2 metadata file to one row per group.

This is used when several libraries (for example, the same DNA extract
sequenced in two different runs) should be analysed as one biological sample.
The feature table is grouped with ``qiime feature-table group``; this script
creates the matching sample metadata.

Columns that are constant within a group are kept. Columns that vary between
libraries of the same group (sequencing run, library prep date, ...) are
dropped, because they no longer describe the combined sample, and reported on
stderr.
"""

import argparse
import csv
import sys
from collections import OrderedDict


def read_metadata(path):
    with open(path, newline="") as handle:
        rows = [row for row in csv.reader(handle, delimiter="\t") if row]
    if not rows:
        raise SystemExit("error: metadata file %s is empty" % path)

    header = rows[0]
    if len(rows) > 1 and rows[1] and rows[1][0] == "#q2:types":
        types = rows[1]
        data = rows[2:]
    else:
        types = None
        data = rows[1:]
    return header, types, data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="input QIIME 2 metadata TSV")
    parser.add_argument("--output", required=True, help="output combined metadata TSV")
    parser.add_argument(
        "--column", required=True, help="metadata column holding the group (e.g. subject)"
    )
    parser.add_argument(
        "--id-column",
        default=None,
        help="sample ID column (default: first column of the input)",
    )
    args = parser.parse_args()

    header, types, data = read_metadata(args.input)
    id_column = args.id_column or header[0]

    if id_column not in header:
        raise SystemExit("error: ID column %r not found in %s" % (id_column, args.input))
    if args.column not in header:
        raise SystemExit(
            "error: group column %r not found in %s" % (args.column, args.input)
        )

    id_index = header.index(id_column)
    group_index = header.index(args.column)
    type_by_column = dict(zip(header, types)) if types else {}

    groups = OrderedDict()
    for row in data:
        if len(row) < len(header):
            row = row + [""] * (len(header) - len(row))
        group = row[group_index].strip()
        if not group:
            raise SystemExit(
                "error: empty group value in column %r for sample %r"
                % (args.column, row[id_index])
            )
        groups.setdefault(group, []).append(row)

    keep = []
    dropped = []
    for index, column in enumerate(header):
        if index in (id_index, group_index):
            continue
        for rows in groups.values():
            values = {row[index].strip() for row in rows if row[index].strip()}
            if len(values) > 1:
                dropped.append(column)
                break
        else:
            keep.append(index)

    output = [[id_column] + [header[index] for index in keep]]
    if types:
        output.append(
            ["#q2:types", "categorical"]
            + [type_by_column.get(header[index], "categorical") for index in keep]
        )

    for group, rows in groups.items():
        values = []
        for index in keep:
            seen = {row[index] for row in rows if row[index].strip()}
            values.append(next(iter(seen)) if seen else "")
        output.append([group] + values)

    with open(args.output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerows(output)

    if dropped:
        print(
            "Dropped metadata columns that differ between libraries of the same "
            "group: " + ", ".join(dropped),
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
