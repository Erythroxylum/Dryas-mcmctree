#!/usr/bin/env python3
"""Filter p28 PHYLIP alignments by occupancy and site-level error indicators.

A site is retained when all of the following are true:
  1. at least --min-called-taxa taxa have a nonmissing call;
  2. it is not an excess-ambiguity site, defined as at least
     --min-heterozygous-taxa IUPAC calls AND an IUPAC fraction of at least
     --heterozygous-fraction-threshold among called taxa;
  3. no more than --max-unambiguous-alleles distinct A/C/G/T bases occur.

Invariant sites are retained whenever they pass the occupancy criterion.
IUPAC ambiguity codes count as called for occupancy and as heterozygous or
ambiguous for the excess-heterozygosity diagnostic.
"""

import argparse
import csv
import gzip
import re
from collections import Counter
from pathlib import Path

ACGT = set(b"ACGT")
MISSING = set(b"N-?")
IUPAC = set(b"RYSWKMBDHV")
ALLOWED = ACGT | MISSING | IUPAC


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alignment-dir", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--min-called-taxa", type=int, default=16)
    parser.add_argument("--min-heterozygous-taxa", type=int, default=4)
    parser.add_argument("--heterozygous-fraction-threshold", type=float, default=0.25)
    parser.add_argument("--max-unambiguous-alleles", type=int, default=2)
    return parser.parse_args()


def chromosome_number(path):
    hit = re.search(r"chr0*([1-9])", path.name, flags=re.IGNORECASE)
    return int(hit.group(1)) if hit else None


def find_alignments(directory):
    candidates = sorted(directory.glob("*.phy"))
    by_chr = {}
    for path in candidates:
        chromosome = chromosome_number(path)
        if chromosome is not None:
            by_chr.setdefault(chromosome, []).append(path)
    files = []
    for chromosome in range(1, 10):
        matches = by_chr.get(chromosome, [])
        if len(matches) != 1:
            raise ValueError(
                f"Expected one chr{chromosome:02d} PHYLIP in {directory}; "
                f"found {len(matches)}"
            )
        files.append(matches[0])
    return files


def read_phylip(path):
    with path.open("rb") as handle:
        header = handle.readline().split()
        if len(header) < 2:
            raise ValueError(f"Invalid PHYLIP header: {path}")
        ntaxa, nsites = map(int, header[:2])
        names, sequences = [], []
        for line_number, raw in enumerate(handle, 2):
            if not raw.strip():
                continue
            fields = raw.split()
            if len(fields) != 2:
                raise ValueError(f"Expected sample and sequence at {path}:{line_number}")
            name = fields[0].decode("utf-8")
            sequence = fields[1].upper().replace(b"U", b"T")
            if len(sequence) != nsites:
                raise ValueError(
                    f"{name} has {len(sequence)} sites in {path}; expected {nsites}"
                )
            unexpected = set(sequence) - ALLOWED
            if unexpected:
                raise ValueError(
                    f"Unexpected characters for {name} in {path}: "
                    f"{sorted(chr(x) for x in unexpected)}"
                )
            names.append(name)
            sequences.append(sequence)
    if len(names) != ntaxa or len(set(names)) != ntaxa:
        raise ValueError(f"{path} declares {ntaxa} taxa; read {len(names)} unique names")
    return names, sequences, nsites


def write_phylip(path, names, sequences):
    nsites = len(sequences[0]) if sequences else 0
    if not sequences or any(len(sequence) != nsites for sequence in sequences):
        raise ValueError("Cannot write an empty or unequal-length alignment")
    with path.open("wb") as handle:
        handle.write(f"{len(names)} {nsites}\n".encode())
        for name, sequence in zip(names, sequences):
            handle.write(name.encode("utf-8") + b"  " + bytes(sequence) + b"\n")


def reason_key(reasons):
    order = ("low_occupancy", "excess_heterozygosity", "multiallelic")
    return ";".join(reason for reason in order if reason in reasons)


def flush_interval(writer, chromosome, start, end, reason):
    if start is not None:
        writer.writerow([chromosome, start, end, end - start + 1, reason])


def filter_alignment(path, outdir, args, interval_writer):
    names, sequences, nsites = read_phylip(path)
    ntaxa = len(names)
    if args.min_called_taxa > ntaxa:
        raise ValueError(f"Minimum called taxa exceeds {ntaxa} in {path}")

    chromosome = f"chr{chromosome_number(path):02d}"
    retained = [bytearray() for _ in sequences]
    counts = Counter(input_sites=nsites, taxa=ntaxa)
    current_reason = None
    interval_start = None

    for index in range(nsites):
        column = [sequence[index] for sequence in sequences]
        called = sum(base not in MISSING for base in column)
        heterozygous = sum(base in IUPAC for base in column)
        heterozygous_fraction = heterozygous / called if called else 0.0
        unambiguous_alleles = len({base for base in column if base in ACGT})

        low_occupancy = called < args.min_called_taxa
        excess_heterozygosity = (
            heterozygous >= args.min_heterozygous_taxa
            and heterozygous_fraction >= args.heterozygous_fraction_threshold
        )
        multiallelic = unambiguous_alleles > args.max_unambiguous_alleles
        reasons = []
        if low_occupancy:
            reasons.append("low_occupancy")
            counts["failed_low_occupancy"] += 1
        if excess_heterozygosity:
            reasons.append("excess_heterozygosity")
            counts["failed_excess_heterozygosity"] += 1
        if multiallelic:
            reasons.append("multiallelic")
            counts["failed_multiallelic"] += 1

        if reasons:
            counts["removed_sites"] += 1
            key = reason_key(reasons)
            counts[f"reason__{key}"] += 1
            if key != current_reason:
                flush_interval(
                    interval_writer, chromosome, interval_start, index, current_reason
                )
                interval_start = index + 1
                current_reason = key
        else:
            flush_interval(
                interval_writer, chromosome, interval_start, index, current_reason
            )
            interval_start = None
            current_reason = None
            counts["retained_sites"] += 1
            if unambiguous_alleles <= 1:
                counts["retained_invariant_sites"] += 1
            else:
                counts["retained_variable_sites"] += 1
            for output, base in zip(retained, column):
                output.append(base)

    flush_interval(interval_writer, chromosome, interval_start, nsites, current_reason)
    if counts["retained_sites"] == 0:
        raise ValueError(f"All sites were removed from {path}")

    stem = re.sub(r"\.trimmed$", "", path.stem, flags=re.IGNORECASE)
    output_path = outdir / f"{stem}.sitefiltered.phy"
    write_phylip(output_path, names, retained)

    # Re-read the output as an independent format and length validation.
    out_names, out_sequences, out_nsites = read_phylip(output_path)
    if out_names != names or out_nsites != counts["retained_sites"]:
        raise RuntimeError(f"Output validation failed: {output_path}")
    if any(len(sequence) != out_nsites for sequence in out_sequences):
        raise RuntimeError(f"Unequal output sequence lengths: {output_path}")

    return {
        "chromosome": chromosome,
        "input_alignment": str(path.resolve()),
        "output_alignment": str(output_path.resolve()),
        "taxa": ntaxa,
        "input_sites": nsites,
        "retained_sites": counts["retained_sites"],
        "removed_sites": counts["removed_sites"],
        "retained_fraction": counts["retained_sites"] / nsites,
        "retained_invariant_sites": counts["retained_invariant_sites"],
        "retained_variable_sites": counts["retained_variable_sites"],
        "failed_low_occupancy": counts["failed_low_occupancy"],
        "failed_excess_heterozygosity": counts["failed_excess_heterozygosity"],
        "failed_multiallelic": counts["failed_multiallelic"],
        "removed_low_occupancy_only": counts["reason__low_occupancy"],
        "removed_excess_heterozygosity_only": counts[
            "reason__excess_heterozygosity"
        ],
        "removed_multiallelic_only": counts["reason__multiallelic"],
        "removed_multiple_reasons": sum(
            value
            for key, value in counts.items()
            if key.startswith("reason__") and ";" in key
        ),
    }


def main():
    args = parse_args()
    if args.min_called_taxa < 1 or args.min_heterozygous_taxa < 1:
        raise ValueError("Taxon-count thresholds must be positive")
    if not 0 <= args.heterozygous_fraction_threshold <= 1:
        raise ValueError("--heterozygous-fraction-threshold must be between 0 and 1")
    if args.max_unambiguous_alleles < 1:
        raise ValueError("--max-unambiguous-alleles must be positive")

    args.outdir.mkdir(parents=True, exist_ok=True)
    files = find_alignments(args.alignment_dir)
    interval_path = args.outdir / "removed_site_intervals.tsv.gz"
    rows = []
    with gzip.open(interval_path, "wt", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ["chromosome", "start_site_1based", "end_site_1based", "sites", "reason"]
        )
        for path in files:
            row = filter_alignment(path, args.outdir, args, writer)
            rows.append(row)
            print(
                f"{row['chromosome']}: retained {row['retained_sites']:,} / "
                f"{row['input_sites']:,} sites ({100 * row['retained_fraction']:.1f}%)"
            )

    summary_path = args.outdir / "site_filter_summary.tsv"
    with summary_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)

    settings_path = args.outdir / "site_filter_settings.tsv"
    settings = [
        ("minimum_called_taxa", args.min_called_taxa),
        ("minimum_heterozygous_taxa", args.min_heterozygous_taxa),
        ("minimum_heterozygous_fraction", args.heterozygous_fraction_threshold),
        ("maximum_unambiguous_ACGT_alleles", args.max_unambiguous_alleles),
        ("missing_characters", "N-?"),
        ("heterozygous_or_ambiguous_IUPAC_characters", "RYSWKMBDHV"),
    ]
    with settings_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["setting", "value"])
        writer.writerows(settings)

    print(f"Wrote filtered alignments and reports to: {args.outdir.resolve()}")


if __name__ == "__main__":
    main()
