#!/usr/bin/env python3
"""Interactive residue mapper for protein multiple sequence alignments.

Dependencies:
- Python 3 standard library only

Usage:
- python3 protein_babelfish.py alignment.fasta
"""

from __future__ import annotations

import argparse
import re
import shutil
import sys
import unicodedata
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Optional

try:
    # Importing readline enables arrow-key line editing for input() on Unix-like systems.
    import readline  # noqa: F401
except ImportError:
    pass

GAP_CHARS = {"-", "."}
CONTEXT_WINDOW = 5
HEADER_TAGS = ("OS", "OX", "GN", "PE", "SV")
HEADER_FIELD_PATTERN = "|".join(HEADER_TAGS)
UNIPROT_HEADER_RE = re.compile(r"^(?P<db>[^|]+)\|(?P<accession>[^|]+)\|(?P<entry>\S+)\s+(?P<tail>.+)$")
HEADER_FIELD_RE = re.compile(
    rf"\b(?P<tag>{HEADER_FIELD_PATTERN})=(?P<value>.*?)(?=\s(?:{HEADER_FIELD_PATTERN})=|$)"
)
SEQUENCE_CHARS_RE = re.compile(r"^[A-Z*.-]+$")
NUMERIC_QUERY_RE = re.compile(r"(?:([A-Za-z]{1,20})(?:\s*-\s*|\s*)|)(\d+)")
STRIP_QUERY_CHARS = " \t\r\n\"'`“”‘’,.;:"
INPUT_TRANSLATION = str.maketrans(
    {
        "\u00A0": " ",
        "\u1680": " ",
        "\u2000": " ",
        "\u2001": " ",
        "\u2002": " ",
        "\u2003": " ",
        "\u2004": " ",
        "\u2005": " ",
        "\u2006": " ",
        "\u2007": " ",
        "\u2008": " ",
        "\u2009": " ",
        "\u200A": " ",
        "\u202F": " ",
        "\u205F": " ",
        "\u3000": " ",
        "\u200B": "",
        "\u200C": "",
        "\u200D": "",
        "\u2060": "",
        "\uFEFF": "",
        "\u2010": "-",
        "\u2011": "-",
        "\u2012": "-",
        "\u2013": "-",
        "\u2014": "-",
        "\u2015": "-",
        "\u2212": "-",
        "\uFE58": "-",
        "\uFE63": "-",
        "\uFF0D": "-",
        "\u204E": "*",
        "\u2217": "*",
        "\u2731": "*",
        "\u274B": "*",
        "\uFF0A": "*",
        "\u2018": "'",
        "\u2019": "'",
        "\u201C": '"',
        "\u201D": '"',
    }
)

AMINO_ACID_CODES = {
    "A": "A",
    "ALA": "A",
    "ALANINE": "A",
    "R": "R",
    "ARG": "R",
    "ARGININE": "R",
    "N": "N",
    "ASN": "N",
    "ASPARAGINE": "N",
    "D": "D",
    "ASP": "D",
    "ASPARTATE": "D",
    "ASPARTICACID": "D",
    "C": "C",
    "CYS": "C",
    "CYSTEINE": "C",
    "U": "U",
    "SEC": "U",
    "SELENOCYSTEINE": "U",
    "Q": "Q",
    "GLN": "Q",
    "GLUTAMINE": "Q",
    "E": "E",
    "GLU": "E",
    "GLUTAMATE": "E",
    "GLUTAMICACID": "E",
    "G": "G",
    "GLY": "G",
    "GLYCINE": "G",
    "H": "H",
    "HIS": "H",
    "HISTIDINE": "H",
    "I": "I",
    "ILE": "I",
    "ISOLEUCINE": "I",
    "J": "J",
    "XLE": "J",
    "L": "L",
    "LEU": "L",
    "LEUCINE": "L",
    "K": "K",
    "LYS": "K",
    "LYSINE": "K",
    "M": "M",
    "MET": "M",
    "METHIONINE": "M",
    "O": "O",
    "PYL": "O",
    "PYRROLYSINE": "O",
    "F": "F",
    "PHE": "F",
    "PHENYLALANINE": "F",
    "P": "P",
    "PRO": "P",
    "PROLINE": "P",
    "S": "S",
    "SER": "S",
    "SERINE": "S",
    "T": "T",
    "THR": "T",
    "THREONINE": "T",
    "W": "W",
    "TRP": "W",
    "TRYPTOPHAN": "W",
    "Y": "Y",
    "TYR": "Y",
    "TYROSINE": "Y",
    "B": "B",
    "ASX": "B",
    "Z": "Z",
    "GLX": "Z",
    "X": "X",
    "XAA": "X",
    "V": "V",
    "VAL": "V",
    "VALINE": "V",
}


@dataclass(frozen=True)
class HeaderInfo:
    database: str
    accession: str
    entry_name: str
    protein_name: str
    species: Optional[str]
    taxid: Optional[str]
    gene: Optional[str]
    protein_existence: Optional[str]
    sequence_version: Optional[str]

    @property
    def accession_label(self) -> str:
        return self.accession or "-"

    @property
    def species_label(self) -> str:
        return self.species or "Unknown species"

    @property
    def gene_label(self) -> str:
        return self.gene or "-"


@dataclass(frozen=True)
class AlignmentRecord:
    header: HeaderInfo
    raw_header: str
    aligned_sequence: str
    ungapped_sequence: str
    residue_to_alignment: tuple[int, ...]
    alignment_to_residue: tuple[Optional[int], ...]

    @property
    def search_text(self) -> str:
        return " ".join(
            [
                self.header.accession,
                self.header.entry_name,
                self.header.gene or "",
                self.header.species or "",
                self.header.protein_name,
                self.raw_header,
            ]
        ).casefold()


@dataclass(frozen=True)
class ResidueQuery:
    position: int
    residue: str
    query_text: str

    @property
    def residue_label(self) -> str:
        return f"{self.residue}{self.position}"


class QueryError(Exception):
    pass


def normalize_text(text: str, *, collapse_whitespace: bool = True) -> str:
    text = unicodedata.normalize("NFKC", text).translate(INPUT_TRANSLATION)
    if collapse_whitespace:
        text = re.sub(r"\s+", " ", text)
    return text.strip()


def prompt(message: str) -> str:
    try:
        return normalize_text(input(message))
    except EOFError as exc:
        raise SystemExit("\nInput ended unexpectedly.") from exc


def clean_query_text(text: str) -> str:
    return normalize_text(text).strip(STRIP_QUERY_CHARS)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Interactively map a residue from one sequence in a FASTA multiple "
            "sequence alignment to every other sequence in the alignment."
        ),
        add_help=False,
    )
    parser.add_argument("-h", "--help", "--h", action="help", help="show this help message and exit")
    parser.add_argument("alignment_fasta", nargs="?", type=Path, help="Path to the aligned FASTA file")
    args = parser.parse_args()
    if args.alignment_fasta is None:
        parser.print_help()
        raise SystemExit(0)
    return args


def parse_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    current_header: Optional[str] = None
    current_sequence: list[str] = []

    try:
        lines = path.read_text(encoding="utf-8-sig").splitlines()
    except UnicodeDecodeError as exc:
        raise QueryError(f"Could not decode alignment file as UTF-8 text: {path}") from exc
    except OSError as exc:
        raise QueryError(f"Could not read alignment file: {path}") from exc

    for line_number, line in enumerate(lines, start=1):
        line = normalize_text(line, collapse_whitespace=False)
        if not line:
            continue
        if line.startswith(">"):
            if current_header is not None:
                records.append((current_header, "".join(current_sequence)))
            current_header = normalize_text(line[1:])
            current_sequence = []
            continue
        if current_header is None:
            raise QueryError(f"Found sequence data before the first FASTA header on line {line_number}.")

        chunk = re.sub(r"\s+", "", line).upper()
        if not SEQUENCE_CHARS_RE.fullmatch(chunk):
            raise QueryError(
                "Unexpected character(s) in sequence data on "
                f"line {line_number}. Allowed characters are letters, '*', '.', and '-'."
            )
        current_sequence.append(chunk)

    if current_header is not None:
        records.append((current_header, "".join(current_sequence)))

    if not records:
        raise QueryError("The alignment file did not contain any FASTA records.")

    lengths = {len(sequence) for _, sequence in records}
    if len(lengths) != 1:
        raise QueryError(
            "All sequences in the input FASTA must already be aligned to the same length."
        )

    return records


def parse_uniprot_header(header: str) -> HeaderInfo:
    prefix_match = UNIPROT_HEADER_RE.match(header)
    if not prefix_match:
        first_token = header.split()[0] if header.split() else "sequence"
        return HeaderInfo(
            database="unknown",
            accession="",
            entry_name=first_token,
            protein_name=header,
            species=None,
            taxid=None,
            gene=None,
            protein_existence=None,
            sequence_version=None,
        )

    tail = prefix_match.group("tail")
    tagged_fields = {
        match.group("tag"): normalize_text(match.group("value"))
        for match in HEADER_FIELD_RE.finditer(tail)
    }
    description = normalize_text(HEADER_FIELD_RE.sub("", tail)) or prefix_match.group("entry")

    return HeaderInfo(
        database=prefix_match.group("db"),
        accession=prefix_match.group("accession"),
        entry_name=prefix_match.group("entry"),
        protein_name=description,
        species=tagged_fields.get("OS"),
        taxid=tagged_fields.get("OX"),
        gene=tagged_fields.get("GN"),
        protein_existence=tagged_fields.get("PE"),
        sequence_version=tagged_fields.get("SV"),
    )


def build_alignment_record(raw_header: str, sequence: str) -> AlignmentRecord:
    residue_to_alignment: list[int] = []
    alignment_to_residue: list[Optional[int]] = []
    ungapped_sequence: list[str] = []
    residue_number = 0

    for alignment_index, character in enumerate(sequence):
        if character in GAP_CHARS:
            alignment_to_residue.append(None)
            continue
        residue_number += 1
        residue_to_alignment.append(alignment_index)
        alignment_to_residue.append(residue_number)
        ungapped_sequence.append(character)

    return AlignmentRecord(
        header=parse_uniprot_header(raw_header),
        raw_header=raw_header,
        aligned_sequence=sequence,
        ungapped_sequence="".join(ungapped_sequence),
        residue_to_alignment=tuple(residue_to_alignment),
        alignment_to_residue=tuple(alignment_to_residue),
    )


def load_alignment(path: Path) -> list[AlignmentRecord]:
    return [build_alignment_record(header, sequence) for header, sequence in parse_fasta(path)]


def prompt_metadata_value(label: str, current_value: Optional[str]) -> str:
    while True:
        suffix = f" [{current_value}]" if current_value else ""
        response = prompt(f"{label}{suffix}: ")
        if response:
            return response
        if current_value:
            return current_value
        print(f"{label} cannot be empty.")


def complete_sequence_metadata(records: list[AlignmentRecord]) -> list[AlignmentRecord]:
    needs_review = any(not record.header.species or not record.header.gene for record in records)
    if not needs_review:
        return records

    print(
        "\nCould not extract complete species and gene metadata from all FASTA headers.\n"
        "Please review the sequences below. Press Enter to keep any value shown in brackets."
    )

    updated_records: list[AlignmentRecord] = []
    for index, record in enumerate(records, start=1):
        print(f"\nSequence {index}")
        print(f"Header: {record.raw_header}")
        species = prompt_metadata_value("Species", record.header.species)
        gene = prompt_metadata_value("Gene name", record.header.gene)
        updated_records.append(
            replace(
                record,
                header=replace(record.header, species=species, gene=gene),
            )
        )

    return updated_records


def normalize_residue_token(token: str) -> Optional[str]:
    cleaned = re.sub(r"[^A-Za-z]", "", normalize_text(token)).upper()
    if not cleaned:
        return None
    return AMINO_ACID_CODES.get(cleaned)


def unique_output_path(base_path: Path) -> Path:
    if not base_path.exists():
        return base_path

    stem = base_path.stem
    suffix = base_path.suffix
    counter = 2
    while True:
        candidate = base_path.with_name(f"{stem}_{counter}{suffix}")
        if not candidate.exists():
            return candidate
        counter += 1


def sanitize_filename(fragment: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", normalize_text(fragment)).strip("_") or "residue_mapping"


def format_context(sequence: str, index: int, *, one_based: bool = False, window: int = CONTEXT_WINDOW) -> str:
    zero_based = index - 1 if one_based else index
    left = sequence[max(0, zero_based - window) : zero_based]
    center = sequence[zero_based]
    right = sequence[zero_based + 1 : zero_based + 1 + window]
    return f"{left}*{center}*{right}"


def clip_text(text: str, width: int) -> str:
    if len(text) <= width:
        return text
    if width <= 3:
        return text[:width]
    return f"{text[: width - 3]}..."


def find_all_matches(sequence: str, motif: str) -> list[int]:
    positions: list[int] = []
    start = 0
    while True:
        match_index = sequence.find(motif, start)
        if match_index == -1:
            return positions
        positions.append(match_index)
        start = match_index + 1


def resolve_numeric_query(raw_query: str, record: AlignmentRecord) -> ResidueQuery:
    match = NUMERIC_QUERY_RE.fullmatch(raw_query)
    if not match:
        raise QueryError(
            "Enter a residue as 235, R235, Arg235, Arg-235, R 235, or a motif such as RFTH*R*ADTV."
        )

    residue_token, position_text = match.groups()
    position = int(position_text)
    sequence_length = len(record.ungapped_sequence)

    if position < 1 or position > sequence_length:
        raise QueryError(
            f"Position {position} is outside the selected sequence length (1-{sequence_length})."
        )

    actual_residue = record.ungapped_sequence[position - 1]
    if residue_token:
        expected_residue = normalize_residue_token(residue_token)
        if expected_residue is None:
            raise QueryError(
                f"'{residue_token}' is not a recognized amino-acid code. Use one-letter or three-letter residue names."
            )
        if expected_residue != actual_residue:
            context = format_context(record.ungapped_sequence, position, one_based=True)
            raise QueryError(
                f"There is no {expected_residue} at position {position} in the selected sequence: {context}"
            )

    return ResidueQuery(
        position=position,
        residue=actual_residue,
        query_text=raw_query.strip(),
    )


def resolve_motif_query(raw_query: str, record: AlignmentRecord) -> ResidueQuery:
    if raw_query.count("*") != 2:
        raise QueryError(
            "Motif queries must mark the residue of interest with asterisks, for example RFTH*R*ADTV."
        )

    left_raw, center_raw, right_raw = raw_query.split("*")
    left = re.sub(r"[^A-Za-z]", "", normalize_text(left_raw)).upper()
    center_token = re.sub(r"[^A-Za-z]", "", normalize_text(center_raw))
    right = re.sub(r"[^A-Za-z]", "", normalize_text(right_raw)).upper()

    if not center_token:
        raise QueryError("The marked residue inside the motif is empty. Use a format such as RFTH*R*ADTV.")

    focus_residue = normalize_residue_token(center_token)
    if focus_residue is None:
        raise QueryError(
            f"'{center_raw.strip()}' is not a recognized amino-acid code inside the motif."
        )

    motif = f"{left}{focus_residue}{right}"
    if not motif:
        raise QueryError("The motif is empty after normalization.")

    matches = find_all_matches(record.ungapped_sequence.upper(), motif)
    if not matches:
        raise QueryError(f"The motif {left}*{focus_residue}*{right} was not found in the selected sequence.")
    if len(matches) > 1:
        positions = ", ".join(str(match + len(left) + 1) for match in matches[:10])
        more = "..." if len(matches) > 10 else ""
        raise QueryError(
            f"The motif {left}*{focus_residue}*{right} matched multiple positions ({positions}{more}). "
            "Please provide a longer unique motif."
        )

    position = matches[0] + len(left) + 1
    return ResidueQuery(
        position=position,
        residue=focus_residue,
        query_text=raw_query.strip(),
    )


def resolve_residue_query(raw_query: str, record: AlignmentRecord) -> ResidueQuery:
    if "*" in raw_query:
        return resolve_motif_query(raw_query, record)
    return resolve_numeric_query(raw_query, record)


def print_sequence_choices(records: list[AlignmentRecord]) -> None:
    terminal_width = shutil.get_terminal_size(fallback=(120, 20)).columns
    rows = [
        {
            "No.": str(index),
            "Entry": record.header.entry_name or record.header.accession_label,
            "Accession": record.header.accession_label,
            "Species": record.header.species_label,
            "Gene": record.header.gene_label,
            "TaxID": record.header.taxid or "-",
            "Protein": record.header.protein_name,
        }
        for index, record in enumerate(records, start=1)
    ]
    fixed_columns = ("No.", "Entry", "Accession", "Species", "Gene", "TaxID")
    fixed_width = sum(
        max(len(column), *(len(row[column]) for row in rows))
        for column in fixed_columns
    )
    separator_width = 3 * (len(rows[0]) - 1)
    protein_width = max(12, terminal_width - fixed_width - separator_width)
    for row in rows:
        row["Protein"] = clip_text(row["Protein"], protein_width)

    print("\nAvailable sequences:")
    print(render_plain_table(rows))
    print("Enter the list number, accession, entry name, gene, or a unique text fragment from the list above.")


def resolve_sequence_choice(response: str, records: list[AlignmentRecord]) -> AlignmentRecord:
    choice = normalize_text(response).strip(STRIP_QUERY_CHARS)
    if not choice:
        raise QueryError("Please choose one of the sequences listed above.")

    if choice.isdigit():
        index = int(choice)
        if 1 <= index <= len(records):
            return records[index - 1]
        raise QueryError(f"Sequence {index} is not in the list. Choose a number between 1 and {len(records)}.")

    lowered_choice = choice.casefold()
    matches = [record for record in records if lowered_choice in record.search_text]

    if not matches:
        raise QueryError(f"No sequence matched '{choice}'.")
    if len(matches) > 1:
        options = ", ".join(record.header.entry_name for record in matches)
        raise QueryError(f"'{choice}' matched more than one sequence: {options}. Please be more specific.")
    return matches[0]


def prompt_for_sequence(records: list[AlignmentRecord]) -> AlignmentRecord:
    while True:
        print_sequence_choices(records)
        response = prompt("\nWhich sequence would you like to convert from? ")
        try:
            return resolve_sequence_choice(response, records)
        except QueryError as exc:
            print(f"\n{exc}")


def prompt_after_query_error() -> str:
    while True:
        response = prompt(
            "Press Enter to try again, type 's' to choose another sequence, or 'q' to quit: "
        ).lower()
        if response in {"", "s", "q"}:
            return response
        print("Please press Enter, or type 's' or 'q'.")


def prompt_for_query(record: AlignmentRecord) -> Optional[ResidueQuery]:
    print(
        f"\nSelected sequence: {record.header.entry_name} | {record.header.accession_label} | "
        f"{record.header.species_label} | {record.header.protein_name}"
    )

    while True:
        raw_query = clean_query_text(prompt("Which residue would you like to convert numbering for? "))
        if not raw_query:
            print("Please enter a residue number or a marked motif, for example 235 or RFTH*R*ADTV.")
            continue

        try:
            return resolve_residue_query(raw_query, record)
        except QueryError as exc:
            print(f"\n{exc}")
            follow_up = prompt_after_query_error()
            if follow_up == "s":
                return None
            if follow_up == "q":
                raise SystemExit(0)


def build_results(records: list[AlignmentRecord], source_record: AlignmentRecord, query: ResidueQuery) -> list[dict[str, str]]:
    source_alignment_index = source_record.residue_to_alignment[query.position - 1]
    results: list[dict[str, str]] = []

    for record in records:
        aligned_character = record.aligned_sequence[source_alignment_index]
        context = format_context(record.aligned_sequence, source_alignment_index)

        if record is source_record:
            equivalent = query.residue_label
            status = "source"
        elif aligned_character in GAP_CHARS:
            equivalent = "nonconserved"
            status = "nonconserved"
        else:
            residue_number = record.alignment_to_residue[source_alignment_index]
            equivalent = f"{aligned_character}{residue_number}"
            status = "conserved" if aligned_character == query.residue else "substituted"

        results.append(
            {
                "Sequence": record.header.entry_name,
                "Species": record.header.species_label,
                "Gene": record.header.gene_label,
                "Accession": record.header.accession_label,
                "Equivalent": equivalent,
                "Status": status,
                "Context": context,
            }
        )

    return results


def render_plain_table(rows: list[dict[str, str]]) -> str:
    columns = list(rows[0].keys())
    widths = {
        column: max(len(column), *(len(row[column]) for row in rows))
        for column in columns
    }

    header = " | ".join(column.ljust(widths[column]) for column in columns)
    separator = "-+-".join("-" * widths[column] for column in columns)
    body = [
        " | ".join(row[column].ljust(widths[column]) for column in columns)
        for row in rows
    ]
    return "\n".join([header, separator, *body])


def escape_markdown(value: str) -> str:
    return value.replace("|", r"\|")


def write_markdown_report(
    rows: list[dict[str, str]],
    alignment_path: Path,
    source_record: AlignmentRecord,
    query: ResidueQuery,
) -> Path:
    filename = sanitize_filename(
        f"residue_mapping_{source_record.header.entry_name}_{query.residue_label}"
    )
    output_path = unique_output_path(Path.cwd() / f"{filename}.md")

    lines = [
        "# Residue Mapping Report",
        "",
        f"- Alignment: `{alignment_path}`",
        (
            f"- Source sequence: `{source_record.header.entry_name}` "
            f"({source_record.header.accession_label}; {source_record.header.species_label})"
        ),
        f"- Source residue: `{query.residue_label}`",
        f"- Query input: `{query.query_text}`",
        "",
        "| " + " | ".join(rows[0].keys()) + " |",
        "| " + " | ".join("---" for _ in rows[0]) + " |",
    ]

    for row in rows:
        lines.append("| " + " | ".join(escape_markdown(value) for value in row.values()) + " |")

    output_path.write_text("\n".join(lines) + "\n")
    return output_path


def prompt_for_restart() -> bool:
    response = prompt("\nWould you like to make another query? [y/N] ").lower()
    return response in {"y", "yes"}


def run_interactive_session(records: list[AlignmentRecord], alignment_path: Path) -> int:
    print(f"Loaded {len(records)} aligned sequences from {alignment_path}")

    while True:
        source_record = prompt_for_sequence(records)
        query = prompt_for_query(source_record)
        if query is None:
            continue

        results = build_results(records, source_record, query)
        markdown_path = write_markdown_report(results, alignment_path, source_record, query)

        print(
            f"\nResolved {query.query_text!r} to {query.residue_label} in "
            f"{source_record.header.entry_name} ({source_record.header.species_label}).\n"
        )
        print(render_plain_table(results))
        print(f"\nMarkdown table written to {markdown_path}")

        if not prompt_for_restart():
            return 0


def main() -> int:
    args = parse_arguments()
    try:
        records = complete_sequence_metadata(load_alignment(args.alignment_fasta))
        return run_interactive_session(records, args.alignment_fasta)
    except QueryError as exc:
        print(exc, file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("\nInterrupted.", file=sys.stderr)
        return 130


if __name__ == "__main__":
    raise SystemExit(main())
