# Protein Babelfish

Small interactive CLI for mapping an aligned residue from one protein sequence in a FASTA alignment to the equivalent aligned position in every other sequence.

## Dependencies

- Python 3.10 or newer
- No third-party Python packages

## Run

```bash
python3 protein_babelfish.py /path/to/alignment.fasta
```

### Example Run
```bash
python protein_babelfish.py --i ryr_alignment.fasta

Loaded 13 aligned sequences from ryr_alignment.fasta

Available sequences:
No. | Entry      | Accession | Species                 | Gene | TaxID | Protein
----+------------+-----------+-------------------------+------+-------+---------------------
1   | RYR_DROME  | Q24498    | Drosophila melanogaster | RyR  | 7227  | Ryanodine receptor
2   | RYR1_RABIT | P11716    | Oryctolagus cuniculus   | RYR1 | 9986  | Ryanodine receptor 1
3   | RYR1_PIG   | P16960    | Sus scrofa              | RYR1 | 9823  | Ryanodine receptor 1
4   | RYR1_HUMAN | P21817    | Homo sapiens            | RYR1 | 9606  | Ryanodine receptor 1
5   | RYR1_MOUSE | E9PZQ0    | Mus musculus            | Ryr1 | 10090 | Ryanodine receptor 1
6   | RYR1_RAT   | F1LMY4    | Rattus norvegicus       | Ryr1 | 10116 | Ryanodine receptor 1
7   | RYR2_MOUSE | E9Q401    | Mus musculus            | Ryr2 | 10090 | Ryanodine receptor 2
8   | RYR2_RAT   | B0LPN4    | Rattus norvegicus       | Ryr2 | 10116 | Ryanodine receptor 2
9   | RYR2_RABIT | P30957    | Oryctolagus cuniculus   | RYR2 | 9986  | Ryanodine receptor 2
10  | RYR2_HUMAN | Q92736    | Homo sapiens            | RYR2 | 9606  | Ryanodine receptor 2
11  | RYR3_MOUSE | A2AGL3    | Mus musculus            | Ryr3 | 10090 | Ryanodine receptor 3
12  | RYR3_RABIT | Q9TS33    | Oryctolagus cuniculus   | RYR3 | 9986  | Ryanodine receptor 3
13  | RYR3_HUMAN | Q15413    | Homo sapiens            | RYR3 | 9606  | Ryanodine receptor 3
Enter the list number, accession, entry name, gene, or a unique text fragment from the list above.

Which sequence would you like to convert from? 4

Selected sequence: RYR1_HUMAN | P21817 | Homo sapiens | Ryanodine receptor 1
Which residue would you like to convert numbering for? R163

Resolved 'R163' to R163 in RYR1_HUMAN (Homo sapiens).

Sequence   | Species                 | Gene | Accession | Equivalent | Status    | Context
-----------+-------------------------+------+-----------+------------+-----------+--------------
RYR_DROME  | Drosophila melanogaster | RyR  | Q24498    | R159       | conserved | EGEKV*R*VGDDL
RYR1_RABIT | Oryctolagus cuniculus   | RYR1 | P11716    | R164       | conserved | EGEKV*R*VGDDL
RYR1_PIG   | Sus scrofa              | RYR1 | P16960    | R164       | conserved | EGEKV*R*VGDDL
RYR1_HUMAN | Homo sapiens            | RYR1 | P21817    | R163       | source    | EGEKV*R*VGDDI
RYR1_MOUSE | Mus musculus            | Ryr1 | E9PZQ0    | R165       | conserved | EGEKV*R*VGDDL
RYR1_RAT   | Rattus norvegicus       | Ryr1 | F1LMY4    | R165       | conserved | EGEKV*R*VGDDL
RYR2_MOUSE | Mus musculus            | Ryr2 | E9Q401    | R176       | conserved | EGEKV*R*VGDDL
RYR2_RAT   | Rattus norvegicus       | Ryr2 | B0LPN4    | R169       | conserved | EGEKV*R*VGDDL
RYR2_RABIT | Oryctolagus cuniculus   | RYR2 | P30957    | R176       | conserved | EGEKV*R*VGDDL
RYR2_HUMAN | Homo sapiens            | RYR2 | Q92736    | R176       | conserved | EGEKV*R*VGDDL
RYR3_MOUSE | Mus musculus            | Ryr3 | A2AGL3    | R166       | conserved | EGEKV*R*IGDDL
RYR3_RABIT | Oryctolagus cuniculus   | RYR3 | Q9TS33    | R166       | conserved | EGEKV*R*IGDDL
RYR3_HUMAN | Homo sapiens            | RYR3 | Q15413    | R166       | conserved | EGEKV*R*IGDDL

Markdown table written to /Users/oc2188/Downloads/residue_mapping_RYR1_HUMAN_R163_2.md

Would you like to make another query? [y/N]
```

## Input notes

- The FASTA must already be a multiple sequence alignment with all sequences the same aligned length.
- Residue queries can be entered as `235`, `R235`, `Arg-235`, `R 235`, or as a marked motif such as `RFTH*R*ADTV`.
- Common pasted formatting artifacts such as smart quotes, non-breaking spaces, Unicode dashes, and Unicode asterisks are stripped/fixed automatically.
- If headers are not UniProt-style, the script will ask for species and gene name for each sequence before starting the lookup.

## Caveat
Written with the assistence of Codex/ChatGPT5.4; Tested against a few different alignments, but there may still be bugs; let me know if you encounter any! I have provided a [sample alignment](ryr_alignment.fasta) for testing as part of the repo.
