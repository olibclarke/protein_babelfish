# Protein Babelfish

Small interactive CLI for mapping an aligned residue from one protein sequence in a FASTA alignment to the equivalent aligned position in every other sequence.

## Dependencies

- Python 3.10 or newer
- No third-party Python packages

## Run

```bash
python3 protein_babelfish.py /path/to/alignment.fasta
```

## Input notes

- The FASTA must already be a multiple sequence alignment with all sequences the same aligned length.
- Residue queries can be entered as `235`, `R235`, `Arg-235`, `R 235`, or as a marked motif such as `RFTH*R*ADTV`.
- Common pasted formatting artifacts such as smart quotes, non-breaking spaces, Unicode dashes, and Unicode asterisks are stripped/fixed automatically.
- If headers are not UniProt-style, the script will ask for species and gene name for each sequence before starting the lookup.

## Caveat
Written with the assistence of Codex/ChatGPT5.4; Tested against a few different structures, but there may still be bugs; let me know if you encounter any! I have provided a [sample alignment](ryr_alignment.fasta) for testing as part of the repo.
