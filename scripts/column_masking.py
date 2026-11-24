"""Generate a perturbed MSA from an AFsample2 run: mask the columns it chose.

AFsample2 records the alignment columns it masked under `X_msa_indexes` in its output
pickle. This script replays that choice onto an .a3m so the same perturbation can be fed
to another pipeline (AF3sample2, in our case).

Two things about a3m make the naive version of this wrong, and both are handled here:

  * the file may open with ColabFold's "#<lengths>\t<cardinalities>" line, which has to
    survive into the output or AlphaFold can no longer tell where one chain ends and the
    next begins;
  * a lowercase letter is an insertion relative to the query and occupies no alignment
    column, so alignment column k is the k-th non-lowercase character of a row, not its
    k-th character. Most rows of a real MSA carry at least one insertion.
"""
import argparse
import os
import pickle
from pathlib import Path


def read_a3m(path):
    """Return (comment_line_or_None, [(header, sequence), ...]).

    Sequences wrapped over several lines are joined back together.
    """
    comment, records, header, chunks = None, [], None, []
    for line in Path(path).read_text().splitlines():
        if not line.strip():
            continue
        if line.startswith("#") and header is None and not records:
            comment = line
        elif line.startswith(">"):
            if header is not None:
                records.append((header, "".join(chunks)))
            header, chunks = line, []
        else:
            chunks.append(line.strip())
    if header is not None:
        records.append((header, "".join(chunks)))
    return comment, records


def extract_pattern_from_pkl(pkl_path):
    """The alignment columns AFsample2 masked, as query-column indices."""
    with open(pkl_path, "rb") as fh:
        data = pickle.load(fh)
    columns = data.get("X_msa_indexes", [])
    return [int(c) for c in columns]


def match_column_offsets(seq):
    """Where each alignment column sits in this row's string.

    Uppercase letters and gaps are alignment columns; lowercase letters are insertions
    and are skipped, which is what keeps a row with insertions in register with the query.
    """
    return [i for i, c in enumerate(seq) if not c.islower()]


def apply_perturbation(records, columns, mask_gaps=False):
    """Mask the given alignment columns in every row but the query."""
    perturbed = [records[0]]          # the query is left untouched, as AFsample2 leaves it
    masked = 0
    for header, seq in records[1:]:
        chars = list(seq)
        offsets = match_column_offsets(seq)
        for col in columns:
            if 0 <= col < len(offsets):
                at = offsets[col]
                if chars[at] != "-" or mask_gaps:
                    chars[at] = "X"
                    masked += 1
        perturbed.append((header, "".join(chars)))
    return perturbed, masked


def write_a3m(path, comment, records):
    with open(path, "w") as fh:
        if comment:
            fh.write(comment + "\n")
        for header, seq in records:
            fh.write(header + "\n" + seq + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Generate a perturbed MSA from an .a3m and an AFsample2 pickle.")
    # both spellings work: the flags are what the README has always documented, the
    # positionals are what the first version of this script accepted
    parser.add_argument("msa_positional", nargs="?", metavar="msa_file",
                        help="Path to the input .a3m MSA file")
    parser.add_argument("pkl_positional", nargs="?", metavar="pkl_file",
                        help="Path to the AFsample2 .pkl holding X_msa_indexes")
    parser.add_argument("--msa_file", help="Path to the input .a3m MSA file")
    parser.add_argument("--pkl_file", help="Path to the AFsample2 .pkl file")
    parser.add_argument("--output", help="Output .a3m (default: <pkl stem>_perturbed_msa.a3m)")
    parser.add_argument("--mask_gaps", action="store_true",
                        help="Also turn gaps in a masked column into X (default: leave gaps)")
    args = parser.parse_args()

    msa_file = args.msa_file or args.msa_positional
    pkl_file = args.pkl_file or args.pkl_positional
    if not msa_file or not pkl_file:
        parser.error("give an MSA and a pickle, either positionally or with "
                     "--msa_file/--pkl_file")

    columns = extract_pattern_from_pkl(pkl_file)
    if not columns:
        print(f"Warning: no X_msa_indexes in {pkl_file}; the MSA will come out unchanged.")

    comment, records = read_a3m(msa_file)
    if not records:
        raise SystemExit(f"No sequences found in {msa_file}")

    query_columns = len(match_column_offsets(records[0][1]))
    out_of_range = [c for c in columns if c >= query_columns]
    if out_of_range:
        print(f"Warning: {len(out_of_range)} requested column(s) lie beyond the query's "
              f"{query_columns} columns and were skipped "
              f"(first one: {out_of_range[0]}). Is this the MSA that AFsample2 ran on?")

    perturbed, masked = apply_perturbation(records, columns, mask_gaps=args.mask_gaps)
    output = args.output or f"{os.path.splitext(os.path.basename(pkl_file))[0]}_perturbed_msa.a3m"
    write_a3m(output, comment, perturbed)

    print(f"{len(records)} sequences, query has {query_columns} alignment columns")
    print(f"masked {len(columns)} column(s) in {len(perturbed) - 1} sequences "
          f"-> {masked} residues replaced by X")
    print(f"Saved perturbed MSA to {output}")


if __name__ == "__main__":
    main()
