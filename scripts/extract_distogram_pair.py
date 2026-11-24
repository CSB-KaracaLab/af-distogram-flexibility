"""Extract the distogram of one residue pair from an AlphaFold prediction.

Reads the distogram logits from a .pkl, .pickle or .npz output, turns them into a
probability distribution over distance bins, and writes that distribution for a single
residue pair as a CSV.

The distance axis follows AlphaFold2's own distogram head: 64 bins whose breaks run from
2.3125 A to 21.6875 A. The first and last bins are open-ended (everything closer than
2.3125 A, everything further than 21.6875 A); they are reported at half a bin width beyond
the outermost break, which is the usual convention.
"""
import argparse
import os
import pickle

import numpy as np
import pandas as pd

FIRST_BREAK, LAST_BREAK, NUM_BINS = 2.3125, 21.6875, 64


def bin_centers():
    breaks = np.linspace(FIRST_BREAK, LAST_BREAK, NUM_BINS - 1)
    width = breaks[1] - breaks[0]
    return np.concatenate([[breaks[0] - width / 2],
                           (breaks[:-1] + breaks[1:]) / 2,
                           [breaks[-1] + width / 2]])


def load_logits(file_path):
    """Find the distogram logits in an AlphaFold output file."""
    if file_path.endswith(".npz"):
        data = np.load(file_path, allow_pickle=True)
        logits = data.get("logits", None)
        if logits is None:
            dgram = data.get("distogram", None)
            if isinstance(dgram, np.ndarray) and dgram.dtype == object:
                dgram = dgram.item()
            if isinstance(dgram, dict) and "logits" in dgram:
                logits = dgram["logits"]
            elif isinstance(dgram, np.ndarray):
                logits = dgram
    elif file_path.endswith((".pkl", ".pickle")):
        with open(file_path, "rb") as fh:
            data = pickle.load(fh)
        logits = data.get("logits", None)
        if logits is None:
            dgram = data.get("distogram", None)
            if isinstance(dgram, dict) and "logits" in dgram:
                logits = dgram["logits"]
            elif isinstance(dgram, np.ndarray):
                logits = dgram
    else:
        raise SystemExit(f"Unsupported file type: {file_path} "
                         "(expected .npz, .pkl or .pickle)")
    if logits is None:
        raise SystemExit("Could not locate 'logits' or 'distogram' in the file.")
    return np.asarray(logits)


def save_distogram_to_csv(file_path, res_a, res_b, temperature=1.0, output=None):
    """Write the distance distribution for one residue pair to a CSV."""
    if not os.path.exists(file_path):
        raise SystemExit(f"File not found: {file_path}")
    # residue numbering is 1-based; letting 0 or a negative through would silently wrap
    # round to the far end of the chain
    for name, value in (("res_a", res_a), ("res_b", res_b)):
        if value < 1:
            raise SystemExit(f"{name}={value}: residue numbers are 1-based, so they "
                             "must be 1 or greater.")

    logits = load_logits(file_path)
    if logits.ndim != 3 or logits.shape[0] != logits.shape[1]:
        raise SystemExit(f"Expected a square L x L x bins distogram, got {logits.shape}.")
    if logits.shape[2] != NUM_BINS:
        print(f"Warning: this distogram has {logits.shape[2]} bins, not {NUM_BINS}. "
              "The distance axis below assumes AlphaFold2's 64-bin head.")

    length = logits.shape[0]
    if res_a > length or res_b > length:
        raise SystemExit(f"Residue {max(res_a, res_b)} is out of bounds: this distogram "
                         f"covers {length} residues. Note that for a complex the "
                         "distogram runs over the chains concatenated, so a residue "
                         "number from the second chain is not its index here.")

    scaled = logits / temperature
    exp = np.exp(scaled - np.max(scaled, axis=-1, keepdims=True))
    probs = exp / np.sum(exp, axis=-1, keepdims=True)
    pair_probs = probs[res_a - 1, res_b - 1, :]

    centers = bin_centers()
    if len(centers) != len(pair_probs):           # a non-AF2 bin count: report indices
        centers = np.arange(len(pair_probs), dtype=float)

    df = pd.DataFrame({"Distance_Angstrom": centers, "Probability": pair_probs})
    if output is None:
        base = os.path.splitext(os.path.basename(file_path))[0]
        output = os.path.join(os.path.dirname(file_path) or ".",
                              f"{base}_res{res_a}_{res_b}.csv")
    df.to_csv(output, index=False)

    expected = float(np.sum(centers * pair_probs))
    print(f"{length} residues, pair {res_a}-{res_b}: "
          f"most likely {centers[int(np.argmax(pair_probs))]:.3f} A, "
          f"expected {expected:.3f} A")
    print(f"Successfully saved distogram data to: {output}")


def main():
    parser = argparse.ArgumentParser(
        description="Extract the distogram of one residue pair from an AlphaFold output.")
    # the flag spelling is what the README documents; the positional one is what the
    # first version of this script accepted
    parser.add_argument("file_positional", nargs="?", metavar="input_file",
                        help="AlphaFold output: .npz, .pkl or .pickle")
    parser.add_argument("res_a", type=int, help="First residue number (1-based)")
    parser.add_argument("res_b", type=int, help="Second residue number (1-based)")
    parser.add_argument("--input_pickle", "--input_file", dest="input_file",
                        help="AlphaFold output: .npz, .pkl or .pickle")
    parser.add_argument("--temperature", type=float, default=1.0,
                        help="Softmax temperature (default: 1.0)")
    parser.add_argument("--output", help="Output CSV (default: next to the input file)")
    args = parser.parse_args()

    file_path = args.input_file or args.file_positional
    if not file_path:
        parser.error("give an input file, either positionally or with --input_pickle")
    save_distogram_to_csv(file_path, args.res_a, args.res_b,
                          temperature=args.temperature, output=args.output)


if __name__ == "__main__":
    main()
