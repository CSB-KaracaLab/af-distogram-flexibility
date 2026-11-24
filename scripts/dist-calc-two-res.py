"""Measure the distance between two residues in a PDB file.

Point --input_pdb at a file or at a directory. A directory is walked recursively and
every .pdb under it is measured, which is how the per-case tables in data/pdb-distances
were built: one ensemble of predictions per folder, one row per structure.

Writes one row per model, so a multi-model PDB gives a distance distribution rather than
a single number.

Defaults to CB, as the paper does. Glycine has no CB, so a pair involving one has to be
measured on CA instead (--atom CA); the script says so rather than reporting zero.

The coordinates are read straight from the PDB's own columns, so this needs nothing but
the standard library. An earlier version drove PyMOL, which measured the same distances
but made the script awkward to install and impossible to pip.
"""
import argparse
import csv
import math
import os
import sys
from pathlib import Path

HEADER = ["PDB file", "Folder", "Residue 1", "Residue 2", "Chain 1", "Chain 2", "Atom",
          "State", "Distance"]


def read_models(path):
    """Coordinates per model: [{(chain, resseq_icode, atom_name): (x, y, z)}, ...].

    Columns are taken at their fixed PDB positions rather than by splitting on spaces,
    which runs together as soon as a coordinate reaches -100. Alternate locations keep
    the first one seen, which is what a viewer shows by default.
    """
    models, current = [], {}
    residue_names = {}
    seen_model_record = False
    for line in path.read_text(errors="ignore").splitlines():
        record = line[:6]
        if record == "MODEL ":
            seen_model_record = True
            current = {}
        elif record == "ENDMDL":
            models.append(current)
            current = {}
        elif record in ("ATOM  ", "HETATM"):
            name = line[12:16].strip()
            chain = line[21].strip() or " "
            resi = line[22:27].strip()            # residue number plus insertion code
            key = (chain, resi, name)
            if key not in current:                # first altloc wins
                try:
                    current[key] = (float(line[30:38]), float(line[38:46]),
                                    float(line[46:54]))
                except ValueError:
                    continue
            residue_names.setdefault((chain, resi), line[17:20].strip())
    if current or not seen_model_record:
        models.append(current)
    return [m for m in models if m], residue_names


def missing_atom(models, residue_names, resi, chain, atom):
    """Why this atom cannot be found, as a sentence - or None if every model has it."""
    key = (chain, str(resi), atom)
    if all(key in model for model in models):
        return None
    name = residue_names.get((chain, str(resi)))
    if name is None:
        return (f"residue {resi} of chain {chain} is not in the file; "
                "check the residue number and the chain")
    reason = f"no {atom} atom for residue {resi} of chain {chain} (residue is {name})"
    if name == "GLY" and atom == "CB":
        return reason + "; glycine has no CB, so measure this pair on CA with --atom CA"
    return reason + "; check the residue number, the chain, and --atom"


def measure(path, args):
    """Every model of one structure, as CSV rows. Returns (rows, problem_or_None)."""
    models, residue_names = read_models(path)
    if not models:
        return [], "no atoms in this file"

    for resi, chain in ((args.residue1, args.chain1), (args.residue2, args.chain2)):
        problem = missing_atom(models, residue_names, resi, chain, args.atom)
        if problem:
            return [], problem

    if args.state and args.state > len(models):
        return [], f"--state {args.state} but this structure has {len(models)} model(s)"
    wanted = range(1, len(models) + 1) if args.state == 0 else [args.state]

    key1 = (args.chain1, str(args.residue1), args.atom)
    key2 = (args.chain2, str(args.residue2), args.atom)
    rows = []
    for state in wanted:
        a = models[state - 1][key1]
        b = models[state - 1][key2]
        distance = math.dist(a, b)
        rows.append([path.name, path.resolve().parent.name, args.residue1, args.residue2,
                     args.chain1, args.chain2, args.atom, state, round(distance, 3)])
    return rows, None


def structures_under(target):
    """The .pdb files to measure: the one given, or every one beneath a directory."""
    path = Path(target)
    if path.is_dir():
        found = sorted(p for p in path.rglob("*.pdb") if p.is_file())
        if not found:
            sys.exit(f"Error: no .pdb files under {target}")
        return found
    if not path.is_file():
        sys.exit(f"Error: {target} is neither a file nor a directory")
    return [path]


def main():
    parser = argparse.ArgumentParser(
        description="Measure a residue pair across one structure or a directory of them.")
    parser.add_argument('--input_pdb', type=str, required=True,
                        help='A PDB file, or a directory to walk for .pdb files')
    parser.add_argument('--residue1', type=int, required=True, help='Residue number 1')
    parser.add_argument('--residue2', type=int, required=True, help='Residue number 2')
    parser.add_argument('--chain1', type=str, default='A', help='Chain ID for residue 1 (default: A)')
    parser.add_argument('--chain2', type=str, default='A', help='Chain ID for residue 2 (default: A)')
    parser.add_argument('--atom', type=str, default='CB',
                        help='Atom name to measure between (default: CB; use CA for glycine)')
    parser.add_argument('--state', type=int, default=0,
                        help='Model to measure, 1-based. 0, the default, measures every model.')
    parser.add_argument('--output_csv', type=str, default='measurements.csv',
                        help='Output CSV file (appended to)')
    parser.add_argument('--quiet', action='store_true',
                        help='Only report the summary, not progress')
    args = parser.parse_args()

    targets = structures_under(args.input_pdb)
    single = len(targets) == 1

    measured, skipped, distances = 0, [], []
    handle = writer = None
    try:
        for n, path in enumerate(targets, start=1):
            rows, problem = measure(path, args)
            if problem:
                if single:
                    sys.exit(f"Error: {path}: {problem}")
                skipped.append((path.name, problem))
                continue
            if writer is None:
                # opened only once there is something to write, so a run that measures
                # nothing does not leave an empty CSV behind
                write_header = not os.path.exists(args.output_csv)
                handle = open(args.output_csv, "a", newline="")
                writer = csv.writer(handle)
                if write_header:
                    writer.writerow(HEADER)
            writer.writerows(rows)
            measured += 1
            distances.extend(row[-1] for row in rows)
            if not single and not args.quiet and n % 1000 == 0:
                print(f"  {n}/{len(targets)} structures", flush=True)
    finally:
        if handle is not None:
            handle.close()

    label = f"{args.atom}{args.residue1}-{args.atom}{args.residue2}"
    print(f"{args.input_pdb}: {measured}/{len(targets)} structure(s) measured, {label}")
    if distances:
        print(f"  {len(distances)} measurement(s): min {min(distances):.3f} A, "
              f"max {max(distances):.3f} A, "
              f"mean {sum(distances) / len(distances):.3f} A")
    if skipped:
        print(f"  skipped {len(skipped)}: {skipped[0][1]}"
              + (f" (and {len(skipped) - 1} more)" if len(skipped) > 1 else ""))
    if measured:
        print(f"  written to {args.output_csv}")


if __name__ == "__main__":
    main()
