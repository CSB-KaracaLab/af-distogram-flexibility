"""Measure the distance between two residues in a PDB file, with PyMOL.

Point --input_pdb at a file or at a directory. A directory is walked recursively and
every .pdb under it is measured, which is how the per-case tables in data/pdb-distances
were built: one ensemble of predictions per folder, one row per structure.

Writes one row per model, so a multi-model PDB gives a distance distribution rather than
a single number.

Defaults to CB, as the paper does. Glycine has no CB, so a pair involving one has to be
measured on CA instead (--atom CA); the script says so rather than reporting zero.
"""
import argparse
import csv
import os
import sys
from pathlib import Path

import pymol

HEADER = ["PDB file", "Folder", "Residue 1", "Residue 2", "Chain 1", "Chain 2", "Atom",
          "State", "Distance"]


def structures_under(target):
    """The .pdb files to measure: the one given, or every one beneath a directory."""
    path = Path(target)
    if path.is_dir():
        found = sorted(p for p in path.rglob("*.pdb") if p.is_file())
        if not found:
            print(f"Error: no .pdb files under {target}", file=sys.stderr, flush=True)
            os._exit(1)
        return found
    if not path.is_file():
        print(f"Error: {target} is neither a file nor a directory",
              file=sys.stderr, flush=True)
        os._exit(1)
    return [path]


def missing_atom(selection, resi, chain, atom):
    """Why this selection is empty, as a sentence - or None if it is fine.

    An empty selection makes cmd.get_distance return 0.0, which looks exactly like a real
    measurement once it is in the CSV, so it has to be caught before measuring.
    """
    if pymol.cmd.count_atoms(selection) > 0:
        return None
    present = pymol.cmd.get_model(f"chain {chain} and resi {resi} and name CA").atom
    name = present[0].resn if present else "not present"
    reason = (f"no {atom} atom for residue {resi} of chain {chain} (residue is {name})")
    if name == "GLY":
        return reason + "; glycine has no CB, so measure this pair on CA with --atom CA"
    return reason + "; check the residue number, the chain, and --atom"


def measure(path, args):
    """Every model of one structure, as CSV rows. Returns (rows, problem_or_None)."""
    pymol.cmd.delete("all")
    pymol.cmd.load(str(path), "subject")
    selection1 = f"chain {args.chain1} and resi {args.residue1} and name {args.atom}"
    selection2 = f"chain {args.chain2} and resi {args.residue2} and name {args.atom}"

    for selection, resi, chain in ((selection1, args.residue1, args.chain1),
                                   (selection2, args.residue2, args.chain2)):
        problem = missing_atom(selection, resi, chain, args.atom)
        if problem:
            return [], problem

    states = pymol.cmd.count_states("subject")
    if args.state and args.state > states:
        return [], f"--state {args.state} but this structure has {states} model(s)"
    wanted = range(1, states + 1) if args.state == 0 else [args.state]

    rows = []
    for state in wanted:
        distance = pymol.cmd.get_distance(selection1, selection2, state=state)
        rows.append([path.name, path.parent.name, args.residue1, args.residue2,
                     args.chain1, args.chain2, args.atom, state, round(distance, 3)])
    return rows, None


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
                        help='Only report the summary, not every structure')
    args = parser.parse_args()

    targets = structures_under(args.input_pdb)
    single = len(targets) == 1

    pymol.finish_launching(['pymol', '-cq'])

    write_header = not os.path.exists(args.output_csv)
    measured, skipped, distances = 0, [], []
    with open(args.output_csv, "a", newline="") as handle:
        writer = csv.writer(handle)
        if write_header:
            writer.writerow(HEADER)
        for n, path in enumerate(targets, start=1):
            rows, problem = measure(path, args)
            if problem:
                if single:
                    print(f"Error: {path}: {problem}", file=sys.stderr, flush=True)
                    os._exit(1)
                skipped.append((path.name, problem))
                continue
            writer.writerows(rows)
            measured += 1
            distances.extend(row[-1] for row in rows)
            if not single and not args.quiet and n % 250 == 0:
                print(f"  {n}/{len(targets)} structures", flush=True)

    label = f"{args.atom}{args.residue1}-{args.atom}{args.residue2}"
    print(f"{args.input_pdb}: {measured}/{len(targets)} structure(s) measured, {label}")
    if distances:
        print(f"  {len(distances)} measurement(s): min {min(distances):.3f} A, "
              f"max {max(distances):.3f} A, "
              f"mean {sum(distances) / len(distances):.3f} A")
    if skipped:
        print(f"  skipped {len(skipped)}: {skipped[0][1]}"
              + (f" (and {len(skipped) - 1} more)" if len(skipped) > 1 else ""))
    print(f"  written to {args.output_csv}")

    pymol.cmd.delete("all")
    pymol.cmd.quit()


if __name__ == "__main__":
    main()
