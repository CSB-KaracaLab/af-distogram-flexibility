"""Measure the distance between two residues in a PDB file, with PyMOL.

Writes one row per model. A PDB holding an ensemble - an AlphaFold sampling run, frames
pulled out of a trajectory - therefore gives a distance distribution rather than a single
number, which is what these measurements are for.

Defaults to CB, as the paper does. Glycine has no CB, so a pair involving one has to be
measured on CA instead (--atom CA); the script says so rather than reporting zero.
"""
import argparse
import csv
import os
import sys

import pymol

parser = argparse.ArgumentParser(
    description="Measure the distance between two residues, one row per model.")
parser.add_argument('--input_pdb', type=str, required=True, help='Input PDB file name')
parser.add_argument('--residue1', type=int, required=True, help='Residue number 1')
parser.add_argument('--residue2', type=int, required=True, help='Residue number 2')
parser.add_argument('--chain1', type=str, default='A', help='Chain ID for residue 1 (default: A)')
parser.add_argument('--chain2', type=str, default='A', help='Chain ID for residue 2 (default: A)')
parser.add_argument('--atom', type=str, default='CB',
                    help='Atom name to measure between (default: CB; use CA for glycine)')
parser.add_argument('--state', type=int, default=0,
                    help='Model to measure, 1-based. 0, the default, measures every model.')
parser.add_argument('--output_csv', type=str, default='measurements.csv', help='Output CSV file')
args = parser.parse_args()

pymol.finish_launching(['pymol', '-cq'])
pymol.cmd.load(args.input_pdb)

selection1 = f"chain {args.chain1} and resi {args.residue1} and name {args.atom}"
selection2 = f"chain {args.chain2} and resi {args.residue2} and name {args.atom}"

# An empty selection makes cmd.get_distance return 0.0, which looks exactly like a real
# measurement in the CSV. Check before measuring instead.
for label, selection, resi, chain in (("1", selection1, args.residue1, args.chain1),
                                      ("2", selection2, args.residue2, args.chain2)):
    found = pymol.cmd.count_atoms(selection)
    if found == 0:
        resn = pymol.cmd.get_model(f"chain {chain} and resi {resi} and name CA").atom
        name = resn[0].resn if resn else "not present"
        # cmd.quit() ends the process before any later line runs, and cmd.quit(1) hangs
        # here, so report the reason and leave immediately with a non-zero status
        print(f"Error: no {args.atom} atom for residue {resi} of chain {chain} "
              f"(residue is {name}) in {args.input_pdb}."
              + ("\nGlycine has no CB; measure this pair on CA with --atom CA."
                 if name == "GLY" else
                 "\nCheck the residue number, the chain, and --atom."),
              file=sys.stderr, flush=True)
        os._exit(1)
    if found > 1:
        print(f"Warning: selection {label} matches {found} atoms in one model "
              f"(alternate locations?); PyMOL will use the first.")

states = pymol.cmd.count_states()
wanted = range(1, states + 1) if args.state == 0 else [args.state]
if args.state and args.state > states:
    print(f"Error: --state {args.state}: {args.input_pdb} has {states} model(s).",
          file=sys.stderr, flush=True)
    os._exit(1)

header = ["PDB file", "Residue 1", "Residue 2", "Chain 1", "Chain 2", "Atom", "State",
          "Distance"]
write_header = not os.path.exists(args.output_csv)
rows = []
for state in wanted:
    distance = pymol.cmd.get_distance(selection1, selection2, state=state)
    rows.append([os.path.basename(args.input_pdb), args.residue1, args.residue2,
                 args.chain1, args.chain2, args.atom, state, round(distance, 3)])

with open(args.output_csv, "a", newline="") as csv_file:
    writer = csv.writer(csv_file)
    if write_header:
        writer.writerow(header)
    writer.writerows(rows)

measured = [r[-1] for r in rows]
print(f"{args.input_pdb}: {len(measured)} model(s), {args.atom}{args.residue1}-"
      f"{args.atom}{args.residue2}")
if len(measured) > 1:
    print(f"  min {min(measured):.3f} A, max {max(measured):.3f} A, "
          f"mean {sum(measured) / len(measured):.3f} A")
else:
    print(f"  {measured[0]:.3f} A")
print(f"  appended {len(rows)} row(s) to {args.output_csv}")

pymol.cmd.delete("all")
pymol.cmd.quit()
