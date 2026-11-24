<img src="logo.png" alt="AF-distogram Logo" width="600"/>


Here we explore **flexibility signals** embedded in **AlphaFold distograms**, with a special focus on **predicting binding-induced hinge motion** interpretation.  

This repository supports the manuscript:
[**"Exploring the Potential of AlphaFold Distograms for Predicting Binding-induced Hinge Motions"**](https://doi.org/10.1002/1873-3468.70297)  
Authors: Büşra Savaş, Ayşe Berçin Barlas, Ezgi Karaca  

---

## 🔬 Research Background

AlphaFold’s distograms encode inter-residue distance distributions that often contain conformational information not visible in the final coordinate predictions. In this project, we investigate how these distograms capture binding-induced hinge motions using the AK2/AIFM1 complex as a test case. We compare distogram outputs across multiple AlphaFold versions (AF2.0 (AF2-ptm), AF2.1, AF2.2, AF2.3, AF3), different enhanced sampling strategies (MinnieFold, AF_cluster, AFsample2, AF3sample2, AF3_cluster), and evaluate their correspondence with molecular dynamics simulations. We also compare structural outputs from AF2- and AF3-based enhanced sampling methods with generative approaches such as BioEmu, Chai-1, and Boltz-2 in terms of sampled hinge distances. Furthermore, we expanded our research for three additional protein systems: IL1R1, MIA40, and EIF2B3 to showcase the ability of AF distograms to pinpoint binding-induced hinge motions. 

---

## 📂 Folder Structure
```
data*/                                        
├── distograms**/                             # Generated pairwise distograms for our case set, including AK2, IL1R1, MIA40, and EIF2B3.
├── pdb-distances**/                          # Calculated pairwise distance measurements.
├── msas/                                     # MSA files (.a3m) used in distogram generation
├── MD/                                       # Molecular dynamics based reference distances

scripts/
├── dist-calc-two-res.py                      # Extracts pairwise distances from PDB files
├── extract_distogram_pair.py                 # Extracts pairwise distograms from .pkl, .pickle or .npz files
└── column_masking.py                         # Generates perturbed MSAs from AFsample2 .pkl outputs and the original MSA
```
\* PDB files omitted from the repo; available via [Zenodo.](https://zenodo.org/records/17638774)

\** Each csv file comes with a corresponding txt file providing a description of the run parameters for each data.
  
## 🧰 Requirements

```
python >= 3.8
numpy
pandas
```

`dist-calc-two-res.py` additionally needs **PyMOL**, which is not installable with pip.
Either use the open-source build:

```bash
conda install -c conda-forge pymol-open-source
```

or run the script with the Python inside an existing PyMOL installation, e.g. on macOS
`/Applications/PyMOL.app/Contents/bin/python3`.

The other two scripts need only numpy and pandas:

```bash
pip install -r requirements.txt
```

---

**Calculate custom pairwise distances**

Extract pairwise distances for any residue pair (e.g., 230–233):
```
python dist-calc-two-res.py --input_pdb /path/to/pdb --residue1 230 --residue2 233
```

One row is written per model, so a multi-model PDB — an ensemble, or frames taken from a
trajectory — gives a distance distribution rather than a single number. Measurements are
taken between CB atoms by default; glycine has no CB, so a pair involving one has to be
measured on CA with `--atom CA`. Use `--chain1` / `--chain2` for an inter-chain pair and
`--state` to pick a single model.

**Retrieve custom pairwise distograms**

Extract pairwise distograms for any residue pair (e.g., 230–233), from .npz, .pickle and
.pkl files:
```
python extract_distogram_pair.py --input_pickle /path/to/pickle 230 233
```

Residue numbers are 1-based and index the distogram directly. For a **complex** the
distogram runs over the chains concatenated, so a residue number taken from the second
chain is not its index here — add the length of the preceding chains first.

**Generate Perturbed MSAs for AF3sample2 pipeline**

Create perturbed MSAs using AFsample2 output pickle files:
```
python column_masking.py --msa_file path/to/input.msa.a3m --pkl_file path/to/afsample2_output.pkl
```

The indices in the pickle are alignment columns of the query. Lowercase letters in an .a3m
are insertions and occupy no alignment column, so they are skipped when the columns are
located in each row; the query sequence and any `#` header line are passed through
unchanged. Gaps falling in a masked column are left as gaps unless `--mask_gaps` is given.

> Each script also accepts its arguments positionally, as earlier versions did.

**Reproducing the deposited tables**

The CSVs under `data/distograms/` and `data/pdb-distances/` collect many runs into one
file — one column per prediction, or one row per structure. The scripts here produce a
single pair from a single file; the per-case tables were assembled from many such calls.

**🧠 Citing This Work**

If you find this repository useful, please consider citing our paper.

**📧 Contact**

For questions or contributions, feel free to contact:

Ezgi Karaca
📩 ezgi.karaca@ibg.edu.tr
