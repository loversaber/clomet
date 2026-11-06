# Cord blood LARRY data from CLADES

Input for `examples/cord_blood_clades.ipynb`: human cord blood LARRY lineage tracing (days 3, 10 and 17), published with

> Gao M, Barile M, Chabra S, Haltalli M, Calderbank EF, Chao Y, Zheng W, Wilson NK, Laurenti E, Göttgens B, Huang Y.
> CLADES: a hybrid NeuralODE-Gillespie approach for unveiling clonal cell fate and differentiation dynamics.
> Nat Commun 16 (2025). https://doi.org/10.1038/s41467-025-63150-6

Raw data: GEO GSE276896. Processed data: figshare https://doi.org/10.6084/m9.figshare.27908142 (`CordBlood_Refine`),
licence CC BY 4.0. The files here are derived from the processed data.

| file | content |
|---|---|
| `clone_counts.csv` | 606 clones x 31 `Day<t>_<cell type>` columns: number of barcoded GFP+ cells of each clone |
| `clades_metaclones.csv` | meta-clone of each clone in CLADES (0–11) |
| `prepare_input.py` | makes the two files from the figshare AnnData and `CordBlood_Obs_Clones.pkl` |

The clone table is made as in CLADES (`demo/CordBlood_Refine/prepare_input.ipynb` in
https://github.com/StatBiomed/clonaltrans): GFP+ cells with a barcode, cells per clone and timepoint_celltype, clones
with at least 3 cells, columns with at least one cell.
