"""Make the clomet input from the CLADES cord blood data (figshare 10.6084/m9.figshare.27908142, CordBlood_Refine).

    python prepare_input.py <folder with the cord blood AnnData and CordBlood_Obs_Clones.pkl> [AnnData file name]

As in CLADES (demo/CordBlood_Refine/prepare_input.ipynb): GFP+ cells with a barcode, cells per clone and
timepoint_celltype, clones with at least 3 cells, columns with at least one cell. Writes clone_counts.csv and
clades_metaclones.csv next to this script.
"""
import sys
from pathlib import Path

import pandas as pd
import anndata as ad

src = Path(sys.argv[1])
h5ad = sys.argv[2] if len(sys.argv) > 2 else 'adata_cord_blood_as_submitted.h5ad'
here = Path(__file__).parent

obs = ad.read_h5ad(src / h5ad, backed='r').obs.copy()
data = pd.read_pickle(src / 'CordBlood_Obs_Clones.pkl')
barcodes = data[data.columns[14:]].fillna(0).loc[obs.index]

pos = obs[(obs.GFP == 'GFP+') & (obs.clones.astype(str) != 'Clone_')]
mat = barcodes.loc[pos.index].assign(time_label_full=pos.time_label_full.astype(str))
mat = mat.groupby('time_label_full').sum().T
mat = mat[mat.sum(axis=1) >= 3]
mat = mat.loc[:, mat.sum() > 0].astype(int)
mat.index.name, mat.columns.name = 'clone', None

clades = (obs.loc[obs['Meta clones'].astype(str) != '-1.0'].groupby('clones', observed=True)['Meta clones']
          .first().astype(str).str.replace(r'\.0$', '', regex=True).reindex(mat.index).rename('clades_meta_clone'))
assert clades.notna().all()

mat.to_csv(here / 'clone_counts.csv')
clades.to_csv(here / 'clades_metaclones.csv', index_label='clone')
print(f'{mat.shape[0]} clones x {mat.shape[1]} columns, {clades.nunique()} CLADES meta-clones')
