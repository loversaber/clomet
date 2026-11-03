<p align="center">
  <img src="docs/logo/clomet_logo_text_light.svg" alt="clomet: CLOnes to METa-clones. Scissors cut a consensus tree into three meta-clones" width="560">
</p>

# clomet

**CLOnes to METa-clones.** Consensus clustering of lineage-traced clones (e.g. LARRY barcodes) by their cell-type
output.

A clone is described by its cell counts in a set of features, for example `timepoint_celltype` columns
(`D3_HSC1`, `D10_Monocyte`, …). clomet groups clones with similar output into meta-clones:

1. **Leiden scan:** log1p of the counts → PCA → k-nearest-neighbour graph → Leiden at many resolutions (0.2–3.1 by
   default). Resolutions that give the same number of clusters are reduced to the one with the best silhouette.
2. **Consensus matrix:** for each pair of clones, the fraction of these partitions in which they share a cluster.
3. **Tree:** average-linkage tree on 1 − consensus, cut into `k` clusters. Cuts over a range of `k` are compared
   (cluster sizes, adjusted Rand index, Sankey of how clones move between cuts) to choose `k`.
4. **Splits (optional):** chosen clusters are split in two, using the best split of their own consensus tree
   (score = separation × √balance).
5. The final groups are numbered 1…n: these are the meta-clones.

The method follows Gao, Barile et al. (Nat Commun 2025).

## Install

```bash
git clone git@github.com:loversaber/clomet.git
cd clomet
pip install -e .            # add [sankey] for the Sankey plot, [test] for the tests
```

## Use

```python
import clomet

# clones: DataFrame, one row per clone, one column per feature (cell counts)
res = clomet.assign(clones, k=14, split_ids=(3, 6, 7, 9, 10))

res.meta_clone      # meta-clone of each clone, '1' ... '19'
res.clusters        # cluster at the tree cut k, before the splits
res.coassoc         # consensus matrix (clones x clones)
res.scan, res.used  # all Leiden resolutions, and the partitions used for the consensus
res.k_summary       # cluster sizes for each cut in k_range
res.splits          # the splits: sizes, balance, separation, score
```

To choose `k` and the clusters to split, look at the diagnostics first:

```python
res = clomet.assign(clones, k=14)                       # no splits yet
clomet.pl.resolution_scan(res.scan, res.used)
clomet.pl.choose_k(res.cuts, res.k_summary, res.k_range, chosen=14)
clomet.pl.cut_ari(res.ari)
clomet.pl.sankey(res.cuts, res.k_range).show()           # needs clomet[sankey]
clomet.pl.cut_clustermaps(res.coassoc, res.tree, res.cuts, [10, 12, 14, 16], 'clustermaps.pdf')
```

and after the splits:

```python
for cid, sp in res.split_details.items():
    clomet.pl.split(cid, sp, res.renumber, n_groups=res.k)
clomet.pl.coassociation(res.coassoc, res.tree, res.meta_clone, colors)   # colors: {meta-clone: colour}
```

Each step is also a function of its own (`leiden_scan`, `consensus`, `cut_tree`, `cut_ari`, `best_split`,
`split_clusters`) if you want to run them one at a time.

### Main parameters

| parameter | default | meaning |
|---|---|---|
| `k` | (required) | number of clusters at the tree cut |
| `split_ids` | `()` | clusters of the cut to split in two |
| `k_range` | `range(4, 21, 2)` | cuts compared in the diagnostics |
| `n_pcs`, `n_neighbors`, `metric` | 40, 25, `'cosine'` | PCA and neighbour graph |
| `resolutions` | 0.2–3.1 by 0.1 | Leiden resolutions |
| `min_k` | 5 | smallest number of clusters for a partition to enter the consensus |
| `seed` | 0 | random seed for PCA and Leiden |

## Tests

```bash
pytest tests
```

`test_small_example` runs on simulated clones. `test_paper_metaclones` checks that clomet reproduces the 19 meta-clones
of the HDBI LARRY paper; it needs the (unpublished) data and is skipped without it.

## Citation

Paper in preparation.
