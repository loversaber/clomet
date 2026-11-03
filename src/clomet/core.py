"""Consensus clustering of clones into meta-clones.

Clones are described by their cell counts in a set of features (e.g. timepoint_celltype columns). The clones are
clustered with Leiden at many resolutions; the consensus matrix (how often two clones share a cluster) gives an
average-linkage tree, which is cut into k clusters; chosen clusters can then be split in two. Final groups are
numbered 1..n and called meta-clones.
"""
from itertools import combinations
from collections import Counter
from types import SimpleNamespace

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import silhouette_score, normalized_mutual_info_score, adjusted_rand_score
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform


def leiden_scan(clones, n_pcs=40, n_neighbors=25, metric='cosine', resolutions=np.round(np.arange(0.2, 3.2, 0.1), 10),
                min_k=5, seed=0, log1p=True):
    """log1p (optional) -> PCA -> kNN -> Leiden at each resolution. Returns the AnnData (one leiden_<r> column per
    resolution), the scan table (k, silhouette, NMI to the previous resolution) and the partitions used (best silhouette
    per k)."""
    x = clones.fillna(0).values.astype(float)
    adata = sc.AnnData(np.log1p(x) if log1p else x, obs=clones.index.to_frame(name='clone_id'))
    sc.tl.pca(adata, n_comps=n_pcs, random_state=seed)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, metric=metric)

    rows, prev = [], None
    for r in resolutions:
        key = f'leiden_{r:.2f}'
        sc.tl.leiden(adata, resolution=r, key_added=key, random_state=seed, flavor='leidenalg')
        lab = adata.obs[key].astype(int).to_numpy()
        k = np.unique(lab).size
        rows.append(dict(resolution=r, k=k,
                         silhouette=silhouette_score(adata.obsm['X_pca'], lab, metric='cosine') if k > 1 else np.nan,
                         nmi_to_prev=normalized_mutual_info_score(prev, lab) if prev is not None else np.nan))
        prev = lab
    scan = pd.DataFrame(rows)
    used = (scan[scan.k >= min_k].sort_values('silhouette', ascending=False).drop_duplicates('k')
            .sort_values('k').reset_index(drop=True))
    return adata, scan, used


def consensus(adata, used):
    """Fraction of the used partitions in which two clones share a cluster, and the average-linkage tree on 1 - it."""
    lab_matrix = np.vstack([adata.obs[f'leiden_{r:.2f}'].astype(int).to_numpy() for r in used.resolution])
    co = np.zeros((adata.n_obs, adata.n_obs), dtype=np.float32)
    for lab in lab_matrix:
        co += (lab[:, None] == lab[None, :]).astype(np.float32)
    co /= len(lab_matrix)
    ids = adata.obs.clone_id.values
    coassoc = pd.DataFrame(co, index=ids, columns=ids)

    dist = 1.0 - co
    np.fill_diagonal(dist, 0.0)
    tree = linkage(squareform(dist, checks=False), method='average')
    return coassoc, tree


def cut_tree(tree, index, k_range):
    """Tree cuts at each k (columns k_<k>) and a summary of the cluster sizes."""
    cuts = pd.DataFrame({f'k_{k}': fcluster(tree, t=k, criterion='maxclust') for k in k_range}, index=index)
    sizes = {k: cuts[f'k_{k}'].value_counts() for k in k_range}
    k_summary = pd.DataFrame({'k': k_range, 'n_clusters': [len(sizes[k]) for k in k_range],
                              'min_size': [sizes[k].min() for k in k_range], 'max_size': [sizes[k].max() for k in k_range],
                              'mean_size': [sizes[k].mean() for k in k_range],
                              'cv_size': [sizes[k].std() / sizes[k].mean() for k in k_range]})
    return cuts, k_summary


def cut_ari(cuts, k_range):
    return pd.DataFrame([[adjusted_rand_score(cuts[f'k_{a}'], cuts[f'k_{b}']) for b in k_range] for a in k_range],
                        index=k_range, columns=k_range)


def best_split(members, coassoc, max_k=15, min_balance=0.15):
    """Best two-group split of `members`: for every cut k = 2..max_k of their own consensus tree and every pair of
    branches, the branches seed groups A and B and the other clones join the seed they are closer to on average.
    Score = separation (mean within - mean between similarity) x sqrt(balance); balance = smaller / larger group."""
    sim = coassoc.loc[members, members]
    d = 1 - sim.values
    np.fill_diagonal(d, 0)
    z = linkage(squareform((d + d.T) / 2), method='average')

    def within(g):
        if len(g) == 1:
            return 1.0
        w = sim.loc[g, g].values.copy()
        np.fill_diagonal(w, np.nan)
        return np.nanmean(w)

    tried = []
    for k in range(2, max_k + 1):
        lab = fcluster(z, t=k, criterion='maxclust')
        for b1, b2 in combinations(list(Counter(lab).keys()), 2):
            seed1 = [c for c, l in zip(members, lab) if l == b1]
            seed2 = [c for c, l in zip(members, lab) if l == b2]
            ga, gb = seed1.copy(), seed2.copy()
            for c in [c for c, l in zip(members, lab) if l not in (b1, b2)]:
                (ga if sim.loc[c, seed1].mean() >= sim.loc[c, seed2].mean() else gb).append(c)
            if len(ga) < len(gb):
                ga, gb = gb, ga
            balance = len(gb) / len(ga)
            separation = (within(ga) + within(gb)) / 2 - sim.loc[ga, gb].values.mean()
            tried.append(dict(k=k, a=ga, b=gb, balance=balance, separation=separation,
                              score=separation * np.sqrt(balance)))
    valid = [t for t in tried if t['balance'] >= min_balance]
    best = max(valid or tried, key=lambda t: t['score'])
    return SimpleNamespace(best=best, tried=tried, z=z, members=members, sim=sim)


def split_clusters(clusters, coassoc, split_ids, max_k=15, min_balance=0.15):
    """Split each cluster in `split_ids` into A (larger) and B, then number all groups 1..n in order
    (1, 2, 3A, 3B, 4, ...). Returns the new labels, the renumbering, a table of the splits and the split details."""
    n = clusters.max()
    groups = [str(c) for c in range(1, n + 1) if c not in split_ids] + [f'{c}{s}' for c in split_ids for s in 'AB']
    renumber = {g: str(i + 1) for i, g in enumerate(sorted(groups, key=lambda g: (int(g.rstrip('AB')), g)))}

    labels = clusters.astype(str).map(renumber).rename('meta_clone')
    rows, details = [], {}
    for cid in split_ids:
        members = clusters[clusters == cid].index.tolist()
        sp = best_split(members, coassoc, max_k, min_balance)
        details[cid] = sp
        labels[sp.best['a']] = renumber[f'{cid}A']
        labels[sp.best['b']] = renumber[f'{cid}B']
        rows.append(dict(meta_clone_14=cid, clones=len(members), splits_tried=len(sp.tried), best_k=sp.best['k'],
                         meta_clone_A=renumber[f'{cid}A'], n_A=len(sp.best['a']), meta_clone_B=renumber[f'{cid}B'],
                         n_B=len(sp.best['b']), balance=sp.best['balance'], separation=sp.best['separation'],
                         score=sp.best['score']))
    assert labels.notna().all()
    return labels, renumber, pd.DataFrame(rows).rename(columns={'meta_clone_14': f'meta_clone_{n}'}), details


def assign(clones, k, split_ids=(), k_range=range(4, 21, 2), **scan_kw):
    """The whole pipeline: `clones` is a clones x features table of counts, `k` the tree cut, `split_ids` the
    clusters of that cut to split in two, `k_range` the cuts to compare; other keywords (n_pcs, n_neighbors, metric,
    resolutions, min_k, seed, log1p) go to leiden_scan. Returns a namespace with adata, scan, used, coassoc, tree, cuts, k_summary, ari,
    clusters (tree cut at k), meta_clone (after the splits), renumber, splits and split_details."""
    k_range = list(k_range)
    adata, scan, used = leiden_scan(clones, **scan_kw)
    coassoc, tree = consensus(adata, used)
    cuts, k_summary = cut_tree(tree, clones.index, k_range)
    clusters = pd.Series(fcluster(tree, t=k, criterion='maxclust'), index=clones.index, name=f'meta_clone_{k}')
    meta_clone, renumber, splits, details = split_clusters(clusters, coassoc, list(split_ids))
    meta_clone.name = 'meta_clone'
    return SimpleNamespace(adata=adata, scan=scan, used=used, coassoc=coassoc, tree=tree, k_range=k_range, cuts=cuts,
                           k_summary=k_summary, ari=cut_ari(cuts, k_range), k=k, clusters=clusters,
                           meta_clone=meta_clone, renumber=renumber, splits=splits, split_details=details)


assign_metaclones = assign
