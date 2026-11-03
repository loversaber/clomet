"""clomet reproduces the 19 meta-clones of the HDBI LARRY paper (Fig 4).

The data are unpublished, so they are not in this repository. Set
    CLOMET_TEST_COUNTS  inputs/clone_counts_scaled.parquet of the manuscript repo
    CLOMET_TEST_LABELS  outputs/fig4_metaclones/Fig4_metaclone_labels.csv of the manuscript repo
to run the test; it is skipped otherwise (about 5 min).
"""
import os

import pandas as pd
import pytest

import clomet

counts_path = os.environ.get('CLOMET_TEST_COUNTS')
labels_path = os.environ.get('CLOMET_TEST_LABELS')
d3_prog = ['D3_HSC1', 'D3_HSC2', 'D3_MPP', 'D3_My-Ly progenitor', 'D3_GMP', 'D3_Monocyte progenitor',
           'D3_neutrophil prognitor', 'D3_EoBasoMast Precursor', 'D3_Mast cell progenitor', 'D3_MEMP', 'D3_CFU-E']


@pytest.mark.skipif(not (counts_path and labels_path), reason='paper data not available')
def test_paper_metaclones():
    counts = pd.read_parquet(counts_path)
    clones = counts[counts[d3_prog].gt(0).any(axis=1)].drop(columns=['D3_Mast cell', 'D3_Megakaryocyte'])
    res = clomet.assign(clones, k=14, split_ids=(3, 6, 7, 9, 10))
    expected = pd.read_csv(labels_path, index_col='clone_id', dtype={'meta_clone': str})
    assert len(res.meta_clone) == 2789
    assert res.meta_clone.nunique() == 19
    pd.testing.assert_series_equal(res.meta_clone.sort_index(), expected.meta_clone.sort_index(), check_names=False)
    pd.testing.assert_series_equal(res.clusters.sort_index(), expected.meta_clone_14.sort_index(),
                                   check_names=False, check_dtype=False)


def test_small_example():
    import numpy as np
    rng = np.random.default_rng(0)
    profiles = np.repeat(np.eye(3) * 20, 40, axis=0) + rng.poisson(1, (120, 3))
    clones = pd.DataFrame(np.hstack([profiles, rng.poisson(1, (120, 9))]),
                          index=[f'c{i}' for i in range(120)], columns=[f'f{j}' for j in range(12)])
    res = clomet.assign(clones, k=3, k_range=range(2, 7), n_pcs=5, n_neighbors=10, min_k=2)
    assert res.meta_clone.nunique() == 3
    truth = pd.Series(np.repeat([0, 1, 2], 40), index=clones.index)
    assert pd.crosstab(truth, res.meta_clone).max(axis=1).min() >= 36


def test_raw_counts_option():
    import numpy as np
    rng = np.random.default_rng(1)
    clones = pd.DataFrame(np.repeat(np.eye(2) * 15, 30, axis=0) + rng.poisson(1, (60, 2)),
                          columns=['a', 'b'], index=[f'c{i}' for i in range(60)])
    clones = pd.concat([clones, pd.DataFrame(rng.poisson(1, (60, 6)), index=clones.index)], axis=1)
    clones.columns = [str(c) for c in clones.columns]
    raw = clomet.leiden_scan(clones, n_pcs=5, n_neighbors=10, min_k=2, log1p=False)[0]
    logged = clomet.leiden_scan(clones, n_pcs=5, n_neighbors=10, min_k=2)[0]
    assert np.allclose(raw.X, clones.values)
    assert np.allclose(logged.X, np.log1p(clones.values))

