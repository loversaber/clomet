"""clomet on simulated clones."""
import numpy as np
import pandas as pd

import clomet


def test_small_example():
    rng = np.random.default_rng(0)
    profiles = np.repeat(np.eye(3) * 20, 40, axis=0) + rng.poisson(1, (120, 3))
    clones = pd.DataFrame(np.hstack([profiles, rng.poisson(1, (120, 9))]),
                          index=[f'c{i}' for i in range(120)], columns=[f'f{j}' for j in range(12)])
    res = clomet.assign(clones, k=3, k_range=range(2, 7), n_pcs=5, n_neighbors=10, min_k=2)
    assert res.meta_clone.nunique() == 3
    truth = pd.Series(np.repeat([0, 1, 2], 40), index=clones.index)
    assert pd.crosstab(truth, res.meta_clone).max(axis=1).min() >= 36


def test_raw_counts_option():
    rng = np.random.default_rng(1)
    clones = pd.DataFrame(np.repeat(np.eye(2) * 15, 30, axis=0) + rng.poisson(1, (60, 2)),
                          columns=['a', 'b'], index=[f'c{i}' for i in range(60)])
    clones = pd.concat([clones, pd.DataFrame(rng.poisson(1, (60, 6)), index=clones.index)], axis=1)
    clones.columns = [str(c) for c in clones.columns]
    raw = clomet.leiden_scan(clones, n_pcs=5, n_neighbors=10, min_k=2, log1p=False)[0]
    logged = clomet.leiden_scan(clones, n_pcs=5, n_neighbors=10, min_k=2)[0]
    assert np.allclose(raw.X, clones.values)
    assert np.allclose(logged.X, np.log1p(clones.values))
