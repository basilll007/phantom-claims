"""Tissue = world state: cells with positions (um), cell types, and expression."""
import json
from dataclasses import dataclass, field
import numpy as np
import anndata as ad
import scipy.sparse as sp
from scipy.spatial import cKDTree

EXCLUDE_TYPES = {"Unlabeled"}


@dataclass
class Tissue:
    name: str
    xy: np.ndarray            # (n, 2) float64 um
    types: np.ndarray         # (n,) str
    counts: sp.csr_matrix     # (n, g) real panel counts
    genes: np.ndarray         # (g,)
    lib: np.ndarray           # (n,) library size (real panel)
    _tree: cKDTree = field(default=None, repr=False)

    @property
    def n(self):
        return self.xy.shape[0]

    @property
    def tree(self):
        if self._tree is None:
            self._tree = cKDTree(self.xy)
        return self._tree

    def type_list(self, min_cells=0):
        u, c = np.unique(self.types, return_counts=True)
        return [t for t, k in zip(u, c) if t not in EXCLUDE_TYPES and k >= min_cells]

    def lognorm(self, counts):
        """log1p(counts / lib * median(lib)); counts: (n, k) dense or sparse."""
        scale = np.median(self.lib) / self.lib
        if sp.issparse(counts):
            m = sp.diags(scale) @ counts
            m = m.tocsr(); m.data = np.log1p(m.data)
            return m
        return np.log1p(counts * scale[:, None])


def load_tissue(name, path=None):
    a = ad.read_h5ad(path or f"data/processed/{name}.h5ad")
    xy = np.c_[a.obs["x"].values, a.obs["y"].values].astype(np.float64)
    X = sp.csr_matrix(a.X)
    lib = np.asarray(X.sum(1)).ravel().astype(np.float64)
    return Tissue(name=name, xy=xy, types=a.obs["cell_type"].astype(str).values, counts=X,
                  genes=np.asarray(a.var_names), lib=lib)


def type_pair_contacts(t: Tissue, r: float, types=None):
    """Ordered type-pair contact counts within r and enrichment vs. random labelling."""
    types = types or t.type_list()
    idx = {k: i for i, k in enumerate(types)}
    code = np.array([idx.get(x, -1) for x in t.types])
    pairs = t.tree.query_pairs(r, output_type="ndarray")
    a, b = code[pairs[:, 0]], code[pairs[:, 1]]
    ok = (a >= 0) & (b >= 0)
    a, b = a[ok], b[ok]
    T = len(types)
    C = np.zeros((T, T))
    np.add.at(C, (a, b), 1)
    np.add.at(C, (b, a), 1)          # symmetric contact counts (ordered pair view)
    np.fill_diagonal(C, np.diag(C) / 2)
    f = np.array([(code == i).mean() for i in range(T)])
    total = ok.sum()
    expected = 2 * total * np.outer(f, f) / max((code >= 0).mean() ** 2, 1e-12)
    return types, C, C / np.maximum(expected, 1e-9)
