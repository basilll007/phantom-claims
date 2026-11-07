"""Expression-only CCC baselines in the CellPhoneDB / LIANA-CellChat spirit (no spatial information).
expr_cpdb: score = (mean_L(A) + mean_R(B)) / 2 ; expr_prod: score = mean_L(A) * mean_R(B).
Both: claim valid iff >= min_frac of A express L and of B express R; p-value from cell-label permutations."""
import numpy as np
import pandas as pd
import torch


def _type_means(X, code, T):
    M = torch.zeros((T, X.shape[1]), device=X.device)
    M.index_add_(0, code, X)
    cnt = torch.bincount(code, minlength=T).clamp(min=1).unsqueeze(1)
    return M / cnt


def score_all(t, expr_lognorm, expr_raw, genes, pairs, n_perm=1000, seed=0, device=None):
    """t: Tissue; expr_*: (n, G) dense arrays for `genes`; pairs: list of (ligand, receptor).
    Returns DataFrame over all (A, B, L, R), A != B."""
    device = device or ("cuda" if torch.cuda.is_available() else "cpu")
    types = t.type_list()
    tix = {k: i for i, k in enumerate(types)}
    code_np = np.array([tix.get(x, -1) for x in t.types])
    keep = code_np >= 0
    T = len(types)
    gi = {g: i for i, g in enumerate(genes)}
    Li = torch.tensor([gi[l] for l, _ in pairs], device=device)
    Ri = torch.tensor([gi[r] for _, r in pairs], device=device)
    X = torch.tensor(np.asarray(expr_lognorm)[keep], dtype=torch.float32, device=device)
    B = torch.tensor((np.asarray(expr_raw)[keep] > 0).astype(np.float32), device=device)
    code = torch.tensor(code_np[keep], device=device)
    M = _type_means(X, code, T)
    F = _type_means(B, code, T)

    def scores(M):
        mL = M[:, Li].T            # (P, T) ligand mean in sender type
        mR = M[:, Ri].T            # (P, T) receptor mean in receiver type
        cp = (mL[:, :, None] + mR[:, None, :]) / 2
        pr = mL[:, :, None] * mR[:, None, :]
        return cp, pr

    cp0, pr0 = scores(M)
    ge_cp = torch.zeros_like(cp0); ge_pr = torch.zeros_like(pr0)
    g = torch.Generator(device=device); g.manual_seed(int(seed))
    for _ in range(n_perm):
        perm = code[torch.randperm(code.shape[0], generator=g, device=device)]
        cp, pr = scores(_type_means(X, perm, T))
        ge_cp += (cp >= cp0 - 1e-12).float(); ge_pr += (pr >= pr0 - 1e-12).float()
    p_cp = ((1 + ge_cp) / (1 + n_perm)).cpu().numpy()
    p_pr = ((1 + ge_pr) / (1 + n_perm)).cpu().numpy()
    cp0, pr0 = cp0.cpu().numpy(), pr0.cpu().numpy()
    Mn, Fn = M.cpu().numpy(), F.cpu().numpy()
    rows = []
    for k, (l, r) in enumerate(pairs):
        for a in range(T):
            for b in range(T):
                if a == b:
                    continue
                rows.append((types[a], types[b], l, r, cp0[k, a, b], p_cp[k, a, b], pr0[k, a, b], p_pr[k, a, b],
                             Mn[a, gi[l]], Mn[b, gi[r]], Fn[a, gi[l]], Fn[b, gi[r]]))
    return pd.DataFrame(rows, columns=["sender", "receiver", "ligand", "receptor", "score_cpdb", "p_cpdb",
                                       "score_prod", "p_prod", "meanL", "meanR", "fracL", "fracR"])


def top_k(df, method, K, min_frac=0.1, alpha=0.05):
    s, p = f"score_{method}", f"p_{method}"
    d = df[(df.fracL >= min_frac) & (df.fracR >= min_frac) & (df[p] < alpha)]
    d = d.sort_values([s, "sender", "receiver", "ligand"], ascending=[False, True, True, True], kind="mergesort")
    return d.head(K)


def run_on_bench(t, inst, n_perm=1000, seed=0):
    C = inst.counts
    Xl = t.lognorm(C).toarray()
    Xr = C.toarray()
    pairs = [(d["ligand"], d["receptor"]) for d in inst.database]
    return score_all(t, Xl, Xr, inst.syn_genes, pairs, n_perm=n_perm, seed=seed)
