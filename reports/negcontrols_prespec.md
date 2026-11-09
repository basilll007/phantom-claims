# Pre-Specification: Anatomical Negative Controls on Real Data (Lane S, S4)

**Status:** Pre-registered before running any negative control evaluation on real tissue data.  
**Date:** 2026-11-08 23:57 CST  
**Branch:** `lane-S`  

---

## 1. Objective
Reviewers and domain experts frequently ask whether spatial CCC gates can avoid false positives on real, non-synthetic spatial transcriptomics data. A critical test is an **anatomical negative control**: testing whether the gate passes proposed interactions between cell types that occupy mutually exclusive anatomical compartments in real tissue sections (breast cancer, pancreatic ductal adenocarcinoma, lung cancer).

Because these cell types reside in segregated tissue compartments, direct juxtacrine/short-range paracrine signaling ($r \le 30\,\mu\text{m}$) should be physically absent. The expected pass rate for such claims is near 0.0.

---

## 2. Selection Rules for Anatomical Negative Type Pairs

An anatomical negative control claim $(A, B, L, R)$ must satisfy three strict rules:

1. **Mutually Exclusive Tissue Compartments (Spatial Segregation):**
   - **Breast Cancer (`breast_rep1`, `breast_rep2`):**
     - Ductal compartment (`DCIS_1`, `DCIS_2`) vs Distant Stroma / Lymphocyte compartment (`B_Cells`, `Mast_Cells`). In situ ductal carcinoma cells are sequestered within ducts, physically segregated from distal stromal leukocytes.
     - Myoepithelial compartment (`Myoepi_ACTA2+`, `Myoepi_KRT15+`) vs Distant Lymphocyte compartment (`B_Cells`).
   - **Pancreatic Cancer (`pancreas`):**
     - Normal acinar parenchyma (`Acinar`) vs invasive tumor compartment (`Epithelial_Tumor`).
     - Endocrine islet cells (`Endocrine`) vs `Epithelial_Tumor`.
   - **Lung Cancer (`lung`):**
     - Airway epithelium (`Ciliated`) vs invasive tumor compartment (`Epithelial_Tumor`).
     - Airway epithelium (`Ciliated`) vs distant plasma cells (`Plasma_Cells`).

2. **Genuine Gene Expression (Non-trivial Expression):**
   - Both ligand $L$ and receptor $R$ must be genuinely expressed on the panel:
     - Ligand $L$ expressed in sender type $A$: cell fraction with raw count $> 0$ must be $\ge 0.05$ (or type mean $\ge 0.02$).
     - Receptor $R$ expressed in receiver type $B$: cell fraction with raw count $> 0$ must be $\ge 0.05$ (or type mean $\ge 0.02$).
   - This rule guarantees that claims are not trivial "zero-expression" pairs, but reflect cases where standard non-spatial tools could propose an interaction based on expression alone.

3. **Curated Literature L-R Pairs:**
   - Pairs $(L, R)$ are drawn from the validated panel L-R sets (`data/processed/lr_pairs_<tissue>.csv`), derived from OmniPath / LIANA.

---

## 3. Exact Selection Code

```python
def select_anatomical_negatives(tissue_obj, lr_pairs_df, min_frac=0.05):
    # Compartment segregation mappings
    segregated_pairs = {
        "breast_rep1": [
            ("DCIS_1", "B_Cells"), ("DCIS_2", "B_Cells"),
            ("DCIS_1", "Mast_Cells"), ("DCIS_2", "Mast_Cells"),
            ("Myoepi_ACTA2+", "B_Cells"), ("Myoepi_KRT15+", "B_Cells"),
        ],
        "breast_rep2": [
            ("DCIS_1", "B_Cells"), ("DCIS_2", "B_Cells"),
            ("DCIS_1", "Mast_Cells"), ("DCIS_2", "Mast_Cells"),
            ("Myoepi_ACTA2+", "B_Cells"), ("Myoepi_KRT15+", "B_Cells"),
        ],
        "pancreas": [
            ("Acinar", "Epithelial_Tumor"), ("Epithelial_Tumor", "Acinar"),
            ("Endocrine", "Epithelial_Tumor"), ("Epithelial_Tumor", "Endocrine"),
            ("Acinar", "Plasma_Cells"), ("Plasma_Cells", "Acinar"),
        ],
        "lung": [
            ("Ciliated", "Epithelial_Tumor"), ("Epithelial_Tumor", "Ciliated"),
            ("Ciliated", "Plasma_Cells"), ("Plasma_Cells", "Ciliated"),
            ("Ciliated", "B_Cells"), ("B_Cells", "Ciliated"),
        ],
    }

    tname = tissue_obj.name
    type_pairs = segregated_pairs.get(tname, [])
    
    # Calculate per-type fraction expressing each gene
    genes = sorted(set(lr_pairs_df.ligand) | set(lr_pairs_df.receptor))
    tix = {g: i for i, g in enumerate(tissue_obj.genes)}
    raw = tissue_obj.counts
    
    # Filter to pairs where L is expressed in A and R is expressed in B
    valid_claims = []
    for sender, receiver in type_pairs:
        s_mask = (tissue_obj.types == sender)
        r_mask = (tissue_obj.types == receiver)
        if s_mask.sum() == 0 or r_mask.sum() == 0:
            continue
        for _, row in lr_pairs_df.iterrows():
            L, R = row.ligand, row.receptor
            if L not in tix or R not in tix:
                continue
            frac_L = float((raw[s_mask, tix[L]] > 0).mean())
            frac_R = float((raw[r_mask, tix[R]] > 0).mean())
            if frac_L >= min_frac and frac_R >= min_frac:
                valid_claims.append({
                    "tissue": tname,
                    "sender": sender,
                    "receiver": receiver,
                    "ligand": L,
                    "receptor": R,
                    "fracL": frac_L,
                    "fracR": frac_R
                })
    return valid_claims
```

---

## 4. Evaluation Protocol and Falsification Rules

1. **Frozen Gate Configuration:**
   $r = 30.0\,\mu\text{m}$, tile $= 250.0\,\mu\text{m}$, $n_{\text{perm}} = 1000$, $\alpha_{\text{half}} = 0.10$, $n_{\text{min}} = 50$, $\rho_{\text{min}} = 0.02$, $q = 0.05$.
   Response variable: receptor expression in receiver cells ($Y_i = \log(1 + R_i)$) as in Phase 1 real-data evaluation.
2. **Metrics:**
   - Pass rate $= (\text{passed claims}) / (\text{total selected anatomical negative claims})$ with 95% Clopper-Pearson or bootstrap confidence interval.
   - Expected outcome: pass rate $\le 0.05$ (ideally 0.0).
3. **Audit of False Passes:**
   - For every claim that passes the gate, we perform a root-cause audit:
     - Check $n_{\text{expressing\_pairs\_within\_r}}$.
     - Check tile spatial distribution and correlation $\rho$.
     - Diagnostic check: examine whether total cell density or library size (sum of transcripts per cell) accounts for the correlation.

