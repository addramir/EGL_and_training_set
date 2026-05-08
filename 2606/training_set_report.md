# L2G Training Set Construction — Release Comparison 25.06 vs 26.03

> Notebooks: `01_EGL_preparation.ipynb` · `02_training_set.ipynb`  
> Releases compared: **25.06** (`2506/`) and **26.03** (`2606/`)  
> OT data: `gs://open-targets-data-releases/{25.06,26.03}/`

---

## Key Parameter Changes Between Releases

| Parameter | 25.06 | 26.03 | Effect |
|-----------|-------|-------|--------|
| Rare variant evidence score threshold | **≥ 0.95** | **≥ 0.75** | +77% more EGL pairs from rare variants |
| STRING interaction score threshold | **≥ 0.80** | **≥ 0.75** | More negatives removed (378k vs 314k pairs) |
| EuroPMC literature evidence | **Skipped** | **Skipped** | No change (skipped in both) |
| ChEMBL filtering | `clinicalPhase` ∈ {3, 4} (numeric) | `clinicalStage` ∈ {PHASE\_3, PHASE\_4, APPROVAL, PREAPPROVAL} | Slightly broader coverage |
| Qualified studies pre-filter | **Yes** (67,751 studies) | **No** | 2606 retains many more loci |
| Replication list | Pre-computed (263,705 CSs) | Computed on-the-fly (≥ 2 GWAS per variant-disease) | 2606 filter is considerably **more lenient** — 654k vs 264k CSs pass |

> ⚠️ **Note for future iterations:** The replication filter in 26.03 is too permissive. Computing it on-the-fly with a simple count ≥ 2 threshold yields 654,390 passing credible sets vs 263,705 in the pre-computed 25.06 list. A stricter criterion — e.g. replication across independent cohorts or ancestries, not just any two GWAS entries — should be used in the next release.

---

## Stage 1 — Effector Gene List (EGL)

### 1a. Rare Variant Evidence

| Data source        | 25.06 (score ≥ 0.95) | 26.03 (score ≥ 0.75) | Δ        |
|--------------------|---------------------:|---------------------:|----------|
| EVA                |               57,746 |              360,991 | +525%    |
| Genomics England   |               31,037 |               42,455 | +37%     |
| UniProt Variants   |               32,368 |               36,120 | +12%     |
| Orphanet           |                6,289 |                7,234 | +15%     |
| UniProt Literature |                6,239 |                6,191 | −1%      |
| Gene2Phenotype     |                2,753 |                3,971 | +44%     |
| ClinGen            |                1,968 |                2,555 | +30%     |
| **Total rows**     |          **138,400** |          **459,517** | **+232%** |
| **Distinct pairs** |           **16,371** |           **29,049** | **+77%** |

The drop in score threshold from 0.95 → 0.75 has the largest impact on EVA, which increases more than 6-fold. EVA dominates the count at 0.75.

### 1b. ChEMBL Clinical Precedence

| Metric | 25.06 | 26.03 | Δ |
|--------|------:|------:|---|
| Filtering logic | `clinicalPhase` ∈ {3.0, 4.0} | `clinicalStage` ∈ {PHASE\_3, PHASE\_4, APPROVAL, PREAPPROVAL} | Broader |
| Distinct disease–target pairs | 25,685 | 25,782 | +97 |

Surprisingly, ChEMBL contributed only **97 net new pairs** between releases despite a broader stage definition. This was unexpected — the hypothesis going into 26.03 was that the updated clinical precedence pipeline would surface meaningfully more drug-target-disease links. The near-identical yield is due to a change in the OT platform's evidence pipeline: in 26.03 ChEMBL evidence is stored in a dedicated `evidence_clinical_precedence/` dataset with a restructured schema (using `clinicalStage` strings and a normalised `score` rather than `clinicalPhase` numerics). The restructuring changed which records are emitted and deduplicated, leaving the net distinct pairs almost unchanged. This is worth revisiting — the expectation was that finer-grained clinical stage annotation would unlock new EGL pairs.

### 1c. Legacy OTG Gold Standard

Identical in both releases (same source JSON file: `otg_gs_230511.json`).

| Confidence | Entries |
|-----------|--------:|
| High       |     629 |
| Medium     |     650 |
| **Distinct pairs after explode** | **812** |

### 1d. Combined EGL

| Source | 25.06 | 26.03 | Δ |
|--------|------:|------:|---|
| OTG Gold Standard | 812 | 812 | — |
| ChEMBL phase 3+ | 25,685 | 25,782 | +97 |
| Rare variants | 16,371 | 29,049 | +12,678 |
| **Combined (distinct)** | **42,288** | **55,003** | **+30%** |

The EGL grows by 30%, driven almost entirely by the lower rare variant score threshold.

---

## Stage 2 — Training Set Construction

### 2.0. Feature Matrix & EGL Join

| Item | 25.06 | 26.03 | Δ |
|------|------:|------:|---|
| Feature matrix total rows | 31,538,858 | 58,050,449 | +84% |
| EGL pairs | 42,288 | 55,003 | +30% |
| EGL-matched credible sets (GSP=1) | 21,799 | 27,364 | +26% |
| Unique gene–disease combos in GSP | 3,090 | 4,026 | +30% |
| FM rows for EGL study loci | 892,206 | 1,123,880 | +26% |
| Protein-coding GSP | 21,778 | 27,293 | +25% |

The 26.03 feature matrix is nearly **twice** the size of 25.06, reflecting a major expansion in credible sets between releases.

Saved: `unfiltered_FM_with_GSP`

---

### Filter 1 — Qualified Studies (25.06 only)

In 25.06 only, the FM was restricted to a curated set of "qualified" GWAS studies before applying replication.

| Item | 25.06 |
|------|------:|
| Qualified studies (GWAS + measurements) | 67,751 |
| FM rows after filter | 525,524 |
| GSP remaining | 12,520 |
| GSP lost | 9,279 (−43%) |

This step does not exist in 26.03, which contributes substantially to the larger downstream counts.

---

### Filter 2 — Replication Across Independent Studies

| Item | 25.06 | 26.03 |
|------|------:|------:|
| Method | Pre-computed replicated CS list | On-the-fly: variant × disease in ≥ 2 GWAS |
| Credible sets passing | **263,705** | **654,390** |
| FM rows after | 524,141 | 654,512 |
| GSP after | 12,883 | 15,999 |

The 26.03 approach is **2.5× more lenient**: 654k vs 264k credible sets pass. This is the main reason the 26.03 training set is larger — and the primary candidate for tightening in the next release.

---

### Filter 3 — Max 2 GSP Genes per Credible Set

| Item | 25.06 | 26.03 |
|------|------:|------:|
| CS with > 2 GSPs (removed) | 142 | 196 |
| CS with 1–2 GSPs (kept) | 12,011 | 14,711 |
| FM rows after | 517,925 | 645,718 |
| GSP after | 12,414 | 15,354 |

---

### Filter 4 — Protein–Protein Interaction (STRING) Filter

| Item | 25.06 | 26.03 |
|------|------:|------:|
| STRING score threshold | **≥ 0.80** | **≥ 0.75** |
| Interaction pairs retained | 314,292 | 378,436 |
| Distinct GSP genes | 394 | 517 |
| GSP genes with no STRING partners | 120 | 191 |
| Negative gene–locus rows flagged | 8,657 | 16,965 |
| FM rows after | 513,018 | 638,307 |
| GSP after | 12,414 *(unchanged)* | 15,354 *(unchanged)* |
| Unique GSP gene–disease pairs | 1,400 | 1,755 |

The lower STRING threshold in 26.03 (0.75 vs 0.80) includes more interaction partners, removing twice as many potential negative ambiguities (16,965 vs 8,657). This is the correct direction — the looser threshold is more conservative about what counts as a clean negative.

---

### Filter 5 — Distance Filter (remove footprint-zero GSPs)

| Item | 25.06 | 26.03 |
|------|------:|------:|
| GSP rows removed (distanceSentinelFootprint == 0) | 222 | 453 |
| FM rows after | 512,796 | 637,854 |
| GSP after | 12,192 | 14,901 |
| Non-GSP rows | 500,604 | 622,953 |
| Protein-coding GSP | 12,192 | 14,881 |
| Non-protein-coding GSP | 0 | 20 |
| Unique GSP gene–disease pairs | 1,377 | 1,713 |

---

### Filter 6 — Protein-Coding Genes Only + Re-apply Max-2-GSP

| Item | 25.06 | 26.03 |
|------|------:|------:|
| FM rows (protein-coding) | 189,240 | 236,144 |
| Study loci with 1–2 protein-coding GSPs | 11,834 | 14,343 |
| FM rows after max-2-GSP re-apply | 184,791 | 227,259 |

---

### Filter 7 — Deduplication Patch

| Item | 25.06 | 26.03 |
|------|------:|------:|
| Rows before | 184,791 | 227,259 |
| Reduction | **28.04%** | **11.76%** |
| Rows removed | 51,821 | 26,723 |

The 25.06 deduplication removed more than twice the fraction of rows. This suggests the 25.06 EGL contained more near-identical entries (likely due to the literature evidence generating repeated disease-gene pairs from many studies on the same locus). Excluding literature from 26.03 onward produces a cleaner EGL.

---

## Final Training Set Comparison

| Metric | 25.06 | 26.03 | Δ |
|--------|------:|------:|---|
| Total rows | **132,970** | **200,536** | +51% |
| Positive (GSP = 1) | 8,520 | 13,060 | +53% |
| Negative (GSP = 0) | 124,450 | 187,476 | +51% |
| Positive rate | 6.4% | 6.5% | ≈ same |
| Unique positive genes | 390 | 505 | +30% |
| Unique positive gene–disease pairs | 1,377 | 1,711 | +24% |

Output files:
- 25.06: `250624_training_set_full_fm.parquet` / `20250625_gentropy_paper_v1.json`
- 26.03: `2603_training_set_full_fm.parquet` / `2603_training_set_full_fm.json`

---

## Complete Filter Funnel

```
                                        25.06           26.03
Feature matrix (full)              31,538,858      58,050,449   (+84%)
                                         │               │
  Restrict to EGL study loci
Unfiltered FM + GSP labels            892,206       1,123,880   (+26%)
  GSP = 1                              21,799          27,364   (+26%)
                                         │               │
  [25.06 only] Qualified studies filter
After qualified studies               525,524               —
  GSP = 1                              12,520               —
                                         │               │
  Replication filter
After replication                     524,141         654,512   (+25%)
  GSP = 1                              12,883          15,999   (+24%)
  Passing credible sets               263,705         654,390  (+148%)
                                         │               │
  Max 2 GSP per credible set
After max-GSP filter                  517,925         645,718   (+25%)
  GSP = 1                              12,414          15,354   (+24%)
                                         │               │
  STRING interaction filter (0.80 / 0.75)
After STRING filter                   513,018         638,307   (+24%)
  GSP = 1                              12,414          15,354   (unchanged)
                                         │               │
  Distance filter (footprint ≠ 0)
After distance filter                 512,796         637,854   (+24%)
  GSP = 1                              12,192          14,901   (+22%)
                                         │               │
  Protein-coding only + max-2-GSP re-apply
After protein-coding filter           184,791         227,259   (+23%)
                                         │               │
  Deduplication (−28% / −12%)
FINAL TRAINING SET                    132,970         200,536   (+51%)
  Positive                              8,520          13,060   (+53%)
  Negative                            124,450         187,476   (+51%)
  Unique genes (positive)               390              505    (+29%)
  Unique gene–disease pairs           1,377            1,711    (+24%)
```

---

## Summary of Changes and Recommendations

| Change | Direction | Assessment |
|--------|-----------|-----------|
| Rare variant threshold 0.95 → 0.75 | More EGL pairs | Increases recall; trades some precision — worth monitoring if precision drops |
| STRING threshold 0.80 → 0.75 | More conservative negatives | Correct direction; removes more ambiguous negatives |
| Dropped qualified-studies filter | More loci pass | Increases volume but loses quality control on GWAS curation |
| Replication: pre-computed → on-the-fly | 2.5× more permissive | **Too lenient.** Next release should use stricter replication criteria (e.g., cross-cohort, cross-ancestry, or a curated study list) |
| FM size doubles (31M → 58M) | More data | Driven by the 26.03 release growth, not a pipeline choice |
| Literature evidence (EuroPMC) | Skipped in both | Consistent; literature adds noise via indirect associations |
