# Result

## Dataset
### Dataset 1
Using proteins from the UniProtKB/Swiss-Prot dataset released on or after January 1, 2025, then removing proteins that appear in the CAFA6 training set to form `baselines/fasta/test.fasta` and `baselines/fasta/test.tsv`. There are 1571 proteins in the full test set. Ontology-specific benchmark sizes are MF 1256, CC 1112, and BP 969.

### Dataset 2
Using human proteins derived from `baselines/fasta/human/human_filtered.fasta` and `baselines/fasta/human/human_filtered.tsv`, then removing sequences that also appear in the CAFA6 training FASTA to form `baselines/fasta/test_human.fasta` and `baselines/fasta/test_human.tsv`. There are 2549 proteins in the full human test set. Ontology-specific benchmark sizes are MF 1639, CC 2336, and BP 1741.

## Comparison
### MF (Dataset 1)
| Method | Fmax | Smin | AUPR | AUC |
| --- | ---: | ---: | ---: | ---: |
| BLAST | 0.522 | 13.847 | 0.494 | 0.807 |
| DeepGO-SE | 0.546 | 14.387 | 0.543 | 0.827 |
| ESM2+MLP | 0.346 | 16.049 | 0.334 | 0.723 |
| Naive | 0.147 | 19.128 | 0.053 | 0.500 |
| ProteinExt3 | 0.385 | 15.198 | 0.315 | 0.648 |
| ProteinExt3 (ESM2) | 0.496 | **13.283** | 0.406 | 0.709 |
| ProteinExt3 (ProtT5) | 0.513 | 13.437 | 0.424 | 0.725 |
| ProteinExt3 (CNN) | 0.097 | 19.289 | 0.064 | 0.510 |
| SPROF-GO | **0.572** | 13.863 | **0.627** | **0.852** |

### MF (Dataset 2)

| Method | Fmax | Smin | AUPR | AUC |
| --- | ---: | ---: | ---: | ---: |
| BLAST | 0.566 | 10.502 | 0.532 | 0.838 |
| DeepGO-SE | 0.625 | 10.876 | 0.539 | 0.838 |
| ESM2+MLP | 0.463 | 12.071 | 0.384 | 0.762 |
| Naive | 0.186 | 14.569 | 0.048 | 0.500 |
| ProteinExt3 | 0.424 | 11.810 | 0.265 | 0.619 |
| SPROF-GO | **0.677** | **9.940** | **0.645** | **0.866** |

### CC (Dataset 1)

| Method | Fmax | Smin | AUPR | AUC |
| --- | ---: | ---: | ---: | ---: |
| BLAST | 0.439 | 8.716 | 0.252 | 0.726 |
| DeepGO-SE | 0.636 | 8.639 | 0.344 | 0.758 |
| ESM2+MLP | 0.631 | 6.839 | 0.306 | 0.751 |
| Naive | 0.378 | 8.512 | 0.067 | 0.500 |
| ProteinExt3 | **0.661** | **6.206** | 0.279 | 0.654 |
| ProteinExt3 (ESM2) | 0.638 | 6.426 | 0.267 | 0.650 |
| ProteinExt3 (ProtT5) | 0.639 | 6.527 | 0.273 | 0.657 |
| ProteinExt3 (CNN) | 0.308 | 8.741 | 0.076 | 0.508 |
| SPROF-GO | 0.650 | 8.517 | **0.379** | **0.778** |

### CC (Dataset 2)

| Method | Fmax | Smin | AUPR | AUC |
| --- | ---: | ---: | ---: | ---: |
| BLAST | 0.603 | 5.792 | **0.396** | **0.853** |
| DeepGO-SE | 0.625 | 8.008 | 0.223 | 0.714 |
| ESM2+MLP | 0.652 | 6.402 | 0.202 | 0.748 |
| Naive | 0.456 | 7.691 | 0.047 | 0.500 |
| ProteinExt3 | **0.690** | **5.694** | 0.252 | 0.644 |
| SPROF-GO | 0.630 | 7.982 | 0.260 | 0.714 | 

### BP (Dataset 1)

| Method | Fmax | Smin | AUPR | AUC |
| --- | ---: | ---: | ---: | ---: |
| BLAST | **0.443** | 20.389 | 0.310 | 0.791 |
| DeepGO-SE | 0.419 | 20.335 | **0.365** | **0.798** |
| ESM2+MLP | 0.346 | 18.282 | 0.242 | 0.724 |
| Naive | 0.095 | 20.421 | 0.044 | 0.500 |
| ProteinExt3 | 0.327 | **17.387** | 0.197 | 0.609 |
| ProteinExt3 (ESM2) | 0.310 | 17.756 | 0.187 | 0.605 |
| ProteinExt3 (ProtT5) | 0.342 | 17.677 | 0.200 | 0.618 |
| ProteinExt3 (CNN) | 0.023 | 20.338 | 0.048 | 0.503 |
| SPROF-GO | 0.434 | 19.857 | 0.252 | 0.642 |

### BP (Dataset 2)

| Method | Fmax | Smin | AUPR | AUC |
| --- | ---: | ---: | ---: | ---: |
| BLAST | **0.524** | **13.101** | **0.442** | **0.881** |
| DeepGO-SE | 0.420 | 18.173 | 0.244 | 0.732 |
| ESM2+MLP | 0.373 | 16.288 | 0.191 | 0.682 |
| Naive | 0.178 | 17.778 | 0.035 | 0.500 |
| ProteinExt3 | 0.264 | 16.002 | 0.131 | 0.565 |
| SPROF-GO | 0.446 | 17.855 | 0.188 | 0.614 |


## Methods
- **BLAST**: Sequence similarity transfer baseline based on BLAST hits.
- **DeepGO-SE**: DeepGO-SE predictions filtered to the evaluated test set.
- **Naive**: GO term frequency-based prediction using CAFA6 training set statistics.
- **ESM2+MLP**: ESM2-650M embeddings + MLP classifier trained on CAFA6 training set.
- **ProteinExt3**: Full fusion model.
- **ProteinExt3 (ESM2)**: ESM2-only sub method.
- **ProteinExt3 (ProtT5)**: ProtT5-only sub method.
- **ProteinExt3 (CNN)**: CNN-only sub method.
- **SPROF-GO**: Sequence-based protein function predictor.

## Appendix

### Dataset 1 Thresholds

| Category | Method | Fmax Threshold | Smin Threshold |
| --- | --- | ---: | ---: |
| MF | BLAST | 0.279 | 0.963 |
| MF | DeepGO-SE | 0.268 | 0.509 |
| MF | ESM2+MLP | 0.018 | 0.040 |
| MF | Naive | 0.010 | 0.409 |
| MF | ProteinExt3 | 0.301 | 0.303 |
| MF | ProteinExt3 (ESM2) | 0.301 | 0.300 |
| MF | ProteinExt3 (ProtT5) | 0.300 | 0.300 |
| MF | ProteinExt3 (CNN) | 0.300 | 0.784 |
| MF | SPROF-GO | 0.157 | 0.323 |
| CC | BLAST | 0.816 | 0.999 |
| CC | DeepGO-SE | 0.589 | 0.945 |
| CC | ESM2+MLP | 0.148 | 0.175 |
| CC | Naive | 0.161 | 0.161 |
| CC | ProteinExt3 | 0.352 | 0.389 |
| CC | ProteinExt3 (ESM2) | 0.303 | 0.432 |
| CC | ProteinExt3 (ProtT5) | 0.342 | 0.556 |
| CC | ProteinExt3 (CNN) | 0.304 | 0.627 |
| CC | SPROF-GO | 0.393 | 0.715 |
| BP | BLAST | 0.946 | 1.000 |
| BP | DeepGO-SE | 0.473 | 0.967 |
| BP | ESM2+MLP | 0.021 | 0.075 |
| BP | Naive | 0.028 | N/A |
| BP | ProteinExt3 | 0.300 | 0.300 |
| BP | ProteinExt3 (ESM2) | 0.304 | 0.310 |
| BP | ProteinExt3 (ProtT5) | 0.300 | 0.337 |
| BP | ProteinExt3 (CNN) | 0.300 | 0.300 |
| BP | SPROF-GO | 0.319 | 0.848 |

### Dataset 2 Thresholds

| Category | Method | Fmax Threshold | Smin Threshold |
| --- | --- | ---: | ---: |
| MF | BLAST | 0.219 | 0.953 |
| MF | DeepGO-SE | 0.193 | 0.537 |
| MF | ESM2+MLP | 0.019 | 0.036 |
| MF | Naive | 0.409 | 0.409 |
| MF | ProteinExt3 | 0.306 | 0.384 |
| MF | SPROF-GO | 0.164 | 0.365 |
| CC | BLAST | 0.675 | 1.000 |
| CC | DeepGO-SE | 0.728 | 0.983 |
| CC | ESM2+MLP | 0.122 | 0.158 |
| CC | Naive | 0.123 | 0.161 |
| CC | ProteinExt3 | 0.300 | 0.367 |
| CC | SPROF-GO | 0.642 | 0.953 |
| BP | BLAST | 0.613 | 0.999 |
| BP | DeepGO-SE | 0.418 | 0.879 |
| BP | ESM2+MLP | 0.026 | 0.037 |
| BP | Naive | 0.028 | 0.028 |
| BP | ProteinExt3 | 0.301 | 0.301 |
| BP | SPROF-GO | 0.331 | 0.798 |
