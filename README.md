# README
## Test Dataset  
The test dataset was constructed from proteins in UniProt that received Gene Ontology (GO) annotations after January 1, 2025, together with additional annotated human proteins to increase coverage.  

To ensure a fair and leakage-free evaluation, several filtering steps were applied:
- **Training overlap removal**: Proteins present in the training set were removed.
- **Homology filtering**: Test proteins were excluded if they had an MMseqs2 hit to any training protein with sequence identity ≥ 30% and alignment coverage ≥ 50%.
- **Label space alignment**: GO annotations were restricted to the label space observed in the training set.
- **Ontology propagation**: All GO annotations were propagated to their ancestor terms using a consistent GO ontology.

After processing, the final test set contains 1,363 proteins with 17,131 valid annotations, an average annotation density of 12.57, and a label space size of 2,062.  

The dataset covers all three GO aspects—Biological Process (BP), Molecular Function (MF), and Cellular Component (CC)—and exhibits a realistic long-tail distribution of GO terms, making it suitable for evaluating model performance under both common and rare functional labels.


## Baselines
### BLAST
``` bash
python baselines/BLAST/blast.py --base-fasta fasta/cafa/cafa.fasta --base-truth fasta/cafa/cafa.tsv --in fasta/test.fasta
```
### DeepGO-SE

### ESM2 + MLP 
```
# Training
python baselines/ESM2+MLP/esm2_mlp.py --mode train --base fasta/cafa/cafa --model-path baselines/ESM2+MLP/model.pt

# Inference
python baselines/ESM2+MLP/esm2_mlp.py --mode predict --in fasta/test --model-path baselines/ESM2+MLP/model.pt
```

### Naive
``` bash
python baselines/Naive/naive.py
``` 

## Evaluation

### Ontology propagation
Both ground-truth and predicted annotations were propagated along the Gene Ontology hierarchy (is_a and part_of relations) using the `go-basic.obo` ontology release. For predicted annotations, each ancestor term received the maximum confidence score among its descendants. Performance was assessed independently for Biological Process (BP), Molecular Function (MF), and Cellular Component (CC).

### Evaluation setting
We report results under an **open label-space** setting: predictions were evaluated against all GO terms present in the propagated ground truth, rather than being restricted to the terms observed in the benchmark annotations. This setting better reflects the open-world nature of protein function prediction, where valid annotations may extend beyond the labels explicitly present in a given benchmark.

### Metrics

Following the CAFA protein-centric evaluation protocol, we assessed performance using
four complementary metrics: (F_{\max}), (S_{\min}), protein-centric AUPR, and
class-centric average AUC.
(F_{\max}) was computed as the maximum F-score over 101 thresholds from 0.00 to 1.00,
where precision was averaged over proteins with at least one prediction at a given
threshold and recall was averaged over all benchmark proteins. (S_{\min}) was defined
as the minimum semantic distance derived from remaining uncertainty and
misinformation, both weighted by information accretion (IA); IA values were obtained
independently from training annotations via an external pickle file. In addition, we
reported protein-centric AUPR by integrating the protein-level precision-recall curve
over the same 101 thresholds, and class-centric average AUC by computing a ROC AUC
for each GO term and averaging over terms with both positive and negative benchmark
examples.
