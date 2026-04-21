# README

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

# Evaluation
``` bash
python comparison/compare.py
```
