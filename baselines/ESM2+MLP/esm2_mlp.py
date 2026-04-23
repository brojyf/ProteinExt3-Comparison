from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset
from transformers import AutoTokenizer, EsmModel
from tqdm import tqdm


def extract_entry_id(header: str) -> str:
    title = header[1:] if header.startswith(">") else header
    parts = title.split("|")
    if len(parts) >= 2:
        return parts[1]
    return title.split()[0]


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records = []
    header = None
    sequence_lines = []

    with path.open("r", encoding="utf-8") as file:
        for line in file:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(sequence_lines)))
                header = line
                sequence_lines = []
                continue

            sequence_lines.append(line)

    if header is not None:
        records.append((header, "".join(sequence_lines)))

    return records


def read_truth_by_pid(path: Path) -> tuple[dict[str, set[str]], dict[str, str]]:
    truth_by_pid: dict[str, set[str]] = defaultdict(set)
    go_to_aspect: dict[str, str] = {}

    with path.open("r", encoding="utf-8") as file:
        for raw_line in file:
            line = raw_line.strip()
            if not line:
                continue

            parts = [part.strip() for part in line.split("\t")]
            if len(parts) < 2:
                continue

            if parts[0] in {"AUTHOR", "MODEL", "KEYWORDS", "END"}:
                continue
            if parts[0] == "EntryID" and len(parts) >= 2 and parts[1] == "term":
                continue

            protein_id = parts[0]
            go_term = parts[1]
            if not protein_id or not go_term.startswith("GO:"):
                continue

            truth_by_pid[protein_id].add(go_term)
            if len(parts) >= 3 and parts[2].lower() in {"c", "f", "p"}:
                go_to_aspect[go_term] = parts[2].lower()

    if not truth_by_pid:
        raise ValueError(f"无法从 {path} 解析任何 GO 注释")

    return truth_by_pid, go_to_aspect


def build_io_paths(base_prefix: str, input_prefix: str, output_path: str) -> dict[str, Path]:
    return {
        "base_fasta": Path(f"{base_prefix}.fasta"),
        "base_truth": Path(f"{base_prefix}.tsv"),
        "input_fasta": Path(f"{input_prefix}.fasta"),
        "input_truth": Path(f"{input_prefix}.tsv"),
        "output_tsv": Path(output_path),
    }


class MLP(nn.Module):
    def __init__(self, input_dim: int, output_dim: int):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(input_dim, 1024),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(1024, 512),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(512, output_dim),
            nn.Sigmoid(),
        )

    def forward(self, x):
        return self.net(x)


def embed_sequences(
    records: list[tuple[str, str]],
    tokenizer: AutoTokenizer,
    esm_model: EsmModel,
    device: torch.device,
    batch_size: int,
) -> tuple[list[str], np.ndarray]:
    protein_ids = [extract_entry_id(header) for header, _ in records]
    sequences = [sequence for _, sequence in records]
    embeddings = []

    total_batches = (len(sequences) + batch_size - 1) // batch_size
    for start in tqdm(range(0, len(sequences), batch_size), total=total_batches, desc="embedding"):
        batch_sequences = sequences[start : start + batch_size]
        inputs = tokenizer(
            batch_sequences,
            return_tensors="pt",
            padding=True,
            truncation=True,
            max_length=1024,
        ).to(device)

        with torch.no_grad():
            outputs = esm_model(**inputs)

        # Mean pooling over residue embeddings.
        mask = inputs["attention_mask"].unsqueeze(-1).float()
        pooled = (outputs.last_hidden_state * mask).sum(1) / mask.sum(1)
        embeddings.append(pooled.cpu().numpy())

    return protein_ids, np.concatenate(embeddings, axis=0)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["train", "predict"], required=True, help="train: 训练并保存模型；predict: 加载模型并预测")
    parser.add_argument("--base", default="fasta/uniprot/250101", help="训练集前缀，会自动使用 .fasta 和 .tsv（train 模式需要）")
    parser.add_argument("--in", dest="input_prefix", default="fasta/test", help="测试集前缀，会自动使用 .fasta（predict 模式需要）")
    parser.add_argument("--out", default="baselines/ESM2+MLP/predictions.tsv")
    parser.add_argument("--model-path", default="baselines/ESM2+MLP/model.pt", help="模型保存/加载路径")
    parser.add_argument("--model-name", default="facebook/esm2_t6_8M_UR50D")
    parser.add_argument("--min-freq", type=int, default=5)
    parser.add_argument("--threshold", type=float, default=0.01)
    parser.add_argument("--epochs", type=int, default=10)
    parser.add_argument("--batch-size", type=int, default=16)
    parser.add_argument("--train-batch-size", type=int, default=128)
    return parser.parse_args()


def run_train(args):
    base_fasta = Path(f"{args.base}.fasta")
    base_truth = Path(f"{args.base}.tsv")
    model_path = Path(args.model_path)

    for path in [base_fasta, base_truth]:
        if not path.exists():
            raise FileNotFoundError(f"找不到文件: {path}")

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    tokenizer = AutoTokenizer.from_pretrained(args.model_name)
    esm_model = EsmModel.from_pretrained(args.model_name).to(device)
    esm_model.eval()

    base_records = read_fasta(base_fasta)
    truth_by_pid, go_to_aspect = read_truth_by_pid(base_truth)

    base_records = [record for record in base_records if extract_entry_id(record[0]) in truth_by_pid]
    base_protein_ids, base_embeddings = embed_sequences(
        base_records,
        tokenizer=tokenizer,
        esm_model=esm_model,
        device=device,
        batch_size=args.batch_size,
    )

    go_counts = Counter()
    for protein_id in base_protein_ids:
        for go_term in truth_by_pid.get(protein_id, set()):
            go_counts[go_term] += 1

    frequent_gos = sorted([go_term for go_term, count in go_counts.items() if count >= args.min_freq])
    if not frequent_gos:
        raise ValueError("没有任何 GO term 满足 min-freq 条件")

    go2idx = {go_term: index for index, go_term in enumerate(frequent_gos)}
    label_matrix = np.zeros((len(base_protein_ids), len(frequent_gos)), dtype=np.float32)

    for row_index, protein_id in enumerate(base_protein_ids):
        for go_term in truth_by_pid.get(protein_id, set()):
            go_index = go2idx.get(go_term)
            if go_index is not None:
                label_matrix[row_index, go_index] = 1.0

    model = MLP(base_embeddings.shape[1], len(frequent_gos)).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)
    criterion = nn.BCELoss()

    dataset = TensorDataset(
        torch.FloatTensor(base_embeddings),
        torch.FloatTensor(label_matrix),
    )
    loader = DataLoader(dataset, batch_size=args.train_batch_size, shuffle=True)

    for epoch in range(args.epochs):
        model.train()
        total_loss = 0.0
        with tqdm(loader, desc=f"epoch {epoch + 1}/{args.epochs}", unit="batch") as pbar:
            for xb, yb in pbar:
                xb = xb.to(device)
                yb = yb.to(device)
                optimizer.zero_grad()
                pred = model(xb)
                loss = criterion(pred, yb)
                loss.backward()
                optimizer.step()
                total_loss += loss.item()
                pbar.set_postfix(loss=f"{total_loss / max(pbar.n, 1):.6f}")

    model_path.parent.mkdir(parents=True, exist_ok=True)
    torch.save({
        "state_dict": model.state_dict(),
        "frequent_gos": frequent_gos,
        "go_to_aspect": {go: go_to_aspect[go] for go in frequent_gos if go in go_to_aspect},
        "embedding_dim": base_embeddings.shape[1],
        "model_name": args.model_name,
    }, model_path)

    print(f"base_fasta: {base_fasta}")
    print(f"base_truth: {base_truth}")
    print(f"model_path: {model_path}")
    print(f"device: {device}")
    print(f"model_name: {args.model_name}")
    print(f"base_proteins: {len(base_protein_ids)}")
    print(f"frequent_go_terms: {len(frequent_gos)}")


def run_predict(args):
    input_fasta = Path(f"{args.input_prefix}.fasta")
    output_tsv = Path(args.out)
    model_path = Path(args.model_path)

    for path in [input_fasta, model_path]:
        if not path.exists():
            raise FileNotFoundError(f"找不到文件: {path}")

    checkpoint = torch.load(model_path, map_location="cpu")
    frequent_gos = checkpoint["frequent_gos"]
    go_to_aspect = checkpoint.get("go_to_aspect", {})
    embedding_dim = checkpoint["embedding_dim"]
    esm_model_name = checkpoint.get("model_name", args.model_name)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    tokenizer = AutoTokenizer.from_pretrained(esm_model_name)
    esm_model = EsmModel.from_pretrained(esm_model_name).to(device)
    esm_model.eval()

    input_records = read_fasta(input_fasta)
    query_protein_ids, query_embeddings = embed_sequences(
        input_records,
        tokenizer=tokenizer,
        esm_model=esm_model,
        device=device,
        batch_size=args.batch_size,
    )

    model = MLP(embedding_dim, len(frequent_gos)).to(device)
    model.load_state_dict(checkpoint["state_dict"])
    model.eval()

    with torch.no_grad():
        predictions = model(torch.FloatTensor(query_embeddings).to(device)).cpu().numpy()

    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    written_rows = 0
    with output_tsv.open("w", encoding="utf-8") as file:
        for row_index, protein_id in enumerate(query_protein_ids):
            for col_index, go_term in enumerate(frequent_gos):
                score = float(predictions[row_index, col_index])
                if score >= args.threshold:
                    aspect = go_to_aspect.get(go_term, "")
                    file.write(f"{protein_id}\t{go_term}\t{score:.6f}\t{aspect}\n")
                    written_rows += 1

    print(f"input_fasta: {input_fasta}")
    print(f"output_tsv: {output_tsv}")
    print(f"model_path: {model_path}")
    print(f"device: {device}")
    print(f"query_proteins: {len(query_protein_ids)}")
    print(f"frequent_go_terms: {len(frequent_gos)}")
    print(f"prediction_rows: {written_rows}")


def main():
    args = parse_args()
    if args.mode == "train":
        run_train(args)
    else:
        run_predict(args)


if __name__ == "__main__":
    main()
