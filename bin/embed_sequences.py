#!/usr/bin/env python3
"""
Compute per-protein embeddings and save to embeddings.npz.

--model none     : amino-acid composition (mean-pooled one-hot), 21-dimensional
--model esm2     : ESM2 650M, mean-pooled over residues (dim 1280)
--model prot_t5  : ProtT5-XL, mean-pooled over residues (dim 1024)

Each entry in the NPZ file is keyed by protein ID and holds a 1-D float32 array.

CPU is a supported mode for the neural models, but --require-gpu turns "no usable
CUDA device" into an error instead of a silent fall back to it. processes/training.nf
passes that flag exactly when the run was launched with -profile gpu: having asked
for a GPU, quietly getting a ~60x slower CPU run is never what was meant, and at
real dataset sizes it walks straight into the scheduler's walltime.
"""

import argparse
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import read_fasta

_AA_VOCAB = "ACDEFGHIKLMNPQRSTVWY"
_AA_IDX = {aa: i for i, aa in enumerate(_AA_VOCAB)}

# Both neural models truncate here, which also caps a batch's padded width.
_MAX_LEN = 1024


def embed_one_hot(seqs):
    dim = len(_AA_VOCAB) + 1  # 20 standard AAs + 1 for unknown/non-standard
    embeddings = {}
    for pid, seq in seqs.items():
        vec = np.zeros(dim, dtype=np.float32)
        for aa in seq:
            vec[_AA_IDX.get(aa, len(_AA_VOCAB))] += 1
        if seq:
            vec /= len(seq)
        embeddings[pid] = vec
    return embeddings


def resolve_device(require_gpu):
    """Pick cuda when it is genuinely usable, else error or warn per require_gpu."""
    import torch

    built_for = torch.version.cuda or "cpu-only"
    if torch.cuda.is_available():
        name = torch.cuda.get_device_name(0)
        print(f"  device: cuda ({name}); torch {torch.__version__} built for CUDA {built_for}", file=sys.stderr)
        return torch.device("cuda")

    detail = f"torch {torch.__version__} (built for CUDA {built_for}) reports no usable CUDA device"
    if require_gpu:
        sys.exit(
            f"ERROR: {detail}, but this run was launched with -profile gpu.\n"
            "  A GPU may well have been allocated -- torch refusing one usually means the\n"
            "  installed wheel targets a newer CUDA than the node's driver supports. Look in\n"
            "  .command.err for a 'driver is too old' warning and install a matching cu12x\n"
            "  build (see environment.yml). Embedding on CPU takes roughly 2.5 s per sequence,\n"
            "  which at real dataset sizes exceeds any reasonable walltime, so this is an\n"
            "  error rather than a fall back. To embed on CPU deliberately, drop -profile gpu."
        )
    print(f"WARNING: {detail}; embedding on CPU at roughly 2.5 s per sequence.", file=sys.stderr)
    return torch.device("cpu")


def batches(tokens, max_tokens, max_batch):
    """Group ids into batches under a padded-token budget, longest sequence first.

    Sorting by length keeps a batch's members close in size, so padding to the
    longest one wastes little; the budget counts PADDED tokens because that is
    what actually reaches the GPU. The id tie-break keeps batching deterministic:
    dict order here follows --fasta, which is a groupTuple()'s completion order.
    """
    batch = []
    for pid in sorted(tokens, key=lambda p: (-tokens[p], p)):
        width = tokens[batch[0]] if batch else tokens[pid]
        if batch and ((len(batch) + 1) * width > max_tokens or len(batch) >= max_batch):
            yield batch
            batch = []
        batch.append(pid)
    if batch:
        yield batch


def residue_mask(attention_mask):
    """attention_mask with the CLS (first) and EOS (last real) positions cleared.

    The same two tokens the one-at-a-time version dropped with [0, 1:-1]. Taking
    EOS from the mask rather than from index -1 is what keeps that slice correct
    once rows are padded to a common width.
    """
    import torch

    mask = attention_mask.clone()
    mask[:, 0] = 0
    eos = attention_mask.sum(1) - 1
    mask[torch.arange(mask.size(0), device=mask.device), eos] = 0
    return mask


def mean_pool(hidden, mask):
    """Mean over the positions mask keeps.

    clamp() covers a sequence with no positions left -- an empty FASTA entry,
    whose tokens are CLS and EOS and nothing else. That used to divide by zero
    and store NaN; it now stores a zero vector.
    """
    m = mask.unsqueeze(-1).to(hidden.dtype)
    return (hidden * m).sum(1) / m.sum(1).clamp(min=1)


def _report(done, total, reported, step=1000):
    if done - reported >= step or done == total:
        print(f"  {done}/{total} embedded", file=sys.stderr)
        return done
    return reported


def embed_esm2(seqs, device, max_tokens, max_batch):
    import torch
    from transformers import AutoModel, AutoTokenizer

    tokenizer = AutoTokenizer.from_pretrained("facebook/esm2_t33_650M_UR50D")
    model = AutoModel.from_pretrained("facebook/esm2_t33_650M_UR50D").to(device).eval()

    # +2 for the CLS and EOS the tokenizer wraps each sequence in.
    tokens = {pid: min(len(seq), _MAX_LEN - 2) + 2 for pid, seq in seqs.items()}
    embeddings = {}
    done = reported = 0
    for batch in batches(tokens, max_tokens, max_batch):
        enc = tokenizer(
            [seqs[pid] for pid in batch],
            return_tensors="pt",
            padding=True,
            truncation=True,
            max_length=_MAX_LEN,
        ).to(device)
        with torch.no_grad():
            hidden = model(**enc).last_hidden_state
        pooled = mean_pool(hidden, residue_mask(enc["attention_mask"]))
        for pid, vec in zip(batch, pooled):
            embeddings[pid] = vec.cpu().float().numpy()
        done += len(batch)
        reported = _report(done, len(seqs), reported)
    return embeddings


def embed_prot_t5(seqs, device, max_tokens, max_batch):
    import torch
    from transformers import T5EncoderModel, T5Tokenizer

    tokenizer = T5Tokenizer.from_pretrained("Rostlab/prot_t5_xl_half_uniref50-enc", do_lower_case=False)
    model = T5EncoderModel.from_pretrained("Rostlab/prot_t5_xl_half_uniref50-enc").to(device).eval()

    # ProtT5 expects space-separated AAs; replace non-standard residues with X
    fmt = {pid: " ".join(re.sub(r"[UZOB]", "X", seq)) for pid, seq in seqs.items()}
    # +1 for the trailing </s>; T5 prepends no CLS.
    tokens = {pid: min(len(seq), _MAX_LEN - 1) + 1 for pid, seq in seqs.items()}
    embeddings = {}
    done = reported = 0
    for batch in batches(tokens, max_tokens, max_batch):
        enc = tokenizer(
            [fmt[pid] for pid in batch],
            return_tensors="pt",
            padding=True,
            truncation=True,
            max_length=_MAX_LEN,
        ).to(device)
        with torch.no_grad():
            hidden = model(**enc).last_hidden_state
        # Pools every real position, the trailing </s> included -- which is what
        # last_hidden_state[0].mean(0) did one sequence at a time. Dropping it
        # here would be defensible but would shift every ProtT5 embedding.
        pooled = mean_pool(hidden, enc["attention_mask"])
        for pid, vec in zip(batch, pooled):
            embeddings[pid] = vec.cpu().float().numpy()
        done += len(batch)
        reported = _report(done, len(seqs), reported)
    return embeddings


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", nargs="+", required=True)
    ap.add_argument("--model", choices=["none", "esm2", "prot_t5"], default="esm2")
    ap.add_argument(
        "--require-gpu",
        action="store_true",
        help="abort instead of falling back to CPU when no CUDA device is usable",
    )
    ap.add_argument("--batch-tokens", type=int, default=16384, help="padded tokens per forward pass")
    ap.add_argument("--max-batch", type=int, default=64, help="sequences per forward pass")
    args = ap.parse_args()

    seqs = {}
    for path in args.fasta:
        seqs.update(read_fasta(path))
    print(f"Embedding {len(seqs)} unique proteins with {args.model}", file=sys.stderr)

    if args.model == "none":
        # Pure numpy, so a GPU is neither used nor asked about even under -profile gpu.
        embeddings = embed_one_hot(seqs)
    else:
        device = resolve_device(args.require_gpu)
        embed = embed_esm2 if args.model == "esm2" else embed_prot_t5
        embeddings = embed(seqs, device, args.batch_tokens, args.max_batch)

    # Sorted, so the archive is byte-reproducible. np.savez writes members in the
    # order given and zeroes their timestamps, so the only thing that made two
    # runs differ was insertion order -- which follows --fasta, and that comes
    # from a groupTuple() whose order is task-completion order. Contents were
    # always identical; this stops a no-regression `diff -r` on results/ from
    # reporting the shared .npz as changed when nothing has.
    np.savez("embeddings.npz", **{pid: embeddings[pid] for pid in sorted(embeddings)})
    print(f"Saved {len(embeddings)} embeddings to embeddings.npz", file=sys.stderr)


if __name__ == "__main__":
    main()
