#!/usr/bin/env python3
"""
Train a RandomForest PPI classifier with limited hyperparameter tuning.

Feature construction: concatenation of the two sorted protein embeddings.
Sorting by protein ID ensures the feature vector is the same regardless of
the order in which the pair appears in the CSV.

Hyperparameter search: 3 pre-defined configs (max_depth 5/10/30, max_samples 0.2)
evaluated by AUROC on val set. The best config is retrained on train+val, then
evaluated on both test sets.

Node-type agnostic: in DDI mode the rows are domain-instance pairs rather than
protein pairs and the embeddings are keyed by instance id, which changes nothing
here. What DDI mode adds is a second scoring level -- when the labelled CSV
carries family1/family2, several rows represent one DDI, so each DDI's example
predictions are averaged and scored again at family level in its own MultiQC
table. Example level asks "can it call this instance pair", family level "can it
call this DDI"; the latter is what makes the numbers comparable with the
downstream benchmark.
"""

import argparse
import os
import sys

import numpy as np
from sklearn.ensemble import RandomForestClassifier

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import load_embeddings, read_family_pairs, read_labelled_csv
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    f1_score,
    matthews_corrcoef,
    precision_score,
    recall_score,
    roc_auc_score,
)

HP_CONFIGS = [
    {"n_estimators": 100, "max_depth": 5, "max_samples": 0.2},
    {"n_estimators": 100, "max_depth": 10, "max_samples": 0.2},
    {"n_estimators": 100, "max_depth": 30, "max_samples": 0.2},
]


def build_X(pairs, labels, embeddings):
    """Return (X, y, kept), dropping pairs where either protein is missing.

    `kept` are the indices of the surviving pairs, so per-row metadata read
    separately from the same CSV -- the DDI family pairs -- can be filtered in
    lockstep instead of silently misaligning with y once a row is dropped.
    """
    rows, y, kept = [], [], []
    for i, ((p1, p2), label) in enumerate(zip(pairs, labels)):
        a, b = (p1, p2) if p1 <= p2 else (p2, p1)
        if a in embeddings and b in embeddings:
            rows.append(np.concatenate([embeddings[a], embeddings[b]]))
            y.append(label)
            kept.append(i)
    skipped = len(pairs) - len(kept)
    if skipped:
        print(f"  skipped {skipped} pairs with missing embeddings", file=sys.stderr)
    return np.array(rows), np.array(y), kept


def compute_metrics(y_true, y_prob):
    y_pred = (y_prob >= 0.5).astype(int)
    return {
        "auroc": roc_auc_score(y_true, y_prob),
        "auprc": average_precision_score(y_true, y_prob),
        "f1": f1_score(y_true, y_pred, zero_division=0),
        "mcc": matthews_corrcoef(y_true, y_pred),
        "precision": precision_score(y_true, y_pred, zero_division=0),
        "recall": recall_score(y_true, y_pred, zero_division=0),
        "accuracy": accuracy_score(y_true, y_pred),
    }


def aggregate_families(y_true, y_prob, fam_pairs):
    """Mean example probability per DDI -- (y, prob, both_labels) at family level.

    Keyed on (family pair, label) rather than the family pair alone: a pair that
    somehow carries both labels then stays two rows and is reported, instead of
    silently collapsing a positive and a negative DDI into one averaged score.
    EXPAND_NEGATIVES samples negative family pairs from outside the positive set,
    so `both_labels` is expected to be empty.
    """
    groups = {}
    for label, prob, fam in zip(y_true, y_prob, fam_pairs):
        groups.setdefault((fam, int(label)), []).append(float(prob))
    keys = sorted(groups)
    y = np.array([label for _, label in keys])
    prob = np.array([float(np.mean(groups[k])) for k in keys])
    both = {fam for fam, _ in keys if (fam, 0) in groups and (fam, 1) in groups}
    return y, prob, both


def write_mqc(results, id_):
    """One table per test split (test_balanced / test_realistic), each
    merging across datasets -- keeps either table from getting cluttered
    with the other's rows, matching the split evaluation of the two test
    sets."""
    cols = ["auroc", "auprc", "f1", "mcc", "precision", "recall", "accuracy"]
    for name, metrics in results:
        table_id = f"classifier_metrics_{name}"
        with open(f"{table_id}_mqc.tsv", "w") as fh:
            fh.write(
                f"# id: '{table_id}'\n"
                f"# section_name: 'Classifier Performance ({name})'\n"
                "# description: 'RandomForest PPI classifier. Hyperparameters tuned on val AUROC (3 configs: max_depth 5/10/30, max_samples 0.2), then retrained on train+val.'\n"
                "# plot_type: 'table'\n"
                "# pconfig:\n"
                f"#     id: '{table_id}_table'\n"
                f"#     title: 'RF Classifier - {name} Test Performance'\n"
                "# headers:\n"
                "#     ID:         {description: 'Dataset ID'}\n"
                "#     AUROC:      {format: '{:.4f}'}\n"
                "#     AUPRC:      {format: '{:.4f}'}\n"
                "#     F1:         {format: '{:.4f}'}\n"
                "#     MCC:        {format: '{:.4f}'}\n"
                "#     Precision:  {format: '{:.4f}'}\n"
                "#     Recall:     {format: '{:.4f}'}\n"
                "#     Accuracy:   {format: '{:.4f}'}\n"
                "Sample\tID\tAUROC\tAUPRC\tF1\tMCC\tPrecision\tRecall\tAccuracy\n"
            )
            # metrics is None for an empty or single-class split -- keep the row so the
            # dataset still appears in the table, with blank cells.
            row = "\t".join("" if metrics is None else f"{metrics[c]:.4f}" for c in cols)
            fh.write(f"{id_}\t{id_}\t{row}\n")


def write_family_mqc(results, id_):
    """DDI mode's second scoring level: one table per test split, family level.

    A separate file rather than extra columns on the example-level table, so the
    PPI-mode output above is byte-identical whether or not this ever runs.
    """
    cols = ["auroc", "auprc", "f1", "mcc", "precision", "recall", "accuracy"]
    for name, metrics, n_ddis in results:
        table_id = f"classifier_metrics_{name}_family"
        with open(f"{table_id}_mqc.tsv", "w") as fh:
            fh.write(
                f"# id: '{table_id}'\n"
                f"# section_name: 'Classifier Performance, DDI level ({name})'\n"
                # No apostrophes anywhere in this description: MultiQC parses these
                # header lines as YAML, and an unescaped ' inside a single-quoted
                # scalar ends it early -- which fails the whole MULTIQC task, not
                # just this section.
                "# description: 'The same predictions as the example-level table for this split, averaged "
                "over the domain-instance examples of each DDI and scored per Pfam family pair. Example "
                "level asks whether an instance pair can be called, family level whether the interaction "
                "can.'\n"
                "# plot_type: 'table'\n"
                "# pconfig:\n"
                f"#     id: '{table_id}_table'\n"
                f"#     title: 'RF Classifier - {name} Test Performance (per DDI)'\n"
                "# headers:\n"
                "#     ID:         {description: 'Dataset ID'}\n"
                "#     DDIs:       {description: 'Family pairs scored', format: '{:,.0f}'}\n"
                "#     AUROC:      {format: '{:.4f}'}\n"
                "#     AUPRC:      {format: '{:.4f}'}\n"
                "#     F1:         {format: '{:.4f}'}\n"
                "#     MCC:        {format: '{:.4f}'}\n"
                "#     Precision:  {format: '{:.4f}'}\n"
                "#     Recall:     {format: '{:.4f}'}\n"
                "#     Accuracy:   {format: '{:.4f}'}\n"
                "Sample\tID\tDDIs\tAUROC\tAUPRC\tF1\tMCC\tPrecision\tRecall\tAccuracy\n"
            )
            row = "\t".join("" if metrics is None else f"{metrics[c]:.4f}" for c in cols)
            fh.write(f"{id_}\t{id_}\t{n_ddis}\t{row}\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--train", required=True)
    ap.add_argument("--val", required=True)
    ap.add_argument("--test_balanced", required=True)
    ap.add_argument("--test_realistic", required=True)
    ap.add_argument("--embeddings", required=True)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--id", required=True, help="Dataset ID, for MultiQC tagging")
    args = ap.parse_args()

    print("Loading embeddings ...", file=sys.stderr)
    embeddings = load_embeddings(args.embeddings)
    print(f"  {len(embeddings)} proteins", file=sys.stderr)

    train_pairs, y_train = read_labelled_csv(args.train)
    val_pairs, y_val = read_labelled_csv(args.val)

    X_train, y_train, _ = build_X(train_pairs, y_train, embeddings)
    X_val, y_val, _ = build_X(val_pairs, y_val, embeddings)

    # Hyperparameter search on val AUROC. A split can legitimately be empty -- a 0
    # split fraction, or every pair dropped for a missing embedding -- so fall back to
    # the first config instead of failing: roc_auc_score raises on an empty y_true.
    if len(X_val) == 0:
        print("Val split is empty -- skipping tuning, using the first config", file=sys.stderr)
        best_cfg, best_auroc = HP_CONFIGS[0], float("nan")
    else:
        print("Tuning hyperparameters ...", file=sys.stderr)
        best_auroc, best_cfg = -1.0, None
        for cfg in HP_CONFIGS:
            clf = RandomForestClassifier(**cfg, random_state=args.seed, n_jobs=-1)
            clf.fit(X_train, y_train)
            auroc = roc_auc_score(y_val, clf.predict_proba(X_val)[:, 1])
            print(f"  {cfg}  →  val AUROC {auroc:.4f}", file=sys.stderr)
            if auroc > best_auroc:
                best_auroc, best_cfg = auroc, cfg

        print(f"Best: {best_cfg}  (val AUROC {best_auroc:.4f})", file=sys.stderr)

    # Retrain on train + val combined. build_X returns a (0,)-shaped array for an empty
    # split, which np.concatenate cannot stack against a 2-D one.
    X_all = np.concatenate([X_train, X_val]) if len(X_val) else X_train
    y_all = np.concatenate([y_train, y_val]) if len(y_val) else y_train
    final_clf = RandomForestClassifier(**best_cfg, random_state=args.seed, n_jobs=-1)
    final_clf.fit(X_all, y_all)

    # Evaluate on both test sets
    results, family_results = [], []
    for name, path in [("test_balanced", args.test_balanced), ("test_realistic", args.test_realistic)]:
        pairs, y_test_raw = read_labelled_csv(path)
        X_test, y_test, kept = build_X(pairs, y_test_raw, embeddings)
        # DDI mode only: rows are instance pairs, so several of them represent one
        # DDI. None here means a PPI-mode CSV with no family level at all.
        fam_pairs = read_family_pairs(path)
        # An empty split, or one that ended up single-class, has no computable AUROC.
        # Report it as a blank row rather than failing the run.
        if len(X_test) == 0 or len(set(y_test.tolist())) < 2:
            reason = "empty" if len(X_test) == 0 else "single-class"
            print(f"{name}: {reason} -- no metrics computed", file=sys.stderr)
            results.append((name, None))
            if fam_pairs is not None:
                # Blank metrics, but still the true DDI count: a single-class split
                # has rows, and reporting 0 there would read as an empty split.
                n_ddis = len({(fam_pairs[i], int(label)) for i, label in zip(kept, y_test)})
                family_results.append((name, None, n_ddis))
            continue
        y_prob = final_clf.predict_proba(X_test)[:, 1]
        metrics = compute_metrics(y_test, y_prob)
        results.append((name, metrics))
        print(
            f"{name}: " + "  ".join(f"{k}={v:.4f}" for k, v in metrics.items()),
            file=sys.stderr,
        )

        if fam_pairs is None:
            continue
        # `kept` realigns the family pairs with the rows that survived build_X.
        y_fam, prob_fam, both = aggregate_families(y_test, y_prob, [fam_pairs[i] for i in kept])
        if both:
            print(
                f"  Warning: {len(both)} family pairs carry both labels, scored as two DDIs each "
                f"(e.g. {sorted(both)[0]})",
                file=sys.stderr,
            )
        # Both classes are present at example level here, and every example keeps
        # its DDI's label, so the family level cannot be single-class -- but the
        # guard costs nothing and keeps the two tables' failure modes identical.
        if len(set(y_fam.tolist())) < 2:
            print(f"{name} (per DDI): single-class -- no metrics computed", file=sys.stderr)
            family_results.append((name, None, len(y_fam)))
            continue
        fam_metrics = compute_metrics(y_fam, prob_fam)
        family_results.append((name, fam_metrics, len(y_fam)))
        print(
            f"{name} (per DDI, {len(y_fam)} family pairs from {len(y_test)} examples): "
            + "  ".join(f"{k}={v:.4f}" for k, v in fam_metrics.items()),
            file=sys.stderr,
        )

    write_mqc(results, args.id)
    if family_results:
        write_family_mqc(family_results, args.id)


if __name__ == "__main__":
    main()
