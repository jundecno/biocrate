from __future__ import annotations

import numpy as np


def _to_1d_array(x, *, name: str, dtype=None) -> np.ndarray:
	arr = np.asarray(x, dtype=dtype)
	if arr.ndim != 1:
		raise ValueError(f"{name} must be a 1-D array, got shape={arr.shape}")
	return arr


def _check_same_length(y_true: np.ndarray, y_pred: np.ndarray) -> None:
	if y_true.shape[0] != y_pred.shape[0]:
		raise ValueError(f"y_true and y_pred must have the same length, got {y_true.shape[0]} and {y_pred.shape[0]}")


def _binary_confusion(y_true, y_pred) -> tuple[int, int, int, int]:
	"""Return (tn, fp, fn, tp) for binary labels in {0, 1}."""
	yt = _to_1d_array(y_true, name="y_true").astype(int)
	yp = _to_1d_array(y_pred, name="y_pred").astype(int)
	_check_same_length(yt, yp)

	y_values = set(np.unique(yt).tolist())
	p_values = set(np.unique(yp).tolist())
	if not y_values.issubset({0, 1}) or not p_values.issubset({0, 1}):
		raise ValueError("binary metrics only support labels in {0, 1}")

	tp = int(np.sum((yt == 1) & (yp == 1)))
	tn = int(np.sum((yt == 0) & (yp == 0)))
	fp = int(np.sum((yt == 0) & (yp == 1)))
	fn = int(np.sum((yt == 1) & (yp == 0)))
	return tn, fp, fn, tp


def accuracy_score(y_true, y_pred) -> float:
	yt = _to_1d_array(y_true, name="y_true")
	yp = _to_1d_array(y_pred, name="y_pred")
	_check_same_length(yt, yp)
	return float(np.mean(yt == yp))


def precision_score(y_true, y_pred, *, zero_division: float = 0.0) -> float:
	_, fp, _, tp = _binary_confusion(y_true, y_pred)
	denom = tp + fp
	return float(tp / denom) if denom > 0 else float(zero_division)


def recall_score(y_true, y_pred, *, zero_division: float = 0.0) -> float:
	_, _, fn, tp = _binary_confusion(y_true, y_pred)
	denom = tp + fn
	return float(tp / denom) if denom > 0 else float(zero_division)


def specificity_score(y_true, y_pred, *, zero_division: float = 0.0) -> float:
	tn, fp, _, _ = _binary_confusion(y_true, y_pred)
	denom = tn + fp
	return float(tn / denom) if denom > 0 else float(zero_division)


def f1_score(y_true, y_pred, *, zero_division: float = 0.0) -> float:
	p = precision_score(y_true, y_pred, zero_division=zero_division)
	r = recall_score(y_true, y_pred, zero_division=zero_division)
	denom = p + r
	return float(2 * p * r / denom) if denom > 0 else float(zero_division)


def roc_auc_score(y_true, y_score) -> float:
	"""AUC-ROC for binary labels in {0, 1} with continuous prediction scores."""
	yt = _to_1d_array(y_true, name="y_true").astype(int)
	ys = _to_1d_array(y_score, name="y_score", dtype=float)
	_check_same_length(yt, ys)

	if not set(np.unique(yt).tolist()).issubset({0, 1}):
		raise ValueError("roc_auc_score only supports labels in {0, 1}")

	n_pos = int(np.sum(yt == 1))
	n_neg = int(np.sum(yt == 0))
	if n_pos == 0 or n_neg == 0:
		raise ValueError("roc_auc_score is undefined when y_true has only one class")

	order = np.argsort(-ys)
	y_sorted = yt[order]

	tps = np.cumsum(y_sorted == 1)
	fps = np.cumsum(y_sorted == 0)

	tpr = np.concatenate(([0.0], tps / n_pos, [1.0]))
	fpr = np.concatenate(([0.0], fps / n_neg, [1.0]))
	return float(np.trapz(tpr, fpr))


def pr_auc_score(y_true, y_score) -> float:
	"""Area under Precision-Recall curve for binary labels in {0, 1}."""
	yt = _to_1d_array(y_true, name="y_true").astype(int)
	ys = _to_1d_array(y_score, name="y_score", dtype=float)
	_check_same_length(yt, ys)

	if not set(np.unique(yt).tolist()).issubset({0, 1}):
		raise ValueError("pr_auc_score only supports labels in {0, 1}")

	n_pos = int(np.sum(yt == 1))
	if n_pos == 0:
		raise ValueError("pr_auc_score is undefined when there is no positive sample")

	order = np.argsort(-ys)
	y_sorted = yt[order]

	tps = np.cumsum(y_sorted == 1)
	fps = np.cumsum(y_sorted == 0)

	recall = tps / n_pos
	precision = tps / (tps + fps)

	recall = np.concatenate(([0.0], recall, [1.0]))
	precision = np.concatenate(([1.0], precision, [precision[-1]]))
	return float(np.trapz(precision, recall))


def mae(y_true, y_pred) -> float:
	yt = _to_1d_array(y_true, name="y_true", dtype=float)
	yp = _to_1d_array(y_pred, name="y_pred", dtype=float)
	_check_same_length(yt, yp)
	return float(np.mean(np.abs(yt - yp)))


def mse(y_true, y_pred) -> float:
	yt = _to_1d_array(y_true, name="y_true", dtype=float)
	yp = _to_1d_array(y_pred, name="y_pred", dtype=float)
	_check_same_length(yt, yp)
	return float(np.mean((yt - yp) ** 2))


def rmse(y_true, y_pred) -> float:
	return float(np.sqrt(mse(y_true, y_pred)))


def r2_score(y_true, y_pred) -> float:
	yt = _to_1d_array(y_true, name="y_true", dtype=float)
	yp = _to_1d_array(y_pred, name="y_pred", dtype=float)
	_check_same_length(yt, yp)

	ss_res = np.sum((yt - yp) ** 2)
	ss_tot = np.sum((yt - np.mean(yt)) ** 2)
	if ss_tot == 0:
		raise ValueError("r2_score is undefined when y_true is constant")
	return float(1 - ss_res / ss_tot)


def mape(y_true, y_pred, *, epsilon: float = 1e-12) -> float:
	"""Mean absolute percentage error in [0, +inf)."""
	yt = _to_1d_array(y_true, name="y_true", dtype=float)
	yp = _to_1d_array(y_pred, name="y_pred", dtype=float)
	_check_same_length(yt, yp)

	denom = np.maximum(np.abs(yt), epsilon)
	return float(np.mean(np.abs((yt - yp) / denom)))


def pearsonr(y_true, y_pred) -> float:
	yt = _to_1d_array(y_true, name="y_true", dtype=float)
	yp = _to_1d_array(y_pred, name="y_pred", dtype=float)
	_check_same_length(yt, yp)

	yt_c = yt - np.mean(yt)
	yp_c = yp - np.mean(yp)
	denom = np.sqrt(np.sum(yt_c ** 2) * np.sum(yp_c ** 2))
	if denom == 0:
		raise ValueError("pearsonr is undefined when one input is constant")
	return float(np.sum(yt_c * yp_c) / denom)


def _rankdata(x: np.ndarray) -> np.ndarray:
	"""Rank data with average rank for ties (1-based ranks)."""
	order = np.argsort(x)
	ranks = np.empty_like(x, dtype=float)

	i = 0
	while i < len(x):
		j = i
		while j + 1 < len(x) and x[order[j + 1]] == x[order[i]]:
			j += 1
		avg_rank = (i + j + 2) / 2.0
		ranks[order[i : j + 1]] = avg_rank
		i = j + 1
	return ranks


def spearmanr(y_true, y_pred) -> float:
	yt = _to_1d_array(y_true, name="y_true", dtype=float)
	yp = _to_1d_array(y_pred, name="y_pred", dtype=float)
	_check_same_length(yt, yp)

	return pearsonr(_rankdata(yt), _rankdata(yp))


def evaluate_binary_classification(y_true, y_pred, y_score=None) -> dict[str, float]:
	"""Common binary classification metrics.

	y_pred is hard label in {0, 1}; y_score is optional probability/logit-like score.
	"""
	metrics = {
		"accuracy": accuracy_score(y_true, y_pred),
		"precision": precision_score(y_true, y_pred),
		"recall": recall_score(y_true, y_pred),
		"specificity": specificity_score(y_true, y_pred),
		"f1": f1_score(y_true, y_pred),
	}
	if y_score is not None:
		metrics["roc_auc"] = roc_auc_score(y_true, y_score)
		metrics["pr_auc"] = pr_auc_score(y_true, y_score)
	return metrics


def evaluate_regression(y_true, y_pred) -> dict[str, float]:
	"""Common regression metrics."""
	return {
		"mae": mae(y_true, y_pred),
		"mse": mse(y_true, y_pred),
		"rmse": rmse(y_true, y_pred),
		"r2": r2_score(y_true, y_pred),
		"mape": mape(y_true, y_pred),
		"pearsonr": pearsonr(y_true, y_pred),
		"spearmanr": spearmanr(y_true, y_pred),
	}
