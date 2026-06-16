import numpy as np
import pandas as pd
import scipy.io
from pathlib import Path

from sklearn.svm import SVC
from sklearn.pipeline import make_pipeline
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    confusion_matrix,
    roc_auc_score,
)


mat_file = Path('LOOCV_nested_features_candidate_edges_fdr.mat')
out_dir = Path('SVM_results_candidate_edges_fdr')
out_dir.mkdir(parents=True, exist_ok=True)


def load_nested_features(mat_path):
    """Load LOOCV nested features saved from MATLAB."""
    try:
        data = scipy.io.loadmat(mat_path)
        X_train_all = np.asarray(data['X_train_all'], dtype=float)
        y_train_all = np.asarray(data['y_train_all'], dtype=int)
        X_test_all = np.asarray(data['X_test_all'], dtype=float)
        y_test_all = np.asarray(data['y_test_all'], dtype=int).ravel()

        edge_count_negative = np.asarray(data['edge_count_negative'], dtype=float)
        edge_count_positive = np.asarray(data['edge_count_positive'], dtype=float)
        edge_count_candidate = np.asarray(data['edge_count_candidate'], dtype=float)

        return (
            X_train_all,
            y_train_all,
            X_test_all,
            y_test_all,
            edge_count_negative,
            edge_count_positive,
            edge_count_candidate
        )

    except Exception as e:
        print('scipy.io.loadmat failed, trying h5py for MATLAB v7.3 file.')
        print('Original error:', e)

        import h5py

        with h5py.File(mat_path, 'r') as f:
            X_train_all = np.array(f['X_train_all'], dtype=float)
            y_train_all = np.array(f['y_train_all'], dtype=int)
            X_test_all = np.array(f['X_test_all'], dtype=float)
            y_test_all = np.array(f['y_test_all'], dtype=int)

            edge_count_negative = np.array(f['edge_count_negative'], dtype=float)
            edge_count_positive = np.array(f['edge_count_positive'], dtype=float)
            edge_count_candidate = np.array(f['edge_count_candidate'], dtype=float)

        if X_train_all.shape[0] == 10:
            X_train_all = np.transpose(X_train_all, (2, 1, 0))

        if X_test_all.shape[0] == 10:
            X_test_all = X_test_all.T

        N = X_test_all.shape[0]

        if y_train_all.shape[0] == N - 1 and y_train_all.shape[1] == N:
            y_train_all = y_train_all.T

        y_test_all = y_test_all.ravel()
        if len(y_test_all) != N:
            y_test_all = y_test_all.reshape(-1)

        if edge_count_negative.shape[0] == 5:
            edge_count_negative = edge_count_negative.T
        if edge_count_positive.shape[0] == 5:
            edge_count_positive = edge_count_positive.T
        if edge_count_candidate.shape[0] == 5:
            edge_count_candidate = edge_count_candidate.T

        return (
            X_train_all,
            y_train_all,
            X_test_all,
            y_test_all,
            edge_count_negative,
            edge_count_positive,
            edge_count_candidate
        )


(
    X_train_all,
    y_train_all,
    X_test_all,
    y_test_all,
    edge_count_negative,
    edge_count_positive,
    edge_count_candidate
) = load_nested_features(mat_file)

X_train_all = np.asarray(X_train_all, dtype=float)
y_train_all = np.asarray(y_train_all, dtype=int)
X_test_all = np.asarray(X_test_all, dtype=float)
y_test_all = np.asarray(y_test_all, dtype=int).ravel()

N = X_test_all.shape[0]

print('Loaded nested LOOCV features')
print('--------------------------------')
print('X_train_all shape:', X_train_all.shape)
print('y_train_all shape:', y_train_all.shape)
print('X_test_all shape:', X_test_all.shape)
print('y_test_all shape:', y_test_all.shape)
print('Number of subjects:', N)

feature_names = [
    'negative_delta', 'negative_theta', 'negative_alpha', 'negative_beta', 'negative_gamma',
    'positive_delta', 'positive_theta', 'positive_alpha', 'positive_beta', 'positive_gamma'
]

bands = ['delta', 'theta', 'alpha', 'beta', 'gamma']

print('\nEmpty-mask summary')
print('--------------------------------')
for i, band in enumerate(bands):
    n_empty_neg = np.sum(edge_count_negative[:, i] == 0)
    n_empty_pos = np.sum(edge_count_positive[:, i] == 0)
    print(f'{band}: negative empty folds = {n_empty_neg}, positive empty folds = {n_empty_pos}')


y_true = []
y_pred = []
y_score = []
coef_all = []

for fold in range(N):
    X_train = np.squeeze(X_train_all[fold, :, :])
    y_train = np.squeeze(y_train_all[fold, :]).astype(int)

    X_test = X_test_all[fold, :].reshape(1, -1)
    y_test = int(y_test_all[fold])

    model = make_pipeline(
        SVC(kernel='linear', C=0.1, probability=False, class_weight='balanced')
    )

    model.fit(X_train, y_train)

    score = model.decision_function(X_test)[0]

    threshold = 0.0
    pred = int(score >= threshold)

    y_true.append(y_test)
    y_pred.append(pred)
    y_score.append(score)

    svm_model = model.named_steps['svc']
    coef_all.append(np.abs(svm_model.coef_[0]))

y_true = np.array(y_true)
y_pred = np.array(y_pred)
y_score = np.array(y_score)
coef_all = np.array(coef_all)

accuracy = accuracy_score(y_true, y_pred)
balanced_acc = balanced_accuracy_score(y_true, y_pred)

cm = confusion_matrix(y_true, y_pred, labels=[0, 1])
tn, fp, fn, tp = cm.ravel()

specificity = tn / (tn + fp) if (tn + fp) > 0 else np.nan
sensitivity = tp / (tp + fn) if (tp + fn) > 0 else np.nan

roc_auc = roc_auc_score(y_true, y_score)

print('\nNested LOOCV SVM performance')
print('--------------------------------')
print(f'Accuracy:          {accuracy:.4f}')
print(f'Balanced accuracy: {balanced_acc:.4f}')
print(f'Sensitivity:       {sensitivity:.4f}')
print(f'Specificity:       {specificity:.4f}')
print(f'AUC:               {roc_auc:.4f}')
print('\nConfusion matrix [HC=0, MCI=1]')
print(cm)

rng = np.random.default_rng(2026)
n_bootstraps = 5000
boot_aucs = []

for _ in range(n_bootstraps):
    indices = rng.integers(0, len(y_true), len(y_true))

    if len(np.unique(y_true[indices])) < 2:
        continue

    boot_auc = roc_auc_score(y_true[indices], y_score[indices])
    boot_aucs.append(boot_auc)

boot_aucs = np.asarray(boot_aucs)

auc_ci_low = np.percentile(boot_aucs, 2.5)
auc_ci_high = np.percentile(boot_aucs, 97.5)

print(f'AUC = {roc_auc:.4f}')
print(f'95% CI = [{auc_ci_low:.4f}, {auc_ci_high:.4f}]')

pred_df = pd.DataFrame({
    'subject_index': np.arange(1, N + 1),
    'true_label': y_true,
    'predicted_label': y_pred,
    'decision_score': y_score
})

pred_df.to_excel(out_dir / 'nested_LOOCV_predictions.xlsx', index=False)

performance_df = pd.DataFrame({
    'metric': [
        'accuracy',
        'balanced_accuracy',
        'sensitivity',
        'specificity',
        'AUC',
        'TN',
        'FP',
        'FN',
        'TP'
    ],
    'value': [
        accuracy,
        balanced_acc,
        sensitivity,
        specificity,
        roc_auc,
        tn,
        fp,
        fn,
        tp
    ]
})

performance_df.to_excel(out_dir / 'nested_LOOCV_performance.xlsx', index=False)

mean_abs_weight = coef_all.mean(axis=0)
std_abs_weight = coef_all.std(axis=0)
se_abs_weight = std_abs_weight / np.sqrt(N)

feature_df = pd.DataFrame({
    'feature': feature_names,
    'mean_abs_weight': mean_abs_weight,
    'se_abs_weight': se_abs_weight
})

feature_df = feature_df.sort_values('mean_abs_weight', ascending=False)
feature_df.to_excel(out_dir / 'nested_LOOCV_feature_weights.xlsx', index=False)

edge_summary = pd.DataFrame()

for i, band in enumerate(bands):
    edge_summary.loc[band, 'candidate_mean'] = np.mean(edge_count_candidate[:, i])
    edge_summary.loc[band, 'candidate_min'] = np.min(edge_count_candidate[:, i])
    edge_summary.loc[band, 'candidate_max'] = np.max(edge_count_candidate[:, i])

    edge_summary.loc[band, 'negative_mean'] = np.mean(edge_count_negative[:, i])
    edge_summary.loc[band, 'negative_min'] = np.min(edge_count_negative[:, i])
    edge_summary.loc[band, 'negative_max'] = np.max(edge_count_negative[:, i])
    edge_summary.loc[band, 'negative_empty_folds'] = np.sum(edge_count_negative[:, i] == 0)

    edge_summary.loc[band, 'positive_mean'] = np.mean(edge_count_positive[:, i])
    edge_summary.loc[band, 'positive_min'] = np.min(edge_count_positive[:, i])
    edge_summary.loc[band, 'positive_max'] = np.max(edge_count_positive[:, i])
    edge_summary.loc[band, 'positive_empty_folds'] = np.sum(edge_count_positive[:, i] == 0)

edge_summary.to_excel(out_dir / 'nested_LOOCV_edge_count_summary.xlsx')
print('\nFiles saved to:', out_dir)
