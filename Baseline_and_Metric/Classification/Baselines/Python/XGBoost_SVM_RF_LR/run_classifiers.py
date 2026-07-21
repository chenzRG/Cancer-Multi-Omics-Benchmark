#!/usr/bin/env python3
"""
run_classifiers.py — Classical classifiers (XGBoost / SVM / RF / LR) on MLOmics data.

Usage:
    python run_classifiers.py --dataset GS-BRCA --version Top \
        --data_path /path/to/Classification_datasets \
        --method XGBoost --output_dir /path/to/output
"""
import argparse, os, sys

METHODS = ['XGBoost', 'SVM', 'RF', 'LR', 'all']

def main():
    parser = argparse.ArgumentParser(description='Classical classifiers')
    parser.add_argument('--dataset',    required=True)
    parser.add_argument('--version',    default='Top')
    parser.add_argument('--data_path',  required=True)
    parser.add_argument('--output_dir', default='.')
    parser.add_argument('--method',     default='all', choices=METHODS)
    parser.add_argument('--test_size',  type=float, default=0.2)
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    root_dir   = os.path.normpath(os.path.join(script_dir, '../../../../..'))
    if root_dir not in sys.path:
        sys.path.insert(0, root_dir)

    import numpy as np
    import pandas as pd
    from sklearn.model_selection import train_test_split
    from sklearn.preprocessing import LabelEncoder
    from sklearn.metrics import accuracy_score, f1_score
    from sklearn.svm import SVC
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.linear_model import LogisticRegression
    import xgboost as xgb
    from utils.data_loader import MLOmicsLoader

    # data_path may be Classification_datasets/ or Main_Dataset/ — normalise to Main_Dataset/
    data_root = args.data_path
    if os.path.basename(data_root) == 'Classification_datasets':
        data_root = os.path.dirname(data_root)
    loader = MLOmicsLoader(data_root)
    omics  = loader.load_omics(args.dataset, version=args.version)
    labels = loader.load_labels(args.dataset, version=args.version)

    X = np.concatenate([df.T.values for df in omics.values()], axis=1).astype('float32')
    mask = ~pd.isna(labels)
    X, y = X[mask], np.array(labels)[mask]
    le = LabelEncoder()
    y_enc = le.fit_transform(y)

    X_tr, X_te, y_tr, y_te = train_test_split(
        X, y_enc, test_size=args.test_size, random_state=42, stratify=y_enc)
    print(f'Train: {X_tr.shape}, Test: {X_te.shape}, Classes: {len(le.classes_)}')

    classifiers = {}
    if args.method in ('XGBoost', 'all'):
        classifiers['XGBoost'] = xgb.XGBClassifier(
            n_estimators=100, max_depth=6, use_label_encoder=False,
            eval_metric='mlogloss', random_state=42, n_jobs=-1)
    if args.method in ('SVM', 'all'):
        classifiers['SVM'] = SVC(kernel='rbf', C=1.0, random_state=42)
    if args.method in ('RF', 'all'):
        classifiers['RF'] = RandomForestClassifier(
            n_estimators=100, random_state=42, n_jobs=-1)
    if args.method in ('LR', 'all'):
        classifiers['LR'] = LogisticRegression(
            max_iter=1000, random_state=42, n_jobs=-1)

    results = []
    for name, clf in classifiers.items():
        clf.fit(X_tr, y_tr)
        y_pred = clf.predict(X_te)
        acc = accuracy_score(y_te, y_pred)
        f1  = f1_score(y_te, y_pred, average='macro')
        print(f'{name}: Accuracy={acc:.4f}  F1-macro={f1:.4f}')
        results.append({'method': name, 'accuracy': acc, 'f1_macro': f1})

    os.makedirs(args.output_dir, exist_ok=True)
    out_file = os.path.join(args.output_dir,
                            f'results_{args.dataset}_{args.version}_Classifiers.csv')
    pd.DataFrame(results).to_csv(out_file, index=False)
    print(f'Results saved to: {out_file}')

if __name__ == '__main__':
    main()
