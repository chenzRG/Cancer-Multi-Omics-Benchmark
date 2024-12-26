from sklearn.metrics import precision_score, normalized_mutual_info_score, adjusted_rand_score

def classification_metrics(y_true, y_pred):
    # Calculate precision, NMI, and ARI scores.
    
    # Parameters:
    # y_true (array-like): True labels
    # y_pred (array-like): Predicted labels or cluster labels
    
    # Returns:
    # dict: Dictionary containing precision, NMI, and ARI scores

    precision = precision_score(y_true, y_pred, average='weighted')
    nmi = normalized_mutual_info_score(y_true, y_pred)
    ari = adjusted_rand_score(y_true, y_pred)
    
    return {
        'precision': precision,
        'nmi': nmi,
        'ari': ari
    }

# Example data
# y_true = [0, 1, 1, 0, 1, 0, 1]
# y_pred = [0, 1, 0, 0, 1, 1, 1]

# Reading data from CSV files
classifcation_result = pd.read_csv('classifcation_result.csv')
y_pred = classifcation_result['y_pred'].values

metrics = classification_metrics(y_true, y_pred)
print(f"Precision: {metrics['precision']}")
print(f"NMI: {metrics['nmi']}")
print(f"ARI: {metrics['ari']}")