import numpy as np
from sklearn.svm import SVC
import matplotlib.pyplot as plt
from sklearn.model_selection import RepeatedKFold

def calculate_feature_importance_cv(X, y, kernel='rbf', C=1.0, n_splits=5, n_repeats=10):
    # Create the SVC classifier with the specified kernel
    clf = SVC(kernel=kernel, C=C, probability=True, max_iter=10000, random_state=42)
    
    # Create the RepeatedKFold object
    rkf = RepeatedKFold(n_splits=n_splits, n_repeats=n_repeats, random_state=42)
    
    # Initialize an array to store the feature importances
    feature_importances = np.zeros(X.shape[1])
    
    # Perform cross-validation
    for train_index, test_index in rkf.split(X):
        # Split the data into training and test sets
        X_train, X_test = X.iloc[train_index], X.iloc[test_index]
        y_train, y_test = y.iloc[train_index], y.iloc[test_index]
        
        # Fit the classifier to the training data
        clf.fit(X_train, y_train)
        
        # Get the decision function for the test data
        original_decision_function = clf.decision_function(X_test)
        
        # Calculate the feature importances
        for i in range(X.shape[1]):
            # Create a copy of the test data
            X_test_perturbed = X_test.copy()
            
            # Perturb the current feature
            np.random.shuffle(X_test_perturbed.iloc[:, i].values)
            
            # Get the decision function for the perturbed test data
            perturbed_decision_function = clf.decision_function(X_test_perturbed)
            
            # Calculate the change in the decision function
            change = np.mean(np.abs(original_decision_function - perturbed_decision_function))
            
            # Update the feature importance
            feature_importances[i] += change
    
    # Average the feature importances over all cross-validation folds and repeats
    feature_importances /= (n_splits * n_repeats)
    
    return feature_importances