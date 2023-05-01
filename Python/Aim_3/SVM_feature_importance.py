import numpy as np
from sklearn.svm import SVC

def calculate_feature_importance(X, y, kernel='rbf', C=1.0):
    # Create the SVC classifier with the specified kernel
    clf = SVC(kernel=kernel, C=C, probability=True, max_iter=10000, random_state=42)
    
    # Fit the classifier to the data
    clf.fit(X, y)
    
    # Get the decision function for the original data
    original_decision_function = clf.decision_function(X)
    
    # Initialize an array to store the feature importances
    feature_importances = np.zeros(X.shape[1])
    
    # Calculate the feature importances
    for i in range(X.shape[1]):
        # Create a copy of the data
        X_perturbed = X.copy()
        
        # Perturb the current feature
        np.random.shuffle(X_perturbed[:, i])
        
        # Get the decision function for the perturbed data
        perturbed_decision_function = clf.decision_function(X_perturbed)
        
        # Calculate the change in the decision function
        change = np.mean(np.abs(original_decision_function - perturbed_decision_function))
        
        # Store the feature importance
        feature_importances[i] = change
    
    return feature_importances