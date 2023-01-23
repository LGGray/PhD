# Import required libraries
import sys
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedShuffleSplit


file = sys.argv[0]

# Read in the expression matrix as a pandas dataframe
df = pd.read_csv(file, delimiter='\t', index_col=0)

# Replace classes with binary label
df['class'] = df['class'].replace({"OA": 0, "RA": 1})

# Split the data into features (X) and target (y)
X = df.iloc[:, :-1]
y = df.iloc[:, -1]

# Create the stratified sampling object
sss = StratifiedShuffleSplit(n_splits=1, test_size=0.2, random_state=42)

# Loop through the splits
for train_index, test_index in sss.split(X, y):
    # Get the training and testing data
    X_train, X_test = X.iloc[train_index], X.iloc[test_index]
    y_train, y_test = y.iloc[train_index], y.iloc[test_index]

# Build the logistical regression model
model = LogisticRegression()
model.fit(X_train, y_train)

# Predict the target for new data
y_pred = model.predict(X_test)
score = model.score(X_test, y_test)
print(score)

print(prediction)