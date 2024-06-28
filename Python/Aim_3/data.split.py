import pyreadr
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split

# Get the file name from the command line
# Use Tcm.Naive.helper.T.cells because all individuals are included
file = 'new_pseudobulk/Tcm.Naive.helper.T.cells.chrX.RDS'

# Read in expression RDS file
df = pyreadr.read_r(file)
df = df[None]
print(df.head())

df['class'] = df['class'].replace({"control": 0, "disease": 1})
# Remove 'Hispanic or Latin American' ancestry
df = df[df['ancestry'] != 'Hispanic or Latin American']

# Iterate ten times and save index for each iteration
train_indexes = []
test_indexes = []
i=1
while i < 11:
    # Split the data into training and testing sets
    df_train, df_test = train_test_split(df, test_size=0.2, shuffle=True, stratify=df[["class", "ancestry"]], random_state=i)
    train_indexes.append(df_train.index)
    test_indexes.append(df_test.index)
    i += 1

# Write index to file
for i in range(0, 10):
    pd.DataFrame(train_indexes[i]).to_csv(f"new_pseudobulk/split_{i+1}/train_index.csv", index=False)
    pd.DataFrame(test_indexes[i]).to_csv(f"new_pseudobulk/split_{i+1}/test_index.csv", index=False)