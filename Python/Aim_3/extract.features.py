import pickle
import pandas as pd
import glob
import os.path

# loop through the models in ML.models and save features as .txt file
for file in glob.glob('ML.models/*.sav'):
    print(file)
    # load the model from disk
    model = pickle.load(open(file, 'rb'))
    # Save numpy.ndarray as txt file
    df = pd.DataFrame(model.get_feature_names_out(), columns = ['Features'])
    df.to_csv('ML.models/features/'+os.path.basename(file).replace('.sav', '.txt'), index = False, sep='\t')