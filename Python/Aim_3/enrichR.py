import json
import requests
import pyreadr
import pandas as pd

# Load the JSON file
with open('figures/all_features.json', 'r') as file:
    all_features = json.load(file)

result_list = []
for celltype in all_features.keys():
    mtx = pyreadr.read_r(f'{celltype}.RDS')
    mtx = mtx[None]
    background = mtx.columns[5:].tolist()
    genes = all_features[celltype]

    # Remove cellCount, age, and ancestry from genes
    genes = [gene for gene in genes if gene not in ['cellCount', 'age', 'ancestry']]

    base_url = "https://maayanlab.cloud/speedrichr"

    # Add the gene list
    res = requests.post(
        base_url + '/api/addList',
        files=dict(
            list=(None, '\n'.join(genes)),
            description=(None, 'gene list'),
        )
    )
    userlist_response = res.json()

    # Add the background list
    res = requests.post(
        base_url + '/api/addbackground',
        data=dict(background='\n'.join(background)),
    )
    background_response = res.json()

    # Perform enrichment analysis with the background
    res = requests.post(
        base_url + '/api/backgroundenrich',
        data=dict(
            userListId=userlist_response.get('userListId'),
            backgroundid=background_response.get('backgroundid'),
            backgroundType="GO_Biological_Process_2023",
        )
    )
    results = res.json()
    
    # Create a dataframe from the results
    df = pd.DataFrame()
    for result in results['GO_Biological_Process_2023']:
        tmp = pd.DataFrame({
            'Rank': [result[0]],
            'Term name': [result[1]],
            'P-value': [result[2]],
            'Odds ratio': [result[3]],
            'Combined score': [result[4]],
            'Overlapping genes': [', '.join(result[5])],
            'Adjusted p-value': [result[6]],
            'Old p-value': [result[7]],
            'Old adjusted p-value': [result[8]],
            'celltype': [celltype]
        })
        df = pd.concat([df, tmp], ignore_index=True)
    result_list.append(df)

# Combine all dataframes in the result_list
final_df = pd.concat(result_list, ignore_index=True)

# Reorder columns to put celltype first
column_order = ['celltype', 'Rank', 'Term name', 'P-value', 'Odds ratio', 'Combined score', 
                'Overlapping genes', 'Adjusted p-value', 'Old p-value', 'Old adjusted p-value']
final_df = final_df[column_order]

# Export the final dataframe to a CSV file
final_df.to_csv('figures/enrichr_results_GOBP.csv', index=False)