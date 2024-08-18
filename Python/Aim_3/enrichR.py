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
    if res.ok:
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
            backgroundType="MSigDB_Hallmark_2020",
        )
    )

    results = res.json()
    
    # Create a dataframe from the results with specific column names
    df = pd.DataFrame(results, columns=[
        'Rank', 'Term name', 'P-value', 'Odds ratio', 'Combined score', 
        'Overlapping genes', 'Adjusted p-value', 'Old p-value', 'Old adjusted p-value'
    ])
    df['celltype'] = celltype
    
    result_list.append(df)

# Combine all dataframes in the result_list
final_df = pd.concat(result_list, ignore_index=True)

# Reorder columns to put celltype first
column_order = ['celltype', 'Rank', 'Term name', 'P-value', 'Odds ratio', 'Combined score', 
                'Overlapping genes', 'Adjusted p-value', 'Old p-value', 'Old adjusted p-value']
final_df = final_df[column_order]

# Export the final dataframe to a CSV file
final_df.to_csv('enrichr_results.csv', index=False)