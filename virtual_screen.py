import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser

fasta = open("human_proteom.fasta")
df_drugs = pd.read_csv('cann_mols.csv')  

fasta_list = list(SimpleFastaParser(fasta))

protein_df = pd.DataFrame(fasta_list, columns=['info', 'sequence'])

def string_splitter(string):
    string = string.split("HUMAN",1)[1]
    string = string.split("OS=Homo sapiens")
    
    return string

def info_parser(dfx):
    df = dfx.copy()
    df = df['info'].str.split('|', expand=True) # split on the "|" character
    df.columns = ['type', 'id', 'info'] # rename the three columns

    return df

def info_pre_processed(dfx):
    df = dfx.copy() # create a copy of df
    list_of_info = list(dfx['info']) # convert column values to list
    list_of_info = [string_splitter(x) for x in list_of_info] # apply string_splitter function to list elements
    df['temp_col'] = list_of_info # create a temporary column from the processed list
    df.drop("info", axis=1, inplace=True) # drop info column
    split_df = pd.DataFrame(df['temp_col'].tolist(), columns=['Protein name', 'info']) # creating this dataframe to merge back onto processed dataframe
    df.drop("temp_col", axis=1, inplace=True) # dropping temp_col
    df = pd.concat([df, split_df], axis=1) # merging both dataframes

    return df
    
def info_processed(dfx):
    df = dfx.copy() # create a copy of the dataframe
    split_df = df['info'].str.split(' ', expand=True) # creating a seperate dataframe to merge with
    df = pd.concat([df, split_df], axis=1) # merge dataframes
    df.drop("info", axis=1, inplace=True) # drop info column
    df.columns = ['type', 'id', 'Protein name', 'drop', 'Species', 'Gene', 'PE', 'Mutation'] # rename columns
    df.drop("drop", axis=1, inplace=True) # drop empty column
    df = pd.concat([df, dfx], axis=1) # merge dataframes
    df.drop("info", axis=1, inplace=True) # drop info column
    df.drop('Species', axis=1, inplace=True) # drop column
    df['Gene'] = df['Gene'].str[3:] # strip first 3 characters
    
    return df

parsed_proteins = INFO_parser(protein_df)

df_targets = parsed_proteins[['id','sequence']].copy()

df_targets['drug_name'] = df_drugs['Name'][0]
df_targets['SMILES'] = df_drugs['SMILES'][0]

df_targets_2 = df_targets.copy()
df_targets_2['drug_name'] = df_drugs['Name'][1]
df_targets_2['SMILES'] = df_drugs['SMILES'][1]

df_targets_3 = df_targets.copy()
df_targets_3['drug_name'] = df_drugs['Name'][2]
df_targets_3['SMILES'] = df_drugs['SMILES'][2]

df_targets_4 = df_targets.copy()
df_targets_4['drug_name'] = df_drugs['Name'][3]
df_targets_4['SMILES'] = df_drugs['SMILES'][3]

df_targets_5 = df_targets.copy()
df_targets_5['drug_name'] = df_drugs['Name'][4]
df_targets_5['SMILES'] = df_drugs['SMILES'][4]

df_targets_6 = df_targets.copy()
df_targets_6['drug_name'] = df_drugs['Name'][5]
df_targets_6['SMILES'] = df_drugs['SMILES'][5]

df_targets_7 = df_targets.copy()
df_targets_7['drug_name'] = df_drugs['Name'][6]
df_targets_7['SMILES'] = df_drugs['SMILES'][6]

df_targets_8 = df_targets.copy()
df_targets_8['drug_name'] = df_drugs['Name'][7]
df_targets_8['SMILES'] = df_drugs['SMILES'][7]

df_targets_9 = df_targets.copy()
df_targets_9['drug_name'] = df_drugs['Name'][8]
df_targets_9['SMILES'] = df_drugs['SMILES'][8]

dfxx = pd.concat([df_targets, 
                  df_targets_2, 
                  df_targets_3, 
                  df_targets_4, 
                  df_targets_5, 
                  df_targets_6, 
                  df_targets_7, 
                  df_targets_8, 
                  df_targets_9], axis=0)

target_name = dfxx.id.tolist()
target = dfxx.sequence.tolist()
drug_name = dfxx.drug_name.tolist()
drug = dfxx.SMILES.tolist()

dfxx.to_csv('proteome_clean.csv')

import pandas as pd
import re

dfxx = pd.read_csv('proteome_clean.csv')

target_name = dfxx.id.tolist()
target = dfxx.sequence.tolist()
drug_name = dfxx.drug_name.tolist()
drug = dfxx.SMILES.tolist()

import pandas as pd
from DeepPurpose import DTI as models

dfxx = pd.read_csv('proteome_clean.csv')

target_name = dfxx.id.tolist()
target = dfxx.sequence.tolist()
drug_name = dfxx.drug_name.tolist()
drug = dfxx.SMILES.tolist()

# Virtual screening using the trained model or pre-trained model 

net = models.model_pretrained(model = 'MPNN_CNN_DAVIS')

_ = models.virtual_screening(drug, target, net, drug_name, target_name)



import pandas as pd
from DeepPurpose import DTI as models
targets = pd.read_csv('targets_pre_screen.csv')

target_name = targets['Protein name'].tolist()
target = targets.sequence.tolist()
drug_name = targets.drug_name.tolist()
drug = targets.SMILES.tolist()

# Virtual screening using the trained model or pre-trained model 

net = models.model_pretrained(model = 'Morgan_AAC_BindingDB_IC50')

_ = models.virtual_screening(drug, target, net, drug_name, target_name)