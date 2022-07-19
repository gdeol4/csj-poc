import streamlit as st

code_import = ('''
# Install biopython and pandas and import the following:
import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser

# Open the fasta file
fasta = open("human_proteome.fasta")

# The SimpleFastaParser returns a generator object and to access all rows, convert it to a list.
fasta_list = list(SimpleFastaParser(fasta))

# Read in a csv file of cannflavin molecules (containing names and SMILES)
df_drugs = pd.read_csv('cann_mols.csv')  

              ''')

code_df_init = ('''
# create dataframe from list with named columns:

protein_df = pd.DataFrame(fasta_list, columns=['info', 'sequence'])
               ''')

code_functions = ('''

def string_splitter(string):
   string = string.split("HUMAN",1)[1] # only keep string contents after "HUMAN"
   string = string.split("OS=Homo sapiens") # only keep string contents after "OS=Homo sapiens"
    
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
    df.columns = ['type', 'id', 'protein_name', 'drop', 'Species', 'Gene', 'PE', 'Mutation'] # rename columns
    df.drop("drop", axis=1, inplace=True) # drop empty column
    df = pd.concat([df, dfx], axis=1) # merge dataframes
    df.drop("info", axis=1, inplace=True) # drop info column
    df.drop('Species', axis=1, inplace=True) # drop column
    df['Gene'] = df['Gene'].str[3:] # strip first 3 characters
    
    return df                
                 
                 ''')

code_pipeline = ('''
clean_proteins =(protein_df.
                 pipe(info_parser).
                 pipe(info_pre_processed).
                 pipe(info_processed))
                ''')

code_lists = ('''
df_targets = clean_proteins[['protein_name','sequence']].copy() # making a copy of the clean dataframe.

df_targets['drug_name'] = df_drugs['Name'][0] # adding a column with a constant value of "cannaflavin a".
df_targets['SMILES'] = df_drugs['SMILES'][0] # adding a column with a constant value of "cannaflavin a" SMILES string.
             ''')

code_make_lists = ('''
target_name = df_targets.protein_name.tolist()
target = df_targets.sequence.tolist()
drug_name = df_targets.drug_name.tolist()
drug = df_targets.SMILES.tolist()
''')

code_save = ('''
df_targets.to_csv('proteome_clean.csv') # save dataframe as csv
''')

code_entire = ('''
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
    df = pd.concat([df, protein_df], axis=1) # merge dataframes
    df.drop("info", axis=1, inplace=True) # drop info column
    df.drop('Species', axis=1, inplace=True) # drop column
    df['Gene'] = df['Gene'].str[3:] # strip first 3 characters
    
    return df

clean_proteins =(protein_df.
                 pipe(info_parser).
                 pipe(info_pre_processed).
                 pipe(info_processed))

df_targets = clean_proteins[['Protein name','sequence']].copy() # making a copy of the clean dataframe.

df_targets['drug_name'] = df_drugs['Name'][0] # adding a column with a constant value of "cannaflavin a".
df_targets['SMILES'] = df_drugs['SMILES'][0] # adding a column with a constant value of "cannaflavin a" SMILES string.

target_name = df_targets.['Protein name'].tolist()
target = df_targets.sequence.tolist()
drug_name = df_targets.drug_name.tolist()
drug = df_targets.SMILES.tolist()

df_targets.to_csv('proteome_clean.csv') # save dataframe as csv
''')

st.title('Preprocessing the human proteome to use with Deep Purpose')
st.header('''Parsing a FASTA file into a dataframe''')

st.write('''Download the human proteome (~20600 protein amino acid sequences) as a FASTA file from Uniprot:
            https://www.uniprot.org/proteomes/UP000005640
         ''')
st.write('''
            Rename the file "human_proteome.fasta"
         ''')

st.write('''
            Move "cann_mols.csv" into the working directory as well.
         ''')

st.subheader('Importing libraries and reading in the data')
st.code(code_import, language = 'python')

st.subheader('Creating a dataframe from a list of tuples')
st.write('The SimpleFastaParser converts each entry in the fasta file into a tuple containing "title" and "sequence" information. The code below creates a dataframe with columns "info" and "sequence": ' )
st.code(code_df_init, language = 'python')

st.subheader('Writing the functions to format the dataframe:')
st.write('There are four functions:')
st.write('1. string_splitter - this function formats string values in a list')
st.write('2. info_parser - this function splits the dataframe into three initial columns')
st.write('3. info_pre_processed - this function extracts and formats protein names')
st.write('4. info_processed - this function extracts gene names, uniprot ids, and other info')
st.code(code_functions, language = 'python')

st.subheader('Structuring the functions as a pipeline')
st.write('The pandas pipe() function allows for building pipelines like below. The functions can also be called individually.')
st.code(code_pipeline, language = 'python')

st.subheader('Preparing Deep Purpose input data')
st.write('The Deep Purpose virtual screening function takes four lists of values as an input. All four lists have to be the same length:')
st.code(code_lists, language = 'python')

st.subheader('Making the lists')
st.code(code_make_lists, language = 'python')

st.subheader('Saving the CSV file:')
st.code(code_save, language = 'python')


st.subheader('Script as a single file:')
st.code(code_entire, language = 'python')