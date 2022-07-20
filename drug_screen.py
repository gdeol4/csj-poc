from DeepPurpose import DTI as models
import pandas as pd

df_targets = pd.read_csv('disease_targets.csv')
df_drugs = pd.read_csv('cann_mols.csv')

target_name = df_targets['Protein name'].tolist()
target = df_targets.sequence.tolist()
drug_name = df_drugs.Name.tolist()
drug = df_drugs.SMILES.tolist()

# Virtual screening using the trained model or pre-trained model 

net = models.model_pretrained(model = 'Morgan_AAC_BindingDB_IC50')

_ = models.virtual_screening(drug, target, net, drug_name, target_name)