#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 18:10:36 2024

@author: patrick
"""
import os

directory = os.path.dirname(os.path.abspath(__file__))
os.chdir(directory)

import pandas as pd
import numpy as np
import trailblazer as tb
import GTExGenie as gtx
import pickle
#%%
#read in human_proteins and gtex data
human_proteins = pickle.load(open('human_proteins_LoD.p','rb'))

gtex = pd.read_csv('/Users/patrick/Documents/CHBE194/Homework4/GTEX_Uniprot/GTEX_data_2020.gct',sep='\t',skiprows=2,index_col=0)

#%%
#use trailblazer to get a data frame of human genes from KEGG
human_genes = tb.getKeggGenes('hsa')

#use GTEx Genie to match the genes listed in gtex to KEGG ID
gtex_kegg_dictionary = gtx.KeggMatch(gtex, human_genes)
    
#%%
#read in attributes for gtex file
attributes = pd.read_excel('/Users/patrick/Documents/CHBE194/Homework4/GTEX_Uniprot/Sample_attributes_v8.xlsx')

#%%

#unique_tissues = list(set(attributes['SMTSD']))
unique_tissues = gtx.tissueTypes(attributes)
unique_tissues = list(unique_tissues)

#%%
#determine specificity of each transcript, where 1 is most specific and 54 is least specific
##This line of code will take several hours to complete, so I saved the output as a pickle

#specificity = gtx.transcriptSpecificity(gtex,attributes)
specificity = pickle.load(open('specificity.pkl','rb'))
#%%
tb.cytoscape(unique_tissues, gtex, attributes, )

#%%
#these genes are gathered from PCA analysis of patient transcriptomic data
predictive_obesity_genes = pickle.load(open('top_1000.pickle','rb'))
#%%

tb.cytoscape([unique_tissues[0]],
             gtex,
             attributes,
             nModules = 1000,
             KEGG_query = 'secreted',
             addOn = {'Specificity':specificity,
                      'Predictive Genes':predictive_obesity_genes})


#%%
tb.cytoscape(unique_tissues[45:],
             gtex,
             attributes,
             KEGG_query='secreted',
             nModules=20000,
             addOn = {'Specificity':specificity})

#%%
script_path = os.path.abspath(__file__)

tb.
#%%
script_directory = os.path.dirname(os.path.abspath(__file__))

#%%
gtex = pd.read_csv('/Users/patrick/Documents/CHBE194/Homework4/GTEX_Uniprot/GTEX_data_2020.gct',sep='\t',skiprows=2,index_col=0)
human_kegg = '/Users/patrick/Documents/CHBE194/Human_KEGG_cytoscape_06182020.p'
human_kegg = pickle.load(open(human_kegg, 'rb'))
human_attributes = human_kegg['df_attributes']
human_attributes_copy = human_attributes.copy()
human_proteins = '/Users/patrick/Documents/CHBE194/human_proteins_LoD.p'
human_proteins = pickle.load(open(human_proteins, 'rb'))
human_genes = tb.getKeggGenes('hsa')
keggMatch = gtx.KeggMatch(gtex, human_genes)


addOn = {'Specificity':specificity,
         'Predictive_Genes':predictive_obesity_genes}

tissueTranscripts = enrichedTranscripts

num_colors = max(specificity)-min(specificity)
colors = tb.colorPalette(color = 'yellow', nColors = num_colors)

for key in list(addOn.keys()):
    if 'specificity'.upper() == key.upper() and len(addOn[key]) == len(tissueTranscripts):
        tissueTranscripts['Specificity'] = addOn[key]
        
        human_attributes_copy['Specificity'] = 0
        for idx1, row_attribute in human_attributes_copy.iterrows():
            node_id = row_attribute['Node_ID']
            if node_id in list(tissueTranscripts['KEGG_ID']):
                match = tissueTranscripts[tissueTranscripts['KEGG_ID']==node_id].iloc[0]
                spec = match.loc['Specificity']
                row_attribute['Node_Color'] = colors[spec-1]
                row_attribute['Specificity'] = spec
        
    else:
        pass
    

     
#%%
merged_df = pd.merge(human_attributes, tissueTranscripts, left_on='Node_ID', right_on='KEGG_ID', how='left')
merged_df = merged_df.drop_duplicates(subset='Node_ID')

#%%

for key in list(addOn.keys()):
    #make colors that match specificity values
    if 'specificity'.upper() in key.upper():
        merged_df[key]
        num_colors = int(merged_df[key].max(skipna=True)-merged_df[key].min(skipna=True))
        colors = tb.colorPalette(color = 'yellow', nColors = num_colors)

        for index, row in merged_df.iterrows():
            #KEGG_ID = row['Node_ID']
            specificity_value = row[key]
            
            if specificity_value != 0 and math.isnan(specificity_value) == False:
                row['Node_Color'] = colors[int(specificity_value)-1]






#%%
TPM_values = [entry['TPM'] for entry in adipose_transcripts if entry['Significant']=='Yes']

#%%
import matplotlib.pyplot as plt

plt.hist(TPM_values,bins = len(TPM_values),color='skyblue')
#plt.xscale('log')
plt.xlim(0,1000)
plt.ylim(0,2000)
plt.ylabel('Transcripts per million (TPM)')
plt.xlabel('Number of Transcripts')
plt.title('Frequency of TPM count in adipose-enriched transcripts')
#plt.xticks([1,10,100,1000,10000])

#%%
TPM_size = []
for entry in adipose_transcripts.values():
    if entry > 1000:
        TPM_size.append('large')
    elif entry > 600:
        TPM_size.append('medium')
    elif entry > 200:
        TPM_size.append('small')
    else:
        TPM_size.append('extra_small')
        
#%%
gtex_adipose_LoD = []

for key in adipose_transcripts:
    if key in gtex.index:
        Gene_Name = gtex['Description'][key]
        HSA_ID = gtex_kegg_dictionary[Gene_Name]
        TPM = adipose_transcripts[key]
        
        
for key in gtex_kegg_dictionary:
    idx = gtex['Description'].index[gtex['Description']==key][0]
    if idx in adipose_transcripts.keys():
        gtex_adipose_dictionary[key] = gtex_kegg_dictionary[key]