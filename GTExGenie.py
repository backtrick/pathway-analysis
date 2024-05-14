#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 30 20:22:54 2024

@author: patrick
"""

def KeggMatch(gtex, human_genes):
    
    """Makes a dictionary matching gtex transcripts (key) to KEGG ID (value)

  Args:
      gtex (DataFrame, Series, list): A DataFrame containing the gtex data
          or a pandas Series object containing GTEx transcripts
          or a list containing GTEx transcripts
      human_genes (DataFrame): A DataFrame containing KEGG IDs, Gene Names, and Description

  Returns:
      gtex_kegg_dictionary: (dictionary) A dictionary in which the 
          key is the gtex transcript and the
          value is the KEGG ID
  """

    import pandas as pd
    
    if type(gtex) == pd.DataFrame:
        gtex_kegg_dictionary = {}
        
        for transcript in gtex['Description']:
            transcript_name = transcript.upper()
            
            for index, gene_name in enumerate(human_genes['Gene_Names']):
                if type(gene_name) == str and transcript_name in gene_name:
                    KEGG_ID = human_genes['KEGG_ID'].iloc[index]
                    break
                else:
                    KEGG_ID = ''
            gtex_kegg_dictionary[transcript] = KEGG_ID
            
        return gtex_kegg_dictionary
    
    elif type(gtex) == list or type(gtex) == pd.Series:
        gtex_kegg_dictionary = {}
        
        for transcript in gtex:
            transcript_name = transcript.upper()
            
            for index, gene_name in enumerate(human_genes['Gene_Names']):
                if type(gene_name) == str and transcript_name in gene_name:
                    KEGG_ID = human_genes['KEGG_ID'].iloc[index]
                    break
                else:
                    KEGG_ID = ''
            gtex_kegg_dictionary[transcript] = KEGG_ID
            
        return gtex_kegg_dictionary
    


def tissueTypes(attributes):
    
    """Returns a set list of all tissue types in gtex attributes

  Args:
      attributes (DataFrame): A DataFrame containing gtex attributes

  Returns:
      tissues: (list) A set list of all tissue types included in GTEx dataset
  """
    tissues = list(set(attributes['SMTSD']))
    
    return tissues
    
    

def tissueFilter(tissueType, gtex, attributes):
    
    """Filter GTEX data to include only samples of a specific tissue type.

  Args:
      tissue_type (string): A string containing the tissue type of interest
      gtex (DataFrame): A DataFrame containing the gtex data
      attributes (DataFrame): A DataFrame containing the gtex attributes

  Returns:
      DataFrame: A DataFrame contaning TPM values for tissue of interest
      
  """
  
    import re
    
    tissueType = re.sub(r'([()])', r'\\\1',tissueType)
  
    tissue_attributes = attributes[attributes['SMTSD'].str.contains(tissueType)]

    tissue_sampleids = list(tissue_attributes['SAMPID'])

    valid_columns = [col for col in tissue_sampleids if col in gtex.columns]

    TPM_tissue = gtex[valid_columns]
    
    return TPM_tissue



def tissueTranscripts(tissue_type, gtex, attributes):
    
    """Counts the number of genes expressed by a specific tissue type.

  Args:
      tissue_type (string): A string containing the tissue type of interest
      attributes (DataFrame): A DataFrame containing the gtex attributes
      gtex (DataFrame): A DataFrame containing the gtex data

  Returns:
      int: The number of genes expressed by the tissue of interest
  """
    
    import scipy.stats as sp
    
    tissue_attributes = attributes[attributes['SMTSD'].str.contains(tissue_type)]

    tissue_sampleids = list(tissue_attributes['SAMPID'])
    
    valid_columns = [col for col in tissue_sampleids if col in gtex.columns]

    TPM_tissue_type = gtex[valid_columns]
    
    #p_values = []
    tissue_transcripts = []

    for index, row in TPM_tissue_type.iterrows():
        curr_dict = {}
        mean_TPM = sum(row) / len(row)
        transcript_name = index
        t_statistic, p_value = sp.ttest_1samp(row, 10, alternative='greater')
        
        if p_value < 0.05/len(TPM_tissue_type): # Bonferroni post-hoc correction
            significant = 1
        else:
            significant = 0
        
        curr_dict = {'Transcript':transcript_name,
                     'P-value':p_value,
                     'TPM':mean_TPM,
                     'Significant':significant}
        tissue_transcripts.append(curr_dict)
        
    return sum([entry['Significant'] for entry in tissue_transcripts])



def uniqueTissues(attributes):
        
    """Provides set of tissue types in GTEx data set

  Args:
      attributes (DataFrame): A DataFrame containing the GTEx attributes

  Returns:
      unique_tissue (set): A set of unique tissue in GTEx data set
  """
    tissues = attributes['SMTSD']
    
    unique_tissues = set(list(tissues))
    
    unique_tissues = list(unique_tissues)
        
    return unique_tissues



def enrichedTranscripts(gtex, TPM_data, TPM_threshold = 10, post_hoc = 'Bonferroni'):
    
    """Identifies GTEx transcripts that are enriched in a data set, typically
        a data set of a specific tissue samples.

  Args:
      gtex (DataFrame): A DataFrame containing the entire GTex data set
      TPM_data (DataFrame): A DataFrame containing the gtex data from a tissue type
      TPM_threshold (int): The transcripts per million threshold (defaults to 10)
      post_hoc (string): The type of post-hoc correction used in multiple comparison
          (defaults to Bonferroni)

  Returns:
      tissue_transcripts (DataFrame): A DataFrame containing columns 'Transcript', 'P-value', 
          'TPM', and 'Significant'
  """
    
    import scipy.stats as sp
    import pandas as pd
    import GTExGenie as gtx
    import trailblazer as tb
    
    human_genes = tb.getKeggGenes('hsa')
    keggMatch = gtx.KeggMatch(gtex, human_genes)
    
    #keggMatch = pd.DataFrame(list(keggMatch.items()),columns=['Gene','KEGG_ID'])
    tissue_transcripts = []

    for index, row in TPM_data.iterrows():
        curr_dict = {}
        mean_TPM = sum(row) / len(row)
        transcript_name = gtex['Description'][index]
        t_statistic, p_value = sp.ttest_1samp(row, TPM_threshold, alternative='greater')
        
        if post_hoc == 'Bonferroni':
            if p_value < 0.05/len(TPM_data): # Bonferroni post-hoc correction
                significant = 'Yes'
            else:
                significant = 'No'
                
        
        curr_dict = {'Transcript':transcript_name,
                     'P-value':p_value,
                     'TPM':mean_TPM,
                     'Significant':significant,
                     'KEGG_ID':keggMatch[transcript_name]}
        tissue_transcripts.append(curr_dict)
    
    tissue_transcripts = pd.DataFrame(tissue_transcripts, index=gtex.index)
    
    #tissue_transcripts = tissue_transcripts[tissue_transcripts['Significant']=='Yes']
    
    return tissue_transcripts



def transcriptSpecificity(gtex, attributes):
    """Determines how specific each gene's expression is relative to other tissues.
        A specificity of 1 means the gene is only expressed in 1 tissue type, whereas
        a specificity of 54 means the gene is expressed in all 54 tissue types.
        
    Args:
        gtex (DataFrame): A DataFrame containing the entire GTex data set  
        attributes (DataFrame): A DataFrame containing the GTEx attributes
        
    
    """

    import GTExGenie as gtx
    import scipy.stats as sp
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
    # Perform operations that trigger the warning here

    
        synth_gtex = gtex.copy()
        tissue_types = gtx.uniqueTissues(attributes)
    
        synth_gtex['Specificity'] = 0
        
        for tissue in tissue_types:
            print(str(tissue_types.index(tissue)/len(tissue_types)) + '% complete')
            tissue_transcripts = gtx.tissueFilter(tissue, gtex, attributes)
            for index, row in tissue_transcripts.iterrows():
                #mean_TPM = sum(gtex.loc[index]) / len(gtex.loc[index])
                t_statistic, p_value = sp.ttest_1samp(tissue_transcripts.loc[index][1:], 10, alternative='greater')
                
                if p_value < 0.05/len(synth_gtex.loc[index]): # Bonferroni post-hoc correction
                    synth_gtex['Specificity'].loc[index] += 1
                else:
                    synth_gtex['Specificity'].loc[index] += 0
    
    return synth_gtex['Specificity']

        