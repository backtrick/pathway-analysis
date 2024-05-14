#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 21:23:40 2024

@author: patrick
"""

def getKeggGenes(species):
    """Gets all genes, their KEGG IDs, and their descriptions from KEGG.

  Args:
      species (string): Species of interest (human = 'hsa', rat = 'rno', mouse = 'mmu')

  Returns:
      DataFrame: A DataFrame containing 3 columns; KEGG_ID, Gene_Names, Description.
  """
    from urllib.request import urlopen
    import pandas as pd
    
    # read in KEGG dataset for all human genes
    link = 'https://rest.kegg.jp/list/'+species

    fid = urlopen(link)
    genes = fid.read().decode('utf-8')
    
    human_genes_lines = genes.split('\n')

    Human_Genes_LoD = []

    for entry in human_genes_lines:
        gene_dict = {}
        if '\t' in entry:
            gene = entry.split('\t')
            KEGG_ID = gene[0]
            second_half = gene[3]
            gene_dict['KEGG_ID'] = KEGG_ID
            if ';' in second_half:
                gene_symbol_list = second_half.split(';')[0]
                gene_dict['Gene_Names'] = gene_symbol_list
                description = second_half.split(';')[1]
                gene_dict['Description'] = description.strip()
            else:
                gene_dict['Description'] = pd.NA
                gene_dict['Gene_Names'] = pd.NA
        else:
            gene_dict['KEGG_ID'] = pd.NA
            gene_dict['Description'] = pd.NA
        Human_Genes_LoD.append(gene_dict)
        
    humanGenes = pd.DataFrame(Human_Genes_LoD)
    
    return humanGenes



def getKeggPathways(species):
    """Gets all pathways from KEGG.

  Args:
      species (string): Species of interest (human = 'hsa', rat = 'rno', mouse = 'mmu')

  Returns:
      DataFrame: A DataFrame containing 2 columns; KEGG_ID, Pathways
  """
    from urllib.request import urlopen
    import pandas as pd

    link = 'https://rest.kegg.jp/link/pathway/'+species

    fid = urlopen(link)
    pathways_raw = fid.read().decode('utf-8')

    # clean data to create a data frame of human genes and their corresponding pathways
    lines = pathways_raw.split('\n')
    data = [line.split('\t') for line in lines]

    pathways = pd.DataFrame(data[0:],columns=['KEGG_ID','Pathway'])
    
    return pathways
    


def uniquePathways(KEGG_IDs):
    """Identifies all pathways linked with a list of KEGG ID genes.

  Args:
      KEGG_IDs (list): A list of KEGG_IDs

  Returns:
      list: A set list of unique pathways.
  """
    import trailblazer
    non_empty_values = [item for item in KEGG_IDs if item]
    species = non_empty_values[0].split(':')[0].strip()
    pathways = trailblazer.getKeggPathways(species)
    
    pathway_dict = {}
    
    for index1, entry1 in enumerate(non_empty_values):
        pathway_vec = []
        for index2, entry2 in enumerate(pathways['KEGG_ID']):
            if entry1 == entry2:
                pathway_vec.append(pathways['Pathway'][index2])
        pathway_dict[entry1] = pathway_vec
    
    all_pathways = list(pathway_dict.values())
    
    flattened_list = [item for sublist in all_pathways for item in sublist]
    uniquePathways = list(set(flattened_list))
    
    return uniquePathways



def pathwayVec(KEGG_IDs):
    """Creates a vector of all pathways associated with each KEGG_ID.

  Args:
      KEGG_IDs (list): A list of KEGG_IDs

  Returns:
      dictionary: A dictionary in which the key is the KEGG ID and the value is
          a vector of associated pathways.
  """
    import trailblazer
    non_empty_values = [item for item in KEGG_IDs if item]
    species = non_empty_values[0].split(':')[0].strip()
    pathways = trailblazer.getKeggPathways(species)
    
    pathway_dict = {}
    
    for index1, entry1 in enumerate(non_empty_values):
        pathway_vec = []
        for index2, entry2 in enumerate(pathways['KEGG_ID']):
            if entry1 == entry2:
                pathway_vec.append(pathways['Pathway'][index2])
        pathway_dict[entry1] = pathway_vec
    
    return pathway_dict



def fisherExact(uniquePathways, testData, Pvalues = 'P-value',Pathways = 'Pathways'):
    """Determines if pathway is significantly enriched, using Fisher Exact Test.

  Args:
      uniquePathways (list): a set list of pathways associated with data set.
      testData (list): a list of dictionaries containing key-value pairs for
          both 'ttest_pvalue' and 'Pathway_Vec'

  Returns:
      list: A list of lists containing pathways and P values.
  """
    
    import scipy.stats as sp
    
    output_LoL = []
    
    for pathway_idx, unique_pathway in enumerate(uniquePathways):
        print(pathway_idx)
        is_sig_in_pathway = 0
        is_sig_not_pathway = 0
        not_sig_in_pathway = 0
        not_sig_not_pathway = 0
    
        for biomolecule_idx, biomolecule_data in testData.iterrows():
            is_sig = 0
    
            ttest_pvalue = biomolecule_data[Pvalues]
            if ttest_pvalue < 0.05/len(testData[Pvalues]):
                is_sig = 1
    
            in_pathway = 0
    
            pathway_vec = biomolecule_data[Pathways]
            if unique_pathway in pathway_vec:
                in_pathway = 1
    
    
            if in_pathway == 1 and is_sig ==1:
                is_sig_in_pathway+=1
            elif in_pathway == 1 and is_sig == 0:
                not_sig_in_pathway +=1
            elif in_pathway == 0 and is_sig == 1:
                is_sig_not_pathway +=1
            elif in_pathway == 0 and is_sig == 0:
                not_sig_not_pathway +=1
    
        contingency_table = [[is_sig_in_pathway, is_sig_not_pathway], [not_sig_in_pathway, not_sig_not_pathway]]
        [odds_ratio, pvalue_fisher] = sp.fisher_exact(contingency_table)
    
        curr_output = [unique_pathway, pvalue_fisher]
        output_LoL.append(curr_output)
        
    return output_LoL
        


import os

def moduleMaker(nModules, 
                nNodes = 5271, 
                moduleSize = 10, 
                network = os.path.dirname(os.path.abspath(__file__))+'/inputs/G_matrix_hsa_06182020.p'):
    
    """Creates 'n' number of Modules in human network interactions, and returns them. 
    
  Args:
      nModules (int): the number of modules you want to create
      nNodes (int): the number of nodes you want to sample from
      network (pickle): a pickle file containing an adjacency matrix of interactions

  Returns:
      biomolecule_modules (list): A list of lists containing modules of genes.
  """
    
    import pickle
    import numpy as np
    from random import randint
    import os
    
    script_directory = os.path.dirname(os.path.abspath(__file__))
    
    human_kegg = pickle.load(open(script_directory+'/inputs/Human_KEGG_cytoscape_06182020.p', 'rb'))
    #human_kegg = pickle.load(open('/Users/patrick/Documents/CHBE194/Human_KEGG_cytoscape_06182020.p', 'rb'))
    #human_network = human_kegg['df_network']
    #human_attributes = human_kegg['df_attributes']
    unique_nodes = human_kegg['Unique_Nodes_Final']
    
    adj_matrix = pickle.load(open(network,'rb'))
    
    List_of_Modules = []
    

    for i in range(nModules):
        
        interval = nModules / 10
        if i % interval == 0:
            print(" "*60, end='', flush=True)
            print(f"\rmoduleMaker is {i/nModules*100:.0f}% complete",end='',flush=True)
            
        module_nodes = [randint(0,nNodes)]
        
        isloop = 1
        loop_count = 0
        
        while isloop == 1:
            rand_node = module_nodes[randint(0, len(module_nodes)-1)]
            
            # Get all the random nodes that are connected to the random node
            
            from_nodes = np.where(adj_matrix[rand_node,:]==1)
            to_nodes = np.where(adj_matrix[:,rand_node]==1)
            
            adj_nodes = []
            for k in from_nodes[0]:
                adj_nodes.append(k)
                
            for k in to_nodes[0]:
                adj_nodes.append(k)
                
            adj_nodes = list(set(adj_nodes))
            
            if len(adj_nodes) > 0:
                rand_adj_node_index = randint(0, len(adj_nodes)-1)
                curr_adj_node = adj_nodes[rand_adj_node_index]
                if curr_adj_node not in module_nodes:
                    module_nodes.append(curr_adj_node)
            else:
                pass
                
            loop_count += 1
            
            if len(module_nodes) >= moduleSize or loop_count > 5000:
                isloop = 0
                
        List_of_Modules.append(module_nodes)
        
    #match index values in each module to index of unique_nodes to identify gene or compound
    biomolecule_modules = []
    for entry in List_of_Modules:
        biomolecule_list = []
        for index, value in enumerate(entry):
            biomolecule_list.append(unique_nodes[value])
        biomolecule_modules.append(biomolecule_list)

    return biomolecule_modules



def pathfinder(biomolecule_modules, 
               biomoleculesAndPathways, 
               human_kegg = os.path.dirname(os.path.abspath(__file__))+'/inputs/Human_KEGG_cytoscape_06182020.p'):
    
    """Counts the number of significant biomolecules in randomly generated list of modules.
    
  Args:
      biomolecule_modules (list of lists): a list of randomly generated modules
      biomoleculesAndPathways (DataFrame): a DataFrame of your data containing 
          columns for 'KEGG_ID', 'Significant', 'P-value', 'Pathways'
  Returns: 
      top_network (DataFrame): a DataFrame of significant Source-Target pairs,
          to be uploaded into CytoScape
      
  """
    import pandas as pd
    import pickle
    
    List_of_Dict = []

    for idx_1, module in enumerate(biomolecule_modules):
        
        interval = len(biomolecule_modules) / 20
        if idx_1 % interval == 0:
            print(" " * 60, end='', flush=True)
            print(f"\rpathfinder is {idx_1/len(biomolecule_modules)*100:.0f}% complete",end='',flush='True')
            
        sig_dict = {}
        count = 0
        for idx_2, biomolecule in enumerate(module):
            if biomolecule in biomoleculesAndPathways['KEGG_ID'].tolist():
                index = biomoleculesAndPathways['KEGG_ID'][biomoleculesAndPathways['KEGG_ID'] == biomolecule].index[0]
                sig = biomoleculesAndPathways['Significant'][index]
                if sig == 'Yes':
                    count += 1
                else:
                    pass
        sig_dict['index'] = idx_1
        sig_dict['significant_biomolecules'] = count
        List_of_Dict.append(sig_dict)
        
    maximum_count = 0
    
    for entry in List_of_Dict:
        if entry['significant_biomolecules'] > maximum_count:
            maximum_count = entry['significant_biomolecules']
            
    top_modules = []
    
    for entry in List_of_Dict:
        if entry['significant_biomolecules'] == maximum_count:
            top_modules.append(biomolecule_modules[entry['index']])
            
    top_biomolecules_list = list(set([item for sublist in top_modules for item in sublist]))
    
    # need to create a source-target data frame to upload into cytoscape
    human_kegg = pickle.load(open(human_kegg, 'rb'))
    human_network = human_kegg['df_network']
    
    top_network_source = human_network[human_network['Source'].isin(top_biomolecules_list)]
    top_network_target = human_network[human_network['Target'].isin(top_biomolecules_list)]

    top_network = pd.concat([top_network_source, top_network_target], axis=0)
    
    return top_network



def colorPalette(color = 'yellow', nColors = 10):
        
    """Creates a color palette to use on cytoscape network. 
    
  Args:
      color (string): a string (yellow, green, blue, red)
      nColors (int): the number of colors desired in the color palette
  Returns: 
      colors (list): a list of Hex codes for specified color, in order from
          darkest to lightest
      
  """
        
    # Generate a list of hex colors between lightest and darkest yellow
    colors = ['#d3d3d3']
    
    if color == 'yellow':
        light_color = '#FFFFCC'
        dark_color = '#FFA500'
        for i in range(1, nColors):
            # Calculate intermediate color based on linear interpolation between lightest and darkest yellow
            r = int(light_color[1:3], 16) - i * (int(light_color[1:3], 16) - int(dark_color[1:3], 16)) // nColors
            g = int(light_color[3:5], 16) - i * (int(light_color[3:5], 16) - int(dark_color[3:5], 16)) // nColors
            b = int(light_color[5:], 16) - i * (int(light_color[5:], 16) - int(dark_color[5:], 16)) // nColors
        
            # Format the RGB values into a hex color string
            color_hex = f"#{r:02X}{g:02X}{b:02X}"
            
            # Append the color to the list
            colors.append(color_hex)
            
    elif color == 'green':
        light_color = '#C6F7C9'
        dark_color = '#01990A'
        for i in range(1, nColors):
            # Calculate intermediate color based on linear interpolation between lightest and darkest yellow
            r = int(light_color[1:3], 16) - i * (int(light_color[1:3], 16) - int(dark_color[1:3], 16)) // nColors
            g = int(light_color[3:5], 16) - i * (int(light_color[3:5], 16) - int(dark_color[3:5], 16)) // nColors
            b = int(light_color[5:], 16) - i * (int(light_color[5:], 16) - int(dark_color[5:], 16)) // nColors
        
            # Format the RGB values into a hex color string
            color_hex = f"#{r:02X}{g:02X}{b:02X}"
            
            # Append the color to the list
            colors.append(color_hex)
    
    elif color == 'blue':
        light_color = '#D3DFF7'
        dark_color = '#0013FD'
        for i in range(1, nColors):
            # Calculate intermediate color based on linear interpolation between lightest and darkest yellow
            r = int(light_color[1:3], 16) - i * (int(light_color[1:3], 16) - int(dark_color[1:3], 16)) // nColors
            g = int(light_color[3:5], 16) - i * (int(light_color[3:5], 16) - int(dark_color[3:5], 16)) // nColors
            b = int(light_color[5:], 16) - i * (int(light_color[5:], 16) - int(dark_color[5:], 16)) // nColors
        
            # Format the RGB values into a hex color string
            color_hex = f"#{r:02X}{g:02X}{b:02X}"
            
            # Append the color to the list
            colors.append(color_hex)
            
    elif color == 'red':
        light_color = '#FFE1E1'
        dark_color = '#FD0000'
        for i in range(1, nColors):
            # Calculate intermediate color based on linear interpolation between lightest and darkest yellow
            r = int(light_color[1:3], 16) - i * (int(light_color[1:3], 16) - int(dark_color[1:3], 16)) // nColors
            g = int(light_color[3:5], 16) - i * (int(light_color[3:5], 16) - int(dark_color[3:5], 16)) // nColors
            b = int(light_color[5:], 16) - i * (int(light_color[5:], 16) - int(dark_color[5:], 16)) // nColors
        
            # Format the RGB values into a hex color string
            color_hex = f"#{r:02X}{g:02X}{b:02X}"
            
            # Append the color to the list
            colors.append(color_hex)
            
    else:
        print('Please choose ''yellow'', ''green'', ''blue'', or ''red''.')
        
    light_to_dark = colors[1:]
    light_to_dark.reverse()
    
    colors = [colors[0]] + light_to_dark
    
    return colors



def cytoscapeAttributes(tissueTranscripts, 
                        gtex,
                        color = 'yellow', 
                        KEGG_query = None, 
                        addOn = None,
                        human_kegg = os.path.dirname(os.path.abspath(__file__))+'/inputs/Human_KEGG_cytoscape_06182020.p', 
                        human_proteins = os.path.dirname(os.path.abspath(__file__))+'/inputs/human_proteins_LoD.p',
                        ):
            
    """Generates a customized DataFrame that includes all variables required to
        stylize a cytoscape network. 
    
  Args:
      tissueTranscripts (DataFrame): a DataFrame containing columns for Transcript,
          TPM, KEGG_ID, P-value
      gtex (DataFrame): a DataFrame containing all samples and transcripts from
          the GTEx data base
      color (string): ''yellow'', ''green'', ''blue'', or ''red''
      KEGG_query (string): a keyword search to include in final attributes file
      addOn (dictionary): a dictionary of data to include in final attributes file; 
          each key should be a string, and each value should be a list. For example,
          {''Specificity'': []} contains the key 'Specificity' and value of a list
          of ints as long as the GTEx data set
      human_kegg (pickle): a pickle file containing a dictionary with df_attributes,
          df_network, and Unique_Nodes_Final
      human_proteins (pickle): a pickle file containing a list of dictionaries
          of human proteins
  Returns: 
      human_attributes_copy (DataFrame): a customized DataFrame that can be used
          to stylize data in cytoscape
      
  """
  
    import pickle
    import pandas as pd
    import trailblazer as tb
    import GTExGenie as gtx
    import math
    import warnings
    
    print(" " * 60, end='', flush=True)
    print(f"\rwriting cytoscape attributes",end='',flush='True')
    
    #gtex = pd.read_csv('/Users/patrick/Documents/CHBE194/Homework4/GTEX_Uniprot/GTEX_data_2020.gct',sep='\t',skiprows=2,index_col=0)
    human_kegg = pickle.load(open(human_kegg, 'rb'))
    human_attributes = human_kegg['df_attributes']
    human_attributes_copy = human_attributes.copy()
    human_proteins = pickle.load(open(human_proteins, 'rb'))
    human_genes = tb.getKeggGenes('hsa')
    keggMatch = gtx.KeggMatch(gtex, human_genes)
    
    for key in list(addOn.keys()):
        if 'specificity'.upper() == key.upper() and len(addOn[key]) == len(tissueTranscripts):
            tissueTranscripts['Specificity'] = addOn[key]
            specificity = addOn[key]
        else:
            pass
    
    #print(tissueTranscripts.columns)
    
    num_colors = max(specificity)-min(specificity)+1
    colors = tb.colorPalette(color = 'yellow', nColors = num_colors)
    
    for key in list(addOn.keys()):
        if 'specificity'.upper() == key.upper() and len(addOn[key]) == len(tissueTranscripts):
            tissueTranscripts['Specificity'] = addOn[key]
            
            human_attributes_copy['Specificity'] = 0
            for idx1, row_attribute in human_attributes_copy.iterrows():
                node_id = row_attribute['Node_ID']
                transcript = row_attribute.iloc[4]
                if node_id in list(tissueTranscripts['KEGG_ID']):
                    match = tissueTranscripts[tissueTranscripts['KEGG_ID']==node_id]
                    match_idx = match.index[0]
                    match = match.loc[match_idx]
                    
                    spec = match.loc['Specificity']
                    human_attributes_copy.loc[idx1,'Node_Color'] = colors[spec]
                    human_attributes_copy.loc[idx1,'Specificity'] = spec
                    
                elif transcript in list(tissueTranscripts['Transcript']):
                    match = tissueTranscripts[tissueTranscripts['Transcript']==transcript]
                    match_idx = match.index[0]
                    match = match.loc[match_idx]
                    
                    spec = match.loc['Specificity']
                    human_attributes_copy.loc[idx1,'Node_Color'] = colors[spec]
                    human_attributes_copy.loc[idx1,'Specificity'] = spec

        #make columns for each add-on gene list
        else:
            human_attributes_copy[key] = False
            for entry in addOn[key]:
                for index, row in human_attributes_copy.iterrows():
                    if entry in list(row.iloc[4:33]):
                        human_attributes_copy.loc[index,key] = True
            
    human_attributes_copy['TPM_size'] = ''
    for index, row in human_attributes_copy.iterrows():
        node_id = row['Node_ID']
        transcript = row['Display_Name_1']
        if node_id in list(tissueTranscripts['KEGG_ID']):
            match = tissueTranscripts[tissueTranscripts['KEGG_ID']==node_id]
            match_idx = match.index[0]
            match = match.loc[match_idx]
            TPM = match.loc['TPM']

            if TPM > 1000:
                human_attributes_copy.loc[index,'TPM_size'] = 'large'
            elif TPM >= 600 and TPM < 1000:
                human_attributes_copy.loc[index,'TPM_size'] = 'medium'
            elif TPM >= 200 and TPM < 600:
                human_attributes_copy.loc[index,'TPM_size'] = 'small'
            else:
                human_attributes_copy.loc[index,'TPM_size'] = 'extra small'
        elif transcript in list(tissueTranscripts['Transcript']):
            match = tissueTranscripts[tissueTranscripts['Transcript']==transcript]
            match_idx = match.index[0]
            match = match.loc[match_idx]
            TPM = match.loc['TPM']

            if TPM > 1000:
                human_attributes_copy.loc[index,'TPM_size'] = 'large'
            elif TPM >= 600 and TPM < 1000:
                human_attributes_copy.loc[index,'TPM_size'] = 'medium'
            elif TPM >= 200 and TPM < 600:
                human_attributes_copy.loc[index,'TPM_size'] = 'small'
            else:
                human_attributes_copy.loc[index,'TPM_size'] = 'extra small'
        else:
            human_attributes_copy.loc[index,'TPM_size'] = 'extra small'
            
    #get uniprot/KEGG conversion
    from urllib.request import urlopen
    
    link = 'https://rest.kegg.jp/conv/hsa/up'
    
    fid = urlopen(link)
    conversions_raw = fid.read().decode('utf-8')
    
    # clean data to create a data frame of human genes and their corresponding pathways
    lines = conversions_raw.split('\n')
    data = [line.split('\t') for line in lines]
    
    conversion = pd.DataFrame(data[0:],columns=['Uniprot_ID','KEGG_ID'])
    
    if KEGG_query == 'secreted':
        secreted_KEGG_IDs = []
        for entry in human_proteins:
            if 'secreted'.upper() in entry['Subcellular'].upper():
                for ID in entry['Accession_ids']:
                    if str('up:'+ID) in conversion['Uniprot_ID'].tolist():
                        HSA_ID = conversion.loc[conversion['Uniprot_ID'] == 'up:'+ID,'KEGG_ID'].iloc[0]
                        secreted_KEGG_IDs.append(HSA_ID)
    
        human_attributes_copy['Secreted'] = False
    
        warnings.filterwarnings('ignore', message="A value is trying to be set on a copy of a slice from a DataFrame")
    
        for index, entry in enumerate(human_attributes_copy['Node_ID']):
            if entry in secreted_KEGG_IDs:
                human_attributes_copy.loc[index,'Secreted'] = True
            else:
                human_attributes_copy.loc[index,'Secreted'] = False
    elif KEGG_query != None:
        keyword_KEGG_IDs = []
        for entry in human_proteins:
            if KEGG_query.upper() in entry['Subcellular'].upper():
                for ID in entry['Accession_ids']:
                    if str('up:'+ID) in conversion['Uniprot_ID'].tolist():
                        HSA_ID = conversion.loc[conversion['Uniprot_ID'] == 'up:'+ID,'KEGG_ID'].iloc[0]
                        keyword_KEGG_IDs.append(HSA_ID)
    
        human_attributes_copy[KEGG_query] = False
    
        for index, entry in enumerate(human_attributes_copy['Node_ID']):
            if entry in keyword_KEGG_IDs:
                human_attributes_copy[KEGG_query].loc[index] = True
            else:
                human_attributes_copy[KEGG_query].loc[index] = False
    
    
    return human_attributes_copy
    


def cytoscape(unique_tissues, 
              gtex, 
              attributes, 
              human_genes = None, 
              specificity = None, 
              KEGG_query = None,
              addOn = None, 
              nModules = 100000, 
              color = 'yellow',
              version = ''):
    
    """Generates and saves a network file that can be uploaded to cytoscape, 
        along with its associated attributes file. 
    
  Args:
      unique_tissues (list): a set list of unique tissues that will be analyzed
      gtex (DataFrame): all GTEx transcriptomic data
      attributes (DataFrame): associated attributes data used to stylize the network
      human_genes (DataFrame): a DataFrame of human genes containing KEGG_ID,
          Gene_Names, and Description
      specificity (Series): a pandas Series object containing gene specificity
      KEGG_query (string): a keyword search to include in final attributes file
      addOn (dictionary): a dictionary of data to include in final attributes file; 
          each key should be a string, and each value should be a list. For example,
          {''Specificity'': []} contains the key 'Specificity' and value of a list
          of ints as long as the GTEx data set
      nModules (int): number of modules desired, defaults to 100k
      color (string): desired color palette of network, defaults to yellow
      version (string): a modifier to the final file names saved in 'outputs'
    """
    
    import GTExGenie as gtx
    import trailblazer as tb
    
    if human_genes == None:
        human_genes = tb.getKeggGenes('hsa')
        
    for tissue in unique_tissues:
        idx = unique_tissues.index(tissue)
        print(f"looping through unique tissues: {idx}/{len(unique_tissues)} complete",flush=True)
        
        TPM_tissue = gtx.tissueFilter(tissue, gtex, attributes)
        if TPM_tissue.shape[1] == 0:
            continue
        
        enriched_transcripts = gtx.enrichedTranscripts(gtex, TPM_tissue)
        tissue_kegg_match = gtx.KeggMatch(enriched_transcripts['Transcript'], human_genes)
        
        unique_tissue_pathways = tb.uniquePathways(list(tissue_kegg_match.values()))
        pathway_vec = tb.pathwayVec(list(tissue_kegg_match.values()))
        
        #make a data frame containing Transcript Name, KEGG_ID, P-value, and Pathway Vec
        tissue_KEGG_IDs = []
        
        for transcript in enriched_transcripts['Transcript']:
            tissue_KEGG_IDs.append(tissue_kegg_match[transcript])
            
        #find pathways for each of those transcripts in tissue_KEGG_IDs
        pathway_vecs = []
        
        for transcript in tissue_KEGG_IDs:
            if transcript == '':
                pathway_vecs.append([])
            else:
                pathway_vecs.append(pathway_vec[transcript])
        
        #make an add-on column with genes of interest
        if addOn != None:
            keys = list(addOn.keys())

            for key in keys:
                if len(addOn[key]) == len(enriched_transcripts):
                    enriched_transcripts[key] = addOn[key]
                else:
                    column = [False] * len(enriched_transcripts)
                    match = False
                    for idx1, gtex_gene in enumerate(gtex.Description):
                        for idx2, entry in enumerate(addOn[key]):
                            if entry.upper() in gtex_gene.upper():
                                column[idx1] = True
                                match = True
                                break
                        if match != True:
                            column[idx1] = False
                    enriched_transcripts[key] = column
                                
        #add on columns to data frame
        #if specificity != None:
            #enriched_transcripts['Specificity'] = specificity
        enriched_transcripts['KEGG_ID'] = tissue_KEGG_IDs
        enriched_transcripts['Pathways'] = pathway_vecs
        
        #generate a list of random modules, based on human network of interactions
        modules = tb.moduleMaker(nModules=nModules)

        #cull a network of source-target pairs, specific for adipose tissue, to upload into cytoscape
        tissue_network = tb.pathfinder(modules,enriched_transcripts)
        
        #make a dataframe that can be uploaded to cytoscape for network analysis
        cytoscape_attributes = tb.cytoscapeAttributes(enriched_transcripts,gtex,addOn=addOn, KEGG_query=KEGG_query)
        
        
        tissue_file_name = tissue.replace(' ','')
        tissue_network.to_csv(os.path.dirname(os.path.abspath(__file__))+'/outputs/'+tissue_file_name+version+'_network.csv', index=False)
        cytoscape_attributes.to_csv(os.path.dirname(os.path.abspath(__file__))+'/outputs/'+tissue_file_name+version+'_attributes.csv')
