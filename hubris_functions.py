import os
import matplotlib.pyplot as plt
import pandas as pd
import networkx as nx
import numpy as np
import urllib
from biomart import BiomartServer
import pickle
import glob
import random
import biorosetta as br
import tqdm
import urllib.request
import gzip
import shutil

def hubris_download_databases(to_download, df_databases):
    
    '''
    Downloads the specified databases and does some pre-processing for the further HUBRIS steps
    Inputs: to_download (list) - list of database names to be re-downloaded
            df_databases - database configuration file

    The files will be downloaded to the file paths specified in the df_databases (db_local_file)
    '''



    for db_name in to_download:
        print(db_name)
        db_row = df_databases[df_databases.db_name==db_name]

        print('Downloading %s from %s to %s' %(db_name,db_row['link'].iloc[0],db_row['db_local_file'].iloc[0]))
        _,response=urllib.request.urlretrieve(db_row['link'].iloc[0], db_row['db_local_file'].iloc[0])

        #print(response.as_string())
        
        if db_name == 'StringDB':
            print('unzipping')
            with gzip.open(db_row['db_local_file'].iloc[0], 'rb') as file_in:
                with open(db_row['db_local_file'].iloc[0]+'_temp', 'wb') as file_out:
                    shutil.copyfileobj(file_in, file_out)
            print('Unzipping_done')
            os.remove(db_row['db_local_file'].iloc[0])
            os.rename(db_row['db_local_file'].iloc[0]+'_temp', db_row['db_local_file'].iloc[0])
            
        if db_name == 'BioGrid':
            
            from zipfile import ZipFile
            print('unzipping')
            with ZipFile(db_row['db_local_file'].iloc[0], 'r') as file_in:
                temp_name=file_in.namelist()[0]
                file_in.extractall(db_row['db_local_file'].iloc[0].split('/')[0])
            print('Unzipping_done')
            os.remove(db_row['db_local_file'].iloc[0])
            os.rename(db_row['db_local_file'].iloc[0].split('/')[0]+'/'+temp_name, db_row['db_local_file'].iloc[0])
            print('processing columns')
            df=pd.read_csv(db_row['db_local_file'].iloc[0],sep=db_row['sep'].iloc[0])
            df=df[(df['Taxid Interactor A']=='taxid:9606') & (df['Taxid Interactor B']=='taxid:9606')]
            df['source']=df.apply(lambda row: row['#ID Interactor A'].split(':')[-1],axis=1)
            df['target']=df.apply(lambda row: row['ID Interactor B'].split(':')[-1],axis=1) 
            df.to_csv(db_row['db_local_file'].iloc[0],sep='\t')
        
        if db_name=='Reactome':
            
            df=pd.read_csv(db_row['db_local_file'].iloc[0],sep=db_row['sep'].iloc[0])
            df['source']=df.apply(lambda row: row['# Interactor 1 uniprot id'].split(':')[-1],axis=1)
            df['target']=df.apply(lambda row: row['Interactor 2 uniprot id'].split(':')[-1],axis=1) 
            df.to_csv(db_row['db_local_file'].iloc[0],sep='\t')
            
        if db_name=='NCBI':
            print('unzipping')
            with gzip.open(db_row['db_local_file'].iloc[0], 'rb') as file_in:
                with open(db_row['db_local_file'].iloc[0]+'_temp', 'wb') as file_out:
                    shutil.copyfileobj(file_in, file_out)
            print('Unzipping_done')
            os.remove(db_row['db_local_file'].iloc[0])
            os.rename(db_row['db_local_file'].iloc[0]+'_temp', db_row['db_local_file'].iloc[0])
            
    return 1
    
    
    
def create_graphs_with_consensus_ids(df_databases,to_merge,df_hgnc,hgnc_column_dict,consensus_id_type='entr',biorosetta_data_path=''):
    
    '''
    This function iterates through the rows of the database configuration table (df_databases), converts the ids
    to the consensus_id_type using Biorosetta (https://pypi.org/project/biorosetta/) and a local file of the
    hgnc mapping (to download a more recent version visit: https://www.genenames.org/download/statistics-and-files/)
    and greates the undirected PPI graph.
    
    Inputs: df_databases (pandas DataFrame) - db configuration table
            consensus_id_type (str; default='entr') - consensus protein id to which all the db nodes are converted
            ID types:
            'ensg' = Ensembl gene ID
            'entr' = NCBI gene ID (entrezgene)
            'symb' = Gene symbol
            'ensp' = Ensembl protein ID
            'hgnc' = HGNC ID
            'unipr' = Uniprot ID
            to_merge: (list) database names to convert and prepare for merging
            df_hgnc: (pandas DataFrame) the hgnc mapping table
            hgnc_column_dict: (dict) a mapping between the column names of the hgnc mapping table and the ID types listed above
    Returns: a dictionary with database names as keys and nx graphs (with converted ids) as values
    '''

    edge_lists={}
    graphs={}
    for i,row in list(df_databases.iterrows()):

        if row.db_name not in to_merge:
            print(row.db_name, 'not in to_merge')
            continue
        print('Processing',row.db_name)
        header=row.header
        if pd.notna(row.header):
            header=row.header.split(', ')
            df=pd.read_csv(row.db_local_file,sep=row.sep,names=header,engine='python')
        else:
            df=pd.read_csv(row.db_local_file,sep=row.sep,engine='python')

        source_column=row.source_column_name
        target_column=row.target_column_name

        #db specific alterations
        if row.db_name=='StringDB':
            df=df[df.experimental>0]
            df[source_column]=[i.split('.')[1] for i in df[source_column]]
            df[target_column]=[i.split('.')[1] for i in df[target_column]]
        if row.db_name=='HIPPIE':
            df[source_column]=pd.to_numeric(df[source_column], errors='coerce')
            df[target_column]=pd.to_numeric(df[target_column], errors='coerce')
            df = df[~df[[source_column,target_column]].isna().any(axis=1)] 
            df[source_column]=df[source_column].astype(int)
            df[target_column]=df[target_column].astype(int)
        if row.db_name=='NCBI':
            df=df[df['#tax_id']==9606]

        #df=df.dropna()

        edge_lists[row.db_name]=set(zip(df[source_column].astype(str),df[target_column].astype(str))).union(zip(df[target_column].astype(str),df[source_column].astype(str)))

        print('id_type:',row.id_type_short)

        if row.id_type_short!=consensus_id_type:

            print('id type of db not identical with consensus id_type, beggining translation')
            node_list=list(set(df[source_column]).union(set(df[target_column])))
            if row.id_conversion == 'biorosetta':
                #idmap = br.IDMapper([br.EnsemblBiomartMapper(),br.HGNCBiomartMapper(),br.MyGeneMapper()]) # Multiple sources
                idmap = br.IDMapper([br.EnsemblBiomartMapper(data_path=biorosetta_data_path+'ensembl.tsv'),                         br.HGNCBiomartMapper(data_path=biorosetta_data_path+'hgnc.tsv')]) 
                translation_list=idmap.convert(node_list,row.id_type_short,consensus_id_type, multi_hits='all')
                translation_dict=dict(zip(node_list,translation_list))
                for k,v in translation_dict.items():
                    if v=='N/A':
                        translation_dict[k]=k+'_'+row.id_type_short
            elif row.id_conversion == 'hgnc':
                df_hgnc_selection=df_hgnc[df_hgnc[hgnc_column_dict[row.id_type_short]].isin(node_list)]
                translation_dict=dict(zip(df_hgnc_selection[hgnc_column_dict[row.id_type_short]],
                                          df_hgnc_selection[hgnc_column_dict[consensus_id_type]].astype(int).astype(str)))
                                          # df_hgnc_selection[hgnc_column_dict[consensus_id_type]].astype(str)))
                for i in set(node_list)-set(translation_dict):
                    translation_dict[i]=i+'_'+row.id_type_short
            print('translation done')
            print('creating Graph')
            G=nx.Graph()
            G.add_edges_from(edge_lists[row.db_name], db=row.db_name)
            print('relabeling nodes')    
            G=relabel_nodes_with_preserving_attributes(G,translation_dict)

        else:

            G=nx.Graph()
            G.add_edges_from(edge_lists[row.db_name], db=row.db_name)

        graphs[row.db_name]=G
    return graphs
    
    
def filter_hubris_based_on_cell_type_specific_gene_list(G_hubris,cell_line,cell_lines, keep_nodes, gene_list_file_path = 'RNASeq_lists/%s_RNASeq_genes_expressed.csv',biorosetta_data_path='biorosetta_data/data/'):
    
    '''
    This function returns a subgraph of the HUBRIS network (G_hubris) based on a pre-calculated list of 
    genes which is generated based on RNASeq data separately for each cell line (gene_list_file_path)
    
    The RNA Seq gene lists contain the gene ids in ensg format so we also do a quick coversion to gene symbols 
    which is the default ID type of HUBRIS
    
    Inputs: G_hubris (networkx Graph)
            cell_line (str): the cell line analyized. "union" and "intersection" are also accepted parameters
            cell_lines (list): in case the "cell_line" option is "union"/"intersection" the list of cell lines
            keep_nodes (list): list of nodes to be kept in the network regardless of expression
            gene_list_file_path (str): path of the gene lists to be used for filtering. The str should contain
                a single '%s' formatting placement for the cell line name. E.g. gene_list_%s where the %s will 
                be automatically replaced by the cell lines name in the function.  
                
    Returns: G_hubris (networkx Graph) a subgraph of the input with only the nodes expressed + keep_nodes kept
                
    '''
    
    if cell_line=='intersection':
        df_rnaseq_nodes=set(pd.read_csv(gene_list_file_path%cell_lines[0])['geneID'])
        for cl in cell_lines[1:]:
            df_rnaseq_nodes=df_rnaseq_nodes.intersection(set(pd.read_csv(gene_list_file_path%cl)['geneID']))
    elif cell_line=='union':
        df_rnaseq_nodes=set(pd.read_csv(gene_list_file_path%cell_lines[0])['geneID'])
        for cl in cell_lines[1:]:
            df_rnaseq_nodes=df_rnaseq_nodes.union(set(pd.read_csv(gene_list_file_path%cl)['geneID']))

    else: 
        df_rnaseq_nodes=list(pd.read_csv(gene_list_file_path%cell_line)['geneID'])
            
    import biorosetta as br
    #idmap = br.IDMapper([br.EnsemblBiomartMapper(),br.HGNCBiomartMapper(),br.MyGeneMapper()]) # Multiple sources
    idmap = br.IDMapper([br.EnsemblBiomartMapper(data_path=biorosetta_data_path+'ensembl.tsv'),                         br.HGNCBiomartMapper(data_path=biorosetta_data_path+'hgnc.tsv')]) 
    translation_list=idmap.convert(list(df_rnaseq_nodes),'ensg','symb', multi_hits='all')
    translation_list+=keep_nodes
  
    return G_hubris.subgraph(translation_list)
    
def str_to_bool(s):
    if s == 'True' or s == '1':
         return True
    elif s == 'False' or s == '0':
         return False
    else:
         raise ValueError('Incorrect input: %s'%s)
         
         
def induced_graph(G,source,target_set):
    
    '''Creates a subgraph from the original G graph and all the nodes that are part of a shortest path between 
    the node-set included in source and target. All equivalent shortest paths are included. 
    
    Returns: G_sub - a subgraph made of the node set defined above generated from the original G
    sp_count_and_len - a 2D array containing the length and number of shortest paths between source and target'''

    sp_count_and_len={}
    subgraph=set()
    sp_count_and_len=np.zeros((2,len(target_set)))

    for i,protein in tqdm.tqdm(enumerate(target_set)):
    #print(i)

        if protein not in G.nodes():
            print(protein, ' not in network')
            continue

        sps=list(nx.shortest_paths.all_shortest_paths(G,source,protein))
        sp_count_and_len[0,i]=len(sps)
        sp_count_and_len[1,i]=len(sps[0])-1

        for sp in sps:
            subgraph=subgraph.union(set(sp))

    return G.subgraph(subgraph), sp_count_and_len   
    
def deep_merge_attributes(a,b):
    '''
    Merges two dictionaries with the following additional operations:
    1. the keys of both dictionaries are preserved 
    2. If the keys intersect the values of the interesecting keys are united in a set
    3. Values of non-intersecting keys are kept in their original form
    
    E.g. :
    
    a={'x':1, 'y':[3,4], 'w':'s'}
    b={'x':2, 'y':[5,6], 'z':7}
    
    Resulting merge should be:
    {'x':{1,2}, 'y':{3,4,5,6}, 'w':'s', 'z':7}
    
    '''

    shared_attribues = set(a).intersection(set(b))

    merged_dict={}

    for k in set(a)-shared_attribues:
        merged_dict[k] = a[k]
    for k in set(b)-shared_attribues:
        merged_dict[k] = b[k]
    for k in shared_attribues:
        if hasattr(a[k],'__iter__'):
            to_merge_a = set(a[k])
        else: 
            to_merge_a = set([a[k]])
        if hasattr(b[k],'__iter__'):
            to_merge_b = set(b[k])
        else: 
            to_merge_b = set([b[k]])
        merged_dict[k] = to_merge_a.union(to_merge_b)
    return merged_dict


def relabel_nodes_with_preserving_attributes(G, translation_dict):
    
    '''
    Relabels nodes of G ased on the translation_dict, however if due to non-unique translations there are
    edge duplicates created, the attributes of the edges are merged (see deep_merge_attributes(a,b))
    
    '''
    
    
    edge_attr_mismatches = []
    G_r = G.__class__()

    for e in list(G.edges()):
        new_edge = (translation_dict[e[0]],translation_dict[e[1]])
        edge_attributes = G.edges[e]
        #it's possible the edge has already been created in the new graph due to non-uniqe translation
        if G_r.has_edge(*new_edge): 
            equivalent_edge_attributes = G_r.edges[new_edge]
            if equivalent_edge_attributes != edge_attributes:
                merged_attributes = deep_merge_attributes(equivalent_edge_attributes,edge_attributes)
                
                nx.set_edge_attributes(G_r, {new_edge:merged_attributes})

        else:
            G_r.add_edge(*new_edge,**edge_attributes)
    
    return G_r    
 

