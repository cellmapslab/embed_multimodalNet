#!/usr/bin/env python
# coding: utf-8

# # Get top 500 most similar nodes

# In[3]:


import argparse
from pathlib import Path


# In[ ]:


parser = argparse.ArgumentParser()
parser.add_argument('--proj_folder')
parser.add_argument('--traitlist')
parser.add_argument('--output')
args = parser.parse_args()


# In[4]:


import networkx as nx
from node2vec import Node2Vec
import os
import pickle
import random
from itertools import chain
import re
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import gensim
from gensim.models import Word2Vec
from random import sample
import seaborn as sns
import numpy as np


# In[5]:


def save_pickle(obj, project_folder, prefix):
    fname = os.path.join(project_folder, '%s.pkl' % prefix)
    with open(fname, 'wb') as fh:
        pickle.dump(obj, fh)

def load_pickle(project_folder, prefix):
    fname = os.path.join(project_folder, '%s.pkl' % prefix)
    with open(fname, 'rb') as fh:
        return pickle.load(fh)


# ## Functions



def get_sp(mg, source, target):
    try:
        sp = nx.shortest_path_length(mg, source, target)
        return sp
    except nx.exception.NetworkXNoPath:
        return 'Inf'


# In[47]:


def get_sim_scores(node, nv, mg):
    
    node = "GO:"+node
    data_list=nv.similar_by_word(node, topn=1000)
    data=pd.DataFrame(data_list, columns=['node2', 'sim'])
    data['node1']=node
    
    neigh_list=[]
    for i in range(data.shape[0]):
        node2=data.iloc[i]['node2']
        neigh_list.append(get_sp(mg, source=node, target=node2))
        
    data['neighbour'] = neigh_list
    #data = data[['node1', 'node2', 'sim', 'neighbour', 'most_least']]
    data = data[['node1', 'node2', 'sim', 'neighbour']]

    data = data.replace("GO:", "", regex=True) 

    return data.values.tolist() 
    


# In[30]:


def sim_score_allruns(trait, projectfolder=args.proj_folder):
    myfiles= ['deepwalk_node_vectors_'+str(i) for i in range(1,101)]
    i=1
    
    df2=[]
    
    for file in myfiles:
        df=[]
        nv=load_pickle(project_folder=projectfolder,prefix=file)
        #nv=nv=gensim.models.KeyedVectors.load_word2vec_format(os.path.join(projectfolder, '%s.txt' % file))
        df=get_sim_scores(node=trait, nv=nv, mg=mg)
        df=[x + [i] for x in df]
       # print(pd.DataFrame(df))
        df2.append(df)
        i += 1
        
    df3=list(chain.from_iterable(df2))
    res = pd.DataFrame(df3, columns=['node1', 'node2', 'sim', 'neighbour','run'])
    return res
    


# ## Load data

# In[ ]:


mg=load_pickle(project_folder=args.proj_folder,prefix='multi_graph')


# In[48]:


traits=pd.read_csv(args.traitlist, sep="\t")


# ## Execute

# In[ ]:


#df = pd.DataFrame()
non_present = 0
for i in traits['gene']:
    print(i)
    try:
        #print(i)
        df1=sim_score_allruns(i)
        #df=df.append(df1)
        df1.to_csv(args.output, header=None, mode="a")
    except KeyError:
        non_present+= 1
        continue
print("Number of unpresent covid genes: ", non_present)


