#!/usr/bin/env python
# coding: utf-8

# # Embedding pkl to txt

# In[ ]:


import argparse
from pathlib import Path


# In[ ]:


parser = argparse.ArgumentParser()
parser.add_argument('--proj_folder')
parser.add_argument('--out_folder')
args = parser.parse_args()


# In[2]:


import networkx as nx
from node2vec import Node2Vec
import matplotlib 
import os
import pickle
import random
import re
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from gensim.models import Word2Vec


# In[3]:


def save_pickle(obj, project_folder, prefix):
    fname = os.path.join(project_folder, '%s.pkl' % prefix)
    with open(fname, 'wb') as fh:
        pickle.dump(obj, fh)

def load_pickle(project_folder, prefix):
    fname = os.path.join(project_folder, '%s.pkl' % prefix)
    with open(fname, 'rb') as fh:
        return pickle.load(fh)


# In[7]:


def save_node_emb(projectfolder, output_folder):
    myfiles= ['deepwalk_node_vectors_'+str(i) for i in range(1,101)]
    
    i=1
    for file in myfiles:
        nv=load_pickle(project_folder=projectfolder,prefix=file)
        nv.save_word2vec_format(os.path.join(output_folder+ '/' +file +'.txt'), binary=False)
        i=i+1
        print(os.path.join(output_folder+ '/' +file +'.txt'))
        print(i)



# In[8]:


save_node_emb(projectfolder=args.proj_folder, output_folder=args.out_folder)


# In[33]:


i=1
os.path.join('../data/embedding_vectors/2021-10-05_deepwalk_node_vectors_' + str(i) +'.txt')


# In[ ]:




