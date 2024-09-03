#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sccoda


# In[1]:


# Setup
import warnings
warnings.filterwarnings("ignore")

import pandas as pd
import pickle as pkl
import matplotlib.pyplot as plt

from sccoda.util import comp_ana as mod
from sccoda.util import cell_composition_data as dat
from sccoda.util import data_visualization as viz

import sccoda.datasets as scd


# In[2]:


cell_counts = scd.haber()


# cellcomp=pd.read_csv("C:/Users/cesar/OneDrive - CRG - Centre de Regulacio Genomica/Collaborations/Alvaro/Alvaro snRNAseq/Clara_Dierssen/inprogress/glia_cell_composition.csv", header=y)

# In[4]:


cellcomp=pd.read_csv("cell_composition.csv", header=0)
cellcomp


# In[5]:


#Convert data to anndata object
data_all = dat.from_pandas(cellcomp, covariate_columns=["Var1"])
# Extract condition from mouse name and add it as an extra column to the covariates
condition=data_all.obs["Var1"].str.replace(r"_[0-9]", "", regex=True)
condition2=condition.str.replace(r"[0-9]", "", regex=True)
condition2
data_all.obs["Condition"] = condition2
print(data_all)


# In[6]:


data_all.to_df()


# In[7]:


viz.boxplots(data_all, feature_name="Condition")


# In[8]:


# Stacked barplot for the levels of "Condition"
ax=viz.stacked_barplot(data_all, feature_name="Condition")
# Get the figure from the AxesSubplot
fig = ax.get_figure()

# Save the plot as a PDF file
fig.savefig("composition.pdf")


# In[9]:


model_data_all = mod.CompositionalAnalysis(data_all, formula="Condition", reference_cell_type="oligop")


# In[10]:


# Run MCMC
sim_results = model_data_all.sample_hmc()


# In[11]:


sim_results.summary()
print(sim_results.credible_effects())
#Only DAM is different


# In[ ]:




