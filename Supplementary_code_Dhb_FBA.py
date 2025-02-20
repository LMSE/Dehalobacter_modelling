#!/usr/bin/env python
# coding: utf-8

# # Supplementary code for FBA simulations
# 
# This code is supplemental to "Modeling reveals metabolic basis of competition among Dehalobacter strains during tandem CF and DCM metabolism"; it is used to run the FBA simulations discussed in the paper, and results in generation of Supplementary Tables S10 and S11.

# ## 1. Import the modules and models

# In[1]:


#import necessary modules
import csv 
import cobra.flux_analysis
import cobra.test
#import cameo
#import plotly
#import scipy
#import os
#from os.path import join
from cobra import Model, Reaction, Metabolite
import pandas as pd


# In[16]:


#import the models
model1 = cobra.io.read_sbml_model('iOB638.xml')
model2 = cobra.io.read_sbml_model('iOB649.xml')

models = [model1,model2]


# In[19]:


#list the amino acid exchange reactions for later
aa_list = ["EX_ala__L_e", "EX_arg__L_e", "EX_asn__L_e","EX_asp__L_e","EX_gln__L_e",
               "EX_glu__L_e","EX_gly_e","EX_his__L_e","EX_ile__L_e","EX_leu__L_e",
               "EX_lys__L_e", "EX_met__L_e","EX_pro__L_e","EX_thr__L_e","EX_trp__L_e",
               "EX_tyr__L_e","EX_val__L_e","EX_phe__L_e", "EX_cys__L_e", "EX_ser__L_e"]

for model in models:
    for reaction_id in aa_list:
        reaction = model.reactions.get_by_id(reaction_id)
        reaction.bounds = (0,1000)


# ## 2. Run FBA simulations

# ### 2.1. Mode 1 (H2/CF) simulations

# #### Mode 1: as is

# In[20]:


#change uptake constraints and hydrogenase inhibition
for model in models:
    #acceptors/donors
    model.reactions.EX_dcm_e.bounds = (-0,1000)
    model.reactions.EX_cf_e.bounds = (-10,1000)
    model.reactions.EX_h2_e.bounds = (-10,1000)

    #carbon
    model.reactions.EX_co2_e.bounds = (-10,1000)
    model.reactions.EX_ac_e.bounds = (-10,1000)
    
    #hydrogenases
    model.reactions.HYDA_ech.bounds = (-0,0)
    model.reactions.HYD_NADH.bounds = (-0,0)

#simulate
sol_1 = model1.optimize()
sol_2 = model2.optimize()


# #### Mode 1A: + amino acids

# In[23]:


#add amino acid uptake
for model in models:
    for reaction_id in aa_list:
        reaction = model.reactions.get_by_id(reaction_id)
        reaction.bounds = (-1,1000)

#simulate
sol_3 = model1.optimize()
sol_4 = model2.optimize()


# #### Mode 1B: + NADH hydrogenase

# In[27]:


for model in models:
    #turn off aa uptake
    for reaction_id in aa_list:
        reaction = model.reactions.get_by_id(reaction_id)
        reaction.bounds = (0,1000)
    #hydrogenases
    model.reactions.HYDA_ech.bounds = (-0,0)
    model.reactions.HYD_NADH.bounds = (-1000,1000)

#simulate
sol_5 = model1.optimize()
sol_6 = model2.optimize()


# #### Mode 1C: + amino acids and NADH hydrogenase

# In[31]:


for model in models:
    #aa uptake
    for reaction_id in aa_list:
        reaction = model.reactions.get_by_id(reaction_id)
        reaction.bounds = (-1,1000)
    #hydrogenases
    model.reactions.HYDA_ech.bounds = (-0,0)
    model.reactions.HYD_NADH.bounds = (-1000,1000)

#simulate
sol_7 = model1.optimize()
sol_8 = model2.optimize()


# ### 2.2. Mode 1 (DCM/CF) simulations

# #### Mode 2: as is

# In[35]:


for model in models:
    #acceptors/donors
    model.reactions.EX_dcm_e.bounds = (-0,1000)
    model.reactions.EX_cf_e.bounds = (-10,1000)
    model.reactions.EX_h2_e.bounds = (-0,1000)
    
    #aa uptake
    for reaction_id in aa_list:
        reaction = model.reactions.get_by_id(reaction_id)
        reaction.bounds = (0,1000)

    #carbon
    model.reactions.EX_co2_e.bounds = (-10,1000)
    model.reactions.EX_ac_e.bounds = (-0,1000)
    
    #hydrogenases
    model.reactions.HYDA_ech.bounds = (-100,0)
    model.reactions.HYD_NADH.bounds = (-0,0)
    
#simulate
sol_9 = model1.optimize()
sol_10 = model2.optimize()


# #### Mode 2A: + amino acids

# In[41]:


for model in models:
    #acceptors/donors
    model.reactions.EX_dcm_e.bounds = (-0,1000)
    model.reactions.EX_cf_e.bounds = (-10,1000)
    model.reactions.EX_h2_e.bounds = (-0,1000)
    
    #aa uptake
    for reaction_id in aa_list:
        reaction = model.reactions.get_by_id(reaction_id)
        reaction.bounds = (-1,1000)

    #carbon
    model.reactions.EX_co2_e.bounds = (-10,1000)
    model.reactions.EX_ac_e.bounds = (-0,1000)
    
    #hydrogenases
    model.reactions.HYDA_ech.bounds = (-100,0)
    model.reactions.HYD_NADH.bounds = (-0,0)
    
#simulate
sol_11 = model1.optimize()
sol_12 = model2.optimize()


# #### Mode 2B: + NADH hydrogenase

# In[46]:


for model in models:
    #acceptors/donors
    model.reactions.EX_dcm_e.bounds = (-0,1000)
    model.reactions.EX_cf_e.bounds = (-10,1000)
    model.reactions.EX_h2_e.bounds = (-0,1000)
    
    #aa uptake
    for reaction_id in aa_list:
        reaction = model.reactions.get_by_id(reaction_id)
        reaction.bounds = (0,1000)

    #hydrogenases
    model.reactions.HYDA_ech.bounds = (-100,0)
    model.reactions.HYD_NADH.bounds = (-100,1000)
    
#simulate
sol_13 = model1.optimize()
sol_14 = model2.optimize()


# #### Mode 2C: + amino acids and NADH hydrogenase

# In[49]:


for model in models:
    #acceptors/donors
    model.reactions.EX_dcm_e.bounds = (-0,1000)
    model.reactions.EX_cf_e.bounds = (-10,1000)
    model.reactions.EX_h2_e.bounds = (-0,1000)
    
    #aa uptake
    for reaction_id in aa_list:
        reaction = model.reactions.get_by_id(reaction_id)
        reaction.bounds = (-1,1000)

    #hydrogenases
    model.reactions.HYDA_ech.bounds = (-100,0)
    model.reactions.HYD_NADH.bounds = (-100,1000)
    
#simulate
sol_15 = model1.optimize()
sol_16 = model2.optimize()


# ### Mode 2 (DCM/CF), with HYDFDN2 inhibited to avoid cycling

# #### Mode 2A: + amino acids

# In[76]:


for model in models:
    #acceptors/donors
    model.reactions.EX_dcm_e.bounds = (-0,1000)
    model.reactions.EX_cf_e.bounds = (-10,1000)
    model.reactions.EX_h2_e.bounds = (-0,1000)
    
    #aa uptake
    for reaction_id in aa_list:
        reaction = model.reactions.get_by_id(reaction_id)
        reaction.bounds = (-1,1000)

    #carbon
    model.reactions.EX_co2_e.bounds = (-10,1000)
    model.reactions.EX_ac_e.bounds = (-0,1000)
    
    #hydrogenases
    model.reactions.HYDA_ech.bounds = (-100,0)
    model.reactions.HYD_NADH.bounds = (-0,0)
    model.reactions.HYDFDN2.bounds = (-100,0)
    
#simulate
sol_17 = model1.optimize()
sol_18 = model2.optimize()


# #### Mode 2B: + NADH hydrogenase

# In[79]:


for model in models:
    #acceptors/donors
    model.reactions.EX_dcm_e.bounds = (-0,1000)
    model.reactions.EX_cf_e.bounds = (-10,1000)
    model.reactions.EX_h2_e.bounds = (-0,1000)
    
    #aa uptake
    for reaction_id in aa_list:
        reaction = model.reactions.get_by_id(reaction_id)
        reaction.bounds = (0,1000)

    #hydrogenases
    model.reactions.HYDA_ech.bounds = (-100,0)
    model.reactions.HYD_NADH.bounds = (-100,1000)
    model.reactions.HYDFDN2.bounds = (-100,0)
    
#simulate
sol_19 = model1.optimize()
sol_20 = model2.optimize()


# #### Mode 2C: + amino acids and NADH hydrogenase

# In[120]:


for model in models:
    #acceptors/donors
    model.reactions.EX_dcm_e.bounds = (-0,1000)
    model.reactions.EX_cf_e.bounds = (-10,1000)
    model.reactions.EX_h2_e.bounds = (-0,1000)
    
    #aa uptake
    for reaction_id in aa_list:
        reaction = model.reactions.get_by_id(reaction_id)
        reaction.bounds = (-1,1000)

    #hydrogenases
    model.reactions.HYDA_ech.bounds = (-100,0)
    model.reactions.HYD_NADH.bounds = (-100,100)
    model.reactions.HYDFDN2.bounds = (-100,0)
    
#simulate
sol_21 = model1.optimize()
sol_22 = model2.optimize()


# In[121]:


model2.summary()


# ### 2.3. Mode 3 (DCM/H+) simulations

# #### Mode 3: as is

# In[110]:


#DCM alone
for model in models:
    #acceptors/donors
    model.reactions.EX_dcm_e.bounds = (-10,1000)
    model.reactions.EX_cf_e.bounds = (-0,1000)
    model.reactions.EX_h2_e.bounds = (-0,1000)

    #carbon
    model.reactions.EX_co2_e.bounds = (-0,1000)
    model.reactions.EX_ac_e.bounds = (-0,1000)
    
    #aa uptake
    for reaction_id in aa_list:
        reaction = model.reactions.get_by_id(reaction_id)
        reaction.bounds = (0,1000)
        
    #hydrogenases
    model.reactions.HYDA_ech.bounds = (-100,0)
    model.reactions.HYD_NADH.bounds = (-0,0)
    model.reactions.HYDFDN2.bounds = (-100,100)
    
#simulate
sol_23 = model1.optimize()
sol_24 = model2.optimize()


# #### Mode 3A: + amino acids (10 mM, 5 mM, 3 mM, 1 mM)

# In[116]:


#test 10mM amino acids
for model in models:
    #acceptors/donors
    model.reactions.EX_dcm_e.bounds = (-10,1000)
    model.reactions.EX_cf_e.bounds = (-0,1000)
    model.reactions.EX_h2_e.bounds = (-0,1000)
    
    #aa uptake
    for reaction_id in aa_list:
        reaction = model.reactions.get_by_id(reaction_id)
        reaction.bounds = (-10,1000)
        
    #hydrogenases
    model.reactions.HYDA_ech.bounds = (-100,0)
    model.reactions.HYD_NADH.bounds = (-0,0)
    model.reactions.HYDFDN2.bounds = (-100,100)
    
#simulate
sol_25 = model1.optimize()
sol_26 = model2.optimize()


#test 5mM amino acids
for model in models:
    #acceptors/donors
    model.reactions.EX_dcm_e.bounds = (-10,1000)
    model.reactions.EX_cf_e.bounds = (-0,1000)
    model.reactions.EX_h2_e.bounds = (-0,1000)
    
    #aa uptake
    for reaction_id in aa_list:
        reaction = model.reactions.get_by_id(reaction_id)
        reaction.bounds = (-5,1000)

#simulate
sol_27 = model1.optimize()
sol_28 = model2.optimize()


#test 3mM amino acids
for model in models:
    #acceptors/donors
    model.reactions.EX_dcm_e.bounds = (-10,1000)
    model.reactions.EX_cf_e.bounds = (-0,1000)
    model.reactions.EX_h2_e.bounds = (-0,1000)
    
    #aa uptake
    for reaction_id in aa_list:
        reaction = model.reactions.get_by_id(reaction_id)
        reaction.bounds = (-3,1000)

#simulate
sol_29 = model1.optimize()
sol_30 = model2.optimize()


#test 1mM amino acids
for model in models:
    #acceptors/donors
    model.reactions.EX_dcm_e.bounds = (-10,1000)
    model.reactions.EX_cf_e.bounds = (-0,1000)
    model.reactions.EX_h2_e.bounds = (-0,1000)
    
    #aa uptake
    for reaction_id in aa_list:
        reaction = model.reactions.get_by_id(reaction_id)
        reaction.bounds = (-1,1000)

#simulate
sol_31 = model1.optimize()
sol_32 = model2.optimize()


# #### Mode 3B: + NADH hydrogenase

# In[112]:


for model in models:
    #acceptors/donors
    model.reactions.EX_dcm_e.bounds = (-10,1000)
    model.reactions.EX_cf_e.bounds = (-0,1000)
    model.reactions.EX_h2_e.bounds = (-0,1000)
    
    #aa uptake
    for reaction_id in aa_list:
        reaction = model.reactions.get_by_id(reaction_id)
        reaction.bounds = (0,1000)
        
    #protons
    model.reactions.HYDA_ech.bounds = (-100,0)
    model.reactions.HYDFDN2.bounds = (-100,0)
    model.reactions.HYD_NADH.bounds = (-100,100)
    
#simulate
sol_33 = model1.optimize()
sol_34 = model2.optimize()


# #### Mode 3C: + NADH hydrogenase and amino acids

# In[113]:


#minimal medium
for model in models:
    #acceptors/donors
    model.reactions.EX_dcm_e.bounds = (-10,1000)
    model.reactions.EX_cf_e.bounds = (-0,1000)
    model.reactions.EX_h2_e.bounds = (-0,1000)
    
    #aa uptake
    for reaction_id in aa_list:
        reaction = model.reactions.get_by_id(reaction_id)
        reaction.bounds = (-1,1000)
        
    #protons
    model.reactions.HYDA_ech.bounds = (-100,0)
    model.reactions.HYDFDN2.bounds = (-100,0)
    model.reactions.HYD_NADH.bounds = (-100,100)
    
#simulate
sol_35 = model1.optimize()
sol_36 = model2.optimize()


# ## 3. Write excel tables

# ### 3.1. FBA summary table (Table S10)

# In[117]:


#list solutions
solutions = [sol_1, sol_2, sol_3, sol_4, sol_5, sol_6, sol_7, sol_8, sol_9, sol_10, sol_11, sol_12, sol_13,
             sol_14, sol_15, sol_16, sol_17, sol_18, sol_19, sol_20, sol_21, sol_22, sol_23, sol_24, sol_25,
             sol_26, sol_27, sol_28, sol_29, sol_30, sol_31, sol_32, sol_33, sol_34, sol_35, sol_36,]

# Initialize an empty DataFrame to store the flux values and reactions
combined_table = pd.DataFrame(columns=['Reaction'])

for index, solution in enumerate(solutions):
    # Extract flux values from the solution object
    flux_values = solution.fluxes

    # Create a DataFrame for the current solution
    solution_table = pd.DataFrame({'Reaction': flux_values.index, 
                                   f'Solution_{index+1}_Flux': flux_values.values})
    
    # Merge the current solution DataFrame with the combined table
    combined_table = pd.merge(combined_table, solution_table, on='Reaction', how='outer')

# Save the combined table to a CSV file
combined_table.to_csv("250214_combined_flux_table.csv", index=False)

# Display the combined table if needed
display(combined_table)


# 

# ### 3.2. Select flux summary table (Table S11)

# In[118]:


# Specify the reactions you want to include in the new table
selected_reactions = ['BIOMASS', 'ATPS4rpp', 'EX_cf_e', 'EX_dcm_e', 'EX_h2_e', 'EX_ac_e',
                       "EX_ala__L_e", "EX_arg__L_e", "EX_asn__L_e","EX_asp__L_e","EX_gln__L_e",
                       "EX_glu__L_e","EX_gly_e","EX_his__L_e","EX_ile__L_e","EX_leu__L_e",
                       "EX_lys__L_e", "EX_met__L_e","EX_pro__L_e","EX_thr__L_e","EX_trp__L_e",
                       "EX_tyr__L_e","EX_val__L_e","EX_phe__L_e", "EX_cys__L_e", "EX_ser__L_e"]

# Filter the combined_table to include only the selected reactions
ATPe_table = combined_table.loc[combined_table['Reaction'].isin(selected_reactions)]

# Save the additional table to a CSV file
ATPe_table.to_csv("250214_combined_flux_select.csv", index=False)

# Display the additional table if needed
display(ATPe_table)


# In[ ]:




