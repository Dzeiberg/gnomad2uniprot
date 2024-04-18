#!/usr/bin/env python
# coding: utf-8

# # Mapping gnomAD v2.1 to uniprot

# In[4]:


import pandas as pd
from collections import defaultdict
from Bio.PDB.Polypeptide import protein_letters_3to1
from tqdm.autonotebook import tqdm


# ## Read Data

# ### Read MANE Summary table to obtain Ensembl protein IDs of MANE-Select transcripts

# In[5]:


mane_summary = pd.read_csv("https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/MANE.GRCh38.v1.3.summary.txt.gz",
                           compression="gzip",sep="\t")


# ### Read in gnomAD v2.1 data (transformed to table using gatk)  

# In[6]:


with open("header.txt") as f:
    header = f.read().split(',')
exomes = pd.read_csv("gnomad.exomes.r2.1.1.sites.missense.vcf.gz",compression="gzip",sep="\t",header=None)
exomes.columns = header


# ### Filter only variants in MANE-Select transcripts

# In[7]:


gnomad_mane = exomes[exomes.HGVSp.str.split(":").str[0].isin(set(mane_summary[mane_summary.MANE_status == "MANE Select"].Ensembl_prot.values))]


# ### Read in uniprot id mapping table for human genes

# In[8]:


id_mapping = pd.read_csv("https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz",
                         compression="gzip",sep='\t',header=None)


# #### Make dictionary mapping Ensembl_prot -> uniprot_acc

# In[9]:


ensembl_prot_2_uniprot = {}
for _,r in id_mapping[id_mapping.loc[:,1] == "Ensembl_PRO"].iterrows():
        ensembl_prot_2_uniprot[r[2]] = r[0]


# #### Identify any Ensembl_prot ids that are missing from this mapping table

# In[10]:


set(gnomad_mane.HGVSp.str.split(":").str[0].values) - set(ensembl_prot_2_uniprot.keys())


# #### Hard code these missing mappings (search Ensembl_prot id on uniprot website)

# In[11]:


ensembl_prot_2_uniprot["ENSP00000221455.3"] = "Q8N2M8"
ensembl_prot_2_uniprot['ENSP00000265310.1'] = "Q9NQA5"
ensembl_prot_2_uniprot['ENSP00000289575.5'] = 'O94956'
ensembl_prot_2_uniprot['ENSP00000365465.2'] = 'Q15517'
ensembl_prot_2_uniprot['ENSP00000369823.3'] = 'Q3MIN7'
ensembl_prot_2_uniprot['ENSP00000405012.2'] = 'Q6ZQQ6'


# In[12]:


df = gnomad_mane.loc[gnomad_mane.Consequence == "missense_variant",["CHROM", "POS", "ID","REF", "ALT", "QUAL","FILTER","AC", "AN","AF","HGVSp"]]
df = df.assign(Ensembl_prot = df.HGVSp.str.split(":").str[0])
df = df.set_index("Ensembl_prot").assign(uniprot_acc = ensembl_prot_2_uniprot).reset_index()
df = df.assign(aa_ref = df.HGVSp.str.split(":").str[1].str[2:5].str.upper().apply(lambda x: protein_letters_3to1[x] if x in protein_letters_3to1 else x),
                aa_alt = df.HGVSp.str.split(":").str[1].str[-3:].str.upper().apply(lambda x: protein_letters_3to1[x] if x in protein_letters_3to1 else x),
                aa_pos = df.HGVSp.str.split(":").str[1].str[5:-3].astype(int))


# In[13]:


from scipy.io import savemat


# In[14]:


savemat("gnomad_mane_mapped.mat",{k:df[k].values for k in df.columns},do_compression=True)


# In[ ]:




