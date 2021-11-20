#!/usr/bin/env python
# coding: utf-8

# In[3]:


from pyopenms import *
seq = AASequence.fromString("VAKA")
V_weight = seq.getMonoWeight()
A_weight = seq.getMonoWeight()
K_weight = seq.getMonoWeight()
A_weight = seq.getMonoWeight()
print("Monoisotopic mass of peptide [VAKA] is",V_weight)

print("The peptide", str(seq), "consists of the following amino acids:")
for aa in seq:
    print(aa.getName(), ":", aa.getMonoWeight())


# In[ ]:




