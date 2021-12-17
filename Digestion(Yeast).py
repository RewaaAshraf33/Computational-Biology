#!/usr/bin/env python
# coding: utf-8

# In[2]:


from pyopenms import *
from urllib.request import urlretrieve
#gh = "https://raw.githubusercontent.com/OpenMS/pyopenms-extra/master"
#urlretrieve (gh + "/src/data/P02769.fasta", "bsa.fasta")

dig = ProteaseDigestion()
dig.getEnzymeName() # Trypsin
bsa = "".join([l.strip() for l in open("uniprot-yourlist_M2021121792C7BAECDB1C5C413EE0E0348724B682374277J.fasta").readlines()[1:]])
bsa = AASequence.fromString(bsa)
# create all digestion products
result = []
dig.digest(bsa, result)
#print(result[4].toString())
#len(result)
for s in result:
    print(s.toString())


# In[ ]:




