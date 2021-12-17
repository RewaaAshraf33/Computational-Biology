#!/usr/bin/env python
# coding: utf-8

# In[1]:


from pyopenms import *
from urllib.request import urlretrieve
#Human 
# Y-ion spectrum hum.
dig = ProteaseDigestion()
dig.getEnzymeName()
bsa = "".join([l.strip() for l in open("uniprot-yourlist_M2021121792C7BAECDB1C5C413EE0E0348724B682373C6DM.fasta").readlines()[1:]])
bsa = AASequence.fromString(bsa)
# create all digestion products
result = []
dig.digest(bsa, result)
peptides = [AASequence.fromString(s.toString()) for s in result]

# Iterate over annotated ions and their masses
for peptide in peptides:
    tsg = TheoreticalSpectrumGenerator()
    spec1 = MSSpectrum()

    # standard behavior is adding b- and y-ions of charge 1
    p = Param()
    p.setValue("add_b_ions", "false")
    p.setValue("add_metainfo", "true")
    tsg.setParameters(p)
    tsg.getSpectrum(spec1, peptide, 1, 1) # charge range 1:1
    print("Spectrum 1 of", peptide, "has", spec1.size(), "peaks.")
    for ion, peak in zip(spec1.getStringDataArrays()[0], spec1):
        print(ion.decode(), "is generated at m/z", peak.getMZ())


# In[2]:


#yeast
# Y-ion spectrum yeast
dig = ProteaseDigestion()
dig.getEnzymeName()
bsa = "".join([l.strip() for l in open("uniprot-yourlist_M2021121792C7BAECDB1C5C413EE0E0348724B682374277J.fasta").readlines()[1:]])
bsa = AASequence.fromString(bsa)
# create all digestion products
result = []
dig.digest(bsa, result)
peptides = [AASequence.fromString(s.toString()) for s in result]

# Iterate over annotated ions and their masses
for peptide in peptides:
    tsg = TheoreticalSpectrumGenerator()
    spec1 = MSSpectrum()

    # standard behavior is adding b- and y-ions of charge 1
    p = Param()
    p.setValue("add_b_ions", "false")
    p.setValue("add_metainfo", "true")
    tsg.setParameters(p)
    tsg.getSpectrum(spec1, peptide, 1, 1) # charge range 1:1
    print("Spectrum 1 of", peptide, "has", spec1.size(), "peaks.")
    for ion, peak in zip(spec1.getStringDataArrays()[0], spec1):
        print(ion.decode(), "is generated at m/z", peak.getMZ())


# In[ ]:




