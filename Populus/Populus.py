import numpy as np
import vcfnp
import os

filename = "comt.chr06.snp.full.final.vcf"
c = vcfnp.calldata_2d(filename, cache=True).view(np.recarray)

G = np.sum(c.genotype.T, axis=0)
G[G == -2] = 3


# In[7]:

(n,m) = G.shape
nb_missing_data = np.count_nonzero(G == 3)
nb_data = n*m
nb_0 = np.count_nonzero(G == 0)
nb_1 = np.count_nonzero(G == 1)
nb_2 = np.count_nonzero(G == 2)

print("nb data : {0}".format(nb_data))
print("nb missing data : {0}".format(nb_missing_data))
print("nb 0 : {0}".format(nb_0))
print("nb 1 : {0}".format(nb_1))
print("nb 2 : {0}".format(nb_2))


# In[8]:

float(nb_missing_data)/nb_data


# In[9]:

nb_data == (nb_missing_data + nb_0 + nb_1 + nb_2)


# In[10]:

mask = np.logical_and.reduce([G <= 3], axis=1)


# In[11]:

G_filtered = G[:,mask[0,:]]


# In[12]:

(n,m) = G_filtered.shape
nb_missing_data = np.count_nonzero(G_filtered == 3)
nb_data = n*m
nb_0 = np.count_nonzero(G_filtered == 0)
nb_1 = np.count_nonzero(G_filtered == 1)
nb_2 = np.count_nonzero(G_filtered == 2)

print("nb data : {0}".format(nb_data))
print("nb missing data : {0}".format(nb_missing_data))
print("nb 0 : {0}".format(nb_0))
print("nb 1 : {0}".format(nb_1))
print("nb 2 : {0}".format(nb_2))


# In[13]:

float(nb_missing_data)/nb_data


# In[14]:

nb_data == (nb_missing_data + nb_0 + nb_1 + nb_2)
