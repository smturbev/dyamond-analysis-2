# %%
from utility import load
import matplotlib.pyplot as plt
import numpy as np

model = "SAM"
region = "TWP"
# %%
temp = load.get_temp(model, region)
z    = load.get_levels(model, region)
pres = load.get_pres(model, region)
# %%
temp
# %%
ind = np.argmin(abs(temp.pfull-225))
# print(temp.pfull[ind].values, z[ind])
print()
# %%
temp_lev = temp[:,ind]
temp_mean = np.nanmean(temp_lev)
temp_perturb = temp_lev - temp_mean

# %%
