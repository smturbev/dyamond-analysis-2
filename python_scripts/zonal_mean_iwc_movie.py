# %%
"""
Plots the total water content (3hrly zonal mean). Outputs png images.
Use ffmpeg to convert to movie.
"""
import matplotlib.pyplot as plt
import numpy as np

from utility import load

model, region = "FV3", "TWP"

if (model=="NICAM") or (model=="SAM"):
    lat = "lat"
elif model=="FV3":
    lat = "grid_yt"
else:
    raise Exception("model not defined")

# %%
# get total water content (saved only for TWP)
q = load.get_twc(model, region)
iwp = load.get_iwp(model, region, ice_only=False)
q = q.where(iwp>=1e-4)
# %%
del iwp
# %%
# calculate zonal mean in g/m3
q_zonal_mean = q.mean(axis=3)*1000 # g/m3
# %%
# get altitude for plotting
z = load.get_levels(model, region)
# %%
# plot for each time step
# for t in range(q.shape[0]-1):
t = 384//16
# cf_mean = np.mean(np.where(np.isnan(q.values),0,1), axis=1)
# print(cf_mean)
plt.figure(figsize=(7,3))
plt.pcolormesh(q_zonal_mean[lat].values, z/1000, 
            (q_zonal_mean[t]), cmap="ocean_r",
            vmin=5e-4, vmax=0.2)
# plt.pcolormesh(q.y_grid.values,np.array([0]*z.shape[0]),cf_mean)
plt.xlabel("Latitude")
plt.ylabel("Height (km)")
plt.ylim([0,20])
plt.title(region+" "+model+" ice + liquid water content (cldy)")
plt.colorbar(label="meridional mean\nwater content (g/m3)")
plt.savefig("../plots/animation/{}_{:0>3}_zonal_mean.png".format(model,t), dpi=120)
plt.close()

# %%
