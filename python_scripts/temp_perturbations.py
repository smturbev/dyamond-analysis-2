# %%
from numpy.core.numeric import False_
from utility import load, util
import matplotlib.pyplot as plt
import numpy as np

model = "FV3"
region = "TWP"
# %%
# load variables
temp = load.get_temp(model, region)
z    = load.get_levels(model, region)
pres = load.get_pres(model, region)
fwc = load.load_frozen(model, region, ice_only=True)
fwp = load.get_iwp(model, region, ice_only=False)[11::12]
# %%
# pres near 210 mb
ind = np.argmin(abs(temp.pfull-210))
print(temp.pfull[ind].values, z[ind])
# %%
# scatter plot of temperature perturbations 
# at near 210 mb
mean_temp = np.nanmean(temp[:,ind,:])
for i in range(100):
    fig, ax = plt.subplots(1,3, figsize=(18,5), sharex=True, sharey=True)
    #temp perturbation
    scat = ax[0].scatter(temp.lon, temp.lat, c=(temp[i,ind,:]-mean_temp),
                         vmin=-2.5,vmax=2.5,cmap="bwr")
    ax[0].set_title("Temp perturbation (K) at "+str(temp.pfull.values[ind])+"mb\nmean="+str(mean_temp)+" K")
    plt.colorbar(scat, ax=ax[0])
    # fwc at same pressure level as temp
    pc = ax[1].pcolormesh(fwc.grid_xt, fwc.grid_yt, np.log10(fwc[i,ind,:,:]), 
                            vmin=np.log10(5e-5), vmax=-1, cmap="gist_earth_r")
    plt.colorbar(pc, ax=ax[1])
    ax[1].set_title("FWC log10(kg/m3) at "+str(temp.pfull.values[ind])+"mb")

    pc = ax[2].pcolormesh(fwc.grid_xt, fwc.grid_yt, np.log10(fwp[i,:,:]), 
                            vmin=-4, vmax=1, cmap="gist_earth_r")
    plt.colorbar(pc, ax=ax[2])
    ax[2].set_title("FWP log10(kg/m2)")

    ax[0].set_xlim([143,153])
    ax[0].set_ylim([-5,5])

    ax[0].grid(True)
    ax[1].grid(True)
    ax[2].grid(True)

    fig.suptitle("t="+str(i)+", "+(util.tstring(i*3+48)))

    plt.savefig("../plots/case_study_fv3/temp/temp_{}_{}_{}.png".format(model.lower(), 
                region.lower(), str(i)), bbox_inches="tight")
