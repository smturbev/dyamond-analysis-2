# %%
import matplotlib.pyplot as plt
import numpy as np

from utility import load, load01deg, util

model = "FV3"
region = "TWP"

DC = False

# %%
# define methods
def plot_wp(t, wp, iwp, lwp, xlim=[None,None], ylim=[None,None], figsize=(10,4), save=False):
    """Returns fig with 3 panels of total water path, frozen and liquid water paths"""
    fig, ax = plt.subplots(1,2, figsize=figsize, sharex=True, sharey=True)

    # pc = ax[0].pcolormesh(twp.grid_xt, twp.grid_yt, np.log10(wp[t]*1000), cmap="gist_earth_r",
    #                         vmax=5, vmin=-1)
    # cbar = plt.colorbar(pc, ax=ax[0])
    # cbar.set_label("log (g/m2)")
    # # ax[0].set_title("total water path >= %s kg/m2"%thres)
    # ax[0].set_title("water path (g/m2)")

    pc = ax[0].pcolormesh(twp.grid_xt, twp.grid_yt, np.log10(iwp[t]*1000), cmap="gist_earth_r",
                            vmax=5, vmin=-1)
    cbar = plt.colorbar(pc, ax=ax[0])
    cbar.set_label("log (g/m2)")
    ax[0].set_title("frozen water path (g/m2)")

    pc = ax[1].pcolormesh(twp.grid_xt, twp.grid_yt, np.log10(lwp[t]*1000), cmap="gist_earth_r",
                            vmax=5, vmin=-1)
    cbar = plt.colorbar(pc, ax=ax[1])
    cbar.set_label("log (g/m2)")
    ax[1].set_title("liquid water path (g/m2)")

    ax[0].grid()
    ax[1].grid()

    ax[0].set_xlim(xlim)
    ax[0].set_ylim(ylim)
    fig.suptitle("time step = "+str(t)+" "+(util.tstring((t*0.25)+48)))
    if save:
        if DC:
            plt.savefig("../plots/fv3/case_study_"+region+"_f{0:03d}.png".format(t), bbox_inches="tight")
        else:
            plt.savefig("../plots/fv3/case_study_"+region+"_ci_f{0:03d}.png".format(t), bbox_inches="tight")
    plt.show()

# %%
# load frozen and liquid water paths
# iwp = load01deg.get_iwp(model, region, ice_only=False)
# lwp = load01deg.get_lwp(model, region)[:,0]
iwp = load.get_iwp(model, region, ice_only=False)
lwp = load.get_lwp(model, region, rain=False)
twp = iwp + lwp
# ttliwp = load.get_ttliwp(model, region)

# %%
# deep convective event with total water path > thres g/m2
# list all times with such events then plot them
if DC:
    thres = 60
    twp_dc = twp.where(twp>=thres) # kg/m2 
    times_dc = []
    for t in range(len(twp.time)):
        if (twp[t]>=thres).any():
            times_dc.append(t)
    print(len(times_dc))
else:
    thres = 200
    times_ci = []
    twp_ci = np.zeros(twp.shape)
    for t in range(len(twp.time)):
        twp_ci[t] = twp[t].where((twp[t]<thres)&(iwp[t].values>=1e-4)&(lwp[t]<1e-4)) # kg/m2 
        if (~np.isnan(twp_ci[t])).any(): 
            times_ci.append(t)
    print(len(times_ci))
# %%
# plot twp,  iwp, lwp 
if DC:
    t = times_dc[0]
    print(t)
    plot_wp(t, twp_dc, iwp, lwp, save=True)
    print("done")
else:
    t = times_ci[0]
    print(t)
    plot_wp(t, twp_ci, iwp, lwp, save=False)
    print("done")

# %%
# time evolution & zoomed
print(115)
for i in range(3500):
    print(i)
    if DC:
        plot_wp(i, twp_dc, iwp, lwp, save=True)
    else:
        plot_wp(i, twp_ci, iwp, lwp, save=True)
# %%
