
# %%
import numpy as np
import matplotlib.pyplot as plt

from utility import load, util

models = ["SAM","NICAM","ICON"] # need to figure out fix for FV3's mismatched grid points
regions = ["SHL","NAU"] # done with TWP

# %%
for r in regions:
    for m in models:
        print("STARTING "+m+" "+r)
        # load olr, alb, w, z
        olr, alb = load.get_olr_alb(m,r)
        w = load.get_w(m,r)
        z = load.get_levels(m,r)

        ind14km = np.argmin(abs(z-14000))
        print(ind14km, z[ind14km])

        # strong updraft w > 0.3 m/s or or 99th percentile
        # weak ascent 0.001 < w < 0.3 m/s or 50th percentile
        # subsidence w < 0.001
        w_thres = [0.3,0.001]
        w_string = ["Strong Ascent", "Weak Ascent", "Subsidence"]
        if len(olr.shape)==4:
            olr = olr[:,0]
        if len(alb.shape)==4:
            alb = alb[:,0]

        for i in range(3):
            print(olr.shape, alb.shape, w[:,ind14km].shape)
            if i==0:
                print(i, ">", w_thres[0])
                olr_binned = np.where(w[:,ind14km]>w_thres[i], olr, np.nan)
                alb_binned = np.where(w[:,ind14km]>w_thres[i], alb, np.nan)
            elif i==1:
                print(i, ">", w_thres[1]," and <=", w_thres[0])
                olr_binned = np.where((w[:,ind14km]>w_thres[1])&(w[:,ind14km]<=w_thres[0]), 
                                olr, np.nan)
                alb_binned = np.where((w[:,ind14km]>w_thres[1])&(w[:,ind14km]<=w_thres[0]),
                                alb, np.nan)
            else:
                print(i, "<", w_thres[1])
                olr_binned = np.where(w[:,ind14km]<=w_thres[1], olr, np.nan)
                alb_binned = np.where(w[:,ind14km]<=w_thres[1], alb, np.nan)

            print(w_string[i], np.nanmean(olr), np.nanmean(alb))
            if (m=="NICAM"):
                time_array = w.time.dt.hour.values
            elif m=="ICON":
                time_array = w.t.dt.hour.values
            elif (m=="SAM"):
                time_array = np.array(list(np.arange(0,24,3))*38)
            else:
                raise Exception("model not defined",m)
            if r=="TWP":
                time_mask = (time_array>=20)
            elif r=="NAU":
                time_mask = np.where((time_array>=22)|(time_array<=2), True, False)
            elif r=="SHL":
                print(time_array.shape, time_array[:5])
                time_mask = np.where((time_array>=11)&(time_array<=15), True, False)
            else:
                raise Exception("region not defined,", r)
            olr_binned = olr_binned[time_mask]
            alb_binned = alb_binned[time_mask]

            fig, ax = plt.subplots(1,1,figsize=(6.5,6))
            util.dennisplot("density", olr_binned.flatten(), alb_binned.flatten(), 
                        ax=ax, model="Native "+m, region=r+"\n"+w_string[i])
            plt.savefig("../plots/w_alb-olr_hist_{}_{}_{}.png".format(m,r,
                        w_string[i].lower().replace(" ","_")),
                        bbox_inches="tight",pad_inches=0.5)
            plt.show()

# %%
