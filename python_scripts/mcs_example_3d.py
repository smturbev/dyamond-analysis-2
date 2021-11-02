#!/usr/bin/env python
# %%
import matplotlib
matplotlib.use("Agg")
import xarray as xr
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from utility import load
import numpy as np
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from utility import analysis_parameters as ap

is_one_by_one = False
t = 26
greys = cm.get_cmap("gist_yarg", 28)
model="FV3"
#%%
transp_grays = greys(range(20))
transp_grays = transp_grays[:,:]
print("greys(range(20))", transp_grays)
# transp = np.zeros(transp_grays.shape)
# for i in range(21):
#     transp[-i,:-1] = transp_grays[-i,:-1]
#     if is_one_by_one:
#         transp[-i,-1] = 0.8-(0.85*i/30) # 1-(np.sqrt(i)/18) # 0.6-(0.85*i/30)
#     else:
#         transp[-i,-1] = 0.8-(0.85*i/30) # 0.3-(np.sqrt(i)/18) # 0.6-(0.85*i/30)
transp_grays[:,-1] = (np.linspace(0.25,0.75,20))**2
print("transp cmap:", transp_grays)
new_cmap = ListedColormap(transp_grays)
del transp_grays
#%%
print("load total water content")
qn = xr.open_dataset(ap.FV3+"FV3_twc_TWP.nc")["iwc"][16:]
print("got it.")
#%%
qn = qn[t] * 1000 # convert to g/m3
print("Converted to g/m3")
print(qn.shape, qn.time.values)
#%%
tstring = str(qn.time.values).split("T")[-1][:5]+" UTC "+str(qn.time.values).split("-")[2][:2]+" Aug 2016"
print(tstring)
if model=="FV3":
    z = load.get_levels(model, "TWP")
else:
    z = load.get_levels(model, "TWP")

if is_one_by_one:
    lat = qn.lat.values
    lon = qn.lon.values
else:
    lat = qn.grid_yt.values
    lon = qn.grid_xt.values
x3d = np.repeat(np.repeat(lat[np.newaxis,:,np.newaxis],\
                          qn.shape[0], axis=0), qn.shape[2], axis=2)
y3d = np.repeat(np.repeat(lon[np.newaxis,np.newaxis,:],\
                          qn.shape[0], axis=0), qn.shape[1], axis=1)
z3d = np.repeat(np.repeat(z[:,np.newaxis, np.newaxis],\
                          qn.shape[1], axis=1), qn.shape[2], axis=2)

print("q shape, lat, lon", qn[:,:,:].shape, x3d.shape, z3d.shape)

if is_one_by_one:
    xskip, yskip, zskip = 1, 1, 1
else:
    xskip, yskip, zskip = 5, 4, 1
in_cloud = (np.where(qn[::zskip,::xskip,::yskip] > 5e-4, True, False)).flatten()
qn = qn[::zskip,::xskip,::yskip].values.flatten()[in_cloud]
print("got incloud and flattened array.\nMaking plot...")
# plt.rcParams["image.cmap"] = new_cmap
# %%
fs = 14
fig = plt.figure(figsize=(10,7))
ax = fig.add_subplot(111, projection="3d")
ax.view_init(30,30) #-220
sc = ax.scatter((x3d[::zskip,::xskip,::yskip]).flatten()[in_cloud], \
                (y3d[::zskip,::xskip,::yskip]).flatten()[in_cloud], \
                (z3d[::zskip,::xskip,::yskip]).flatten()[in_cloud]/1000,\
                c=(qn), \
                edgecolors='face', vmax=0.6, vmin=0.1,\
                s=4*xskip, cmap=new_cmap, depthshade=True) #vmax = np.max(qn[:,:,:])/1.5,
# if not(is_one_by_one):
#     ax.plot([  0,  0,  -1,  -1,  0,  0, 0, -1, -1, 0],
#             [147,148,148,147,147,147,148,148,147,147],
#             [ 20 ,20, 20, 20, 20,  0, 0, 0, 0, 0], c='r')
#     ax.plot([    0,  0,  0,  0,  -1,  -1,  -1,  -1,  -1],
#             [147,147,148,148,148,148,148,147,147],
#             [ 0 , 20, 20, 0,  0,  20, 20, 20, 0], c='r')
ax.set_zlim(0,20)
ax.set_ylabel("\nLongitude ($^\circ$E)", fontsize=fs)
ax.set_xlabel("\nLatitude ($^\circ$N)", fontsize=fs)
ax.set_zlabel("Height (km)", fontsize=fs)
fig.suptitle(tstring, fontsize=fs)
if is_one_by_one:
    ax.set_yticks(np.arange(147,148,0.2))
    ax.set_yticklabels([147,None,147.4,None,147.8,None])
    ax.set_xticks(np.arange(-1,0,0.2))
    ax.set_xticklabels([-1,None,-0.6,None,-0.2,None])

cbar = plt.colorbar(sc, ax=ax, shrink=0.6, extend="max")
cbar.set_label("water content (g m$^{-3}$)", fontsize=fs-2)
ax.w_xaxis.set_pane_color((.5,.58,1.0))
ax.w_yaxis.set_pane_color((.5,.58,1.0))
ax.w_zaxis.set_pane_color((.2,.2,1.0))
ax.tick_params(labelsize=fs-2)
cbar.ax.tick_params(labelsize=fs-2)
if is_one_by_one:
    plt.savefig(f"../plots/fv3/cloud3d_{t}.png", transparent=False, dpi=200, bbox_inches="tight", pad_inches=0.2)
else:
    ax.set_xlim([-5,5])
    ax.set_ylim([143,153])
    ax.invert_xaxis()
    # ax.invert_yaxis()
    plt.savefig(f"../plots/fv3/cloud3d_10x10_long_{t}.png", transparent=False, dpi=200, bbox_inches="tight", pad_inches=0.2)
plt.close()
print("plot saved in ../plots/fv3/")


# %%

# %%
