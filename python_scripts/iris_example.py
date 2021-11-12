# %%
# practice with iris and tobac
import iris.plot as iplt
import iris.quickplot as qplt
import iris
import matplotlib.pyplot as plt
from scipy.ndimage.measurements import label
import xarray as xr
import numpy as np
import tobac
import pandas as pd

import utility.analysis_parameters as ap
from utility import load, load01deg, util
RUN = True
thres = 1e-4

# %%
# olr = load01deg.get_olr("FV3", "TWP")[96*2:960]
# olr_cube = olr.to_iris()

# %%
iwp_cube = load01deg.get_iwp("FV3","TWP")[96*2:400].to_iris()

# %%
t = 100
# for t in range(448,600,4):
plt.figure(figsize=(6,6))
qplt.pcolormesh(iwp_cube[t,:,:], cmap="viridis_r",
            vmin=1e-7, vmax=1)
plt.gca().coastlines()
plt.xticks(np.arange(143,154,2))
plt.yticks(np.arange(-5,6,2))
plt.title(util.tstring(t*0.25))
plt.savefig(f"../plots/fv3/iwp_cube_plot_{t}.jpg")
plt.show()

# %%
## olr coord name change
# olr_cube.coord("latitude").standard_name  = "projection_x_coordinate"
# olr_cube.coord("longitude").standard_name = "projection_y_coordinate"

# %%
# Feature Detection
## Determine temporal and spatial sampling:
dxy,dt = tobac.get_spacings(iwp_cube, grid_spacing=10)

# %%
# Dictionary containing keyword arguments for feature detection step (Keywords could also be given directly in the function call).
if RUN:
    parameters_features={}
    parameters_features['position_threshold']='weighted_diff'
    parameters_features['sigma_threshold']=0.5
    parameters_features['min_num']=4
    parameters_features['target']='maximum'
    if thres==1e-2:
        parameters_features['threshold']=[1,1e-2]
    elif thres==1e-4:
        parameters_features['threshold']=[1,1e-2,thres]
    print("thresholds", parameters_features['threshold'])
# %%
## Perform feature detection:
if RUN:
    print('starting feature detection')
    Features=tobac.feature_detection_multithreshold(iwp_cube,dxy,**parameters_features)
    print(Features)
    Features.to_hdf('../plots/fv3/Features.h5','table')
    print('feature detection performed and saved')

# %%
# Segmentation step
## Dictionary containing keyword options for the segmentation step:
if RUN:
    parameters_segmentation={}
    parameters_segmentation['target']='maximum'
    parameters_segmentation['method']='watershed'
    parameters_segmentation['threshold']=thres

# %%
## Perform segmentation and save results:
if RUN:
    print('Starting segmentation based on IWP.')
    Mask_OLR,Features_OLR=tobac.segmentation_2D(Features,iwp_cube,dxy,**parameters_segmentation)
    print('segmentation OLR performed, start saving results to files')
    iris.save([Mask_OLR],'../plots/fv3/Mask_Segmentation_IWP.nc',zlib=True,complevel=4)                
    Features_OLR.to_hdf('../plots/fv3/Features_IWP.h5','table')
    print('segmentation IWP performed and saved')
# %%
# Trajectory linking
## Arguments for trajectory linking:
if RUN:
    parameters_linking={}
    parameters_linking['v_max']=20
    parameters_linking['stubs']=2
    parameters_linking['order']=1
    parameters_linking['extrapolate']=1
    parameters_linking['memory']=0
    parameters_linking['adaptive_stop']=0.2
    parameters_linking['adaptive_step']=0.95
    parameters_linking['subnetwork_size']=100
    parameters_linking['method_linking']= 'predict'

    ## Perform linking and save results to file:
    Track=tobac.linking_trackpy(Features,iwp_cube,dt=dt,dxy=dxy,**parameters_linking)
    Track.to_hdf('../plots/fv3/Track.h5','table')

# %%
# Visualization
## Set extent of maps created in the following cells:
axis_extent=[143,153,-5,5]
# %%
# load saved track analysis
if not(RUN):
    Track    = pd.read_hdf("../plots/fv3/Track.h5", "table")
    Features = pd.read_hdf("../plots/fv3/Features.h5", "table")
    Mask_OLR = iris.load_cube("../plots/fv3/Mask_Segmentation_IWP.nc")
    Features_OLR = pd.read_hdf("../plots/fv3/Features_IWP.h5", "table")
    # Track.rename(columns = {'projection_x_coordinate':'longitude',
    #             'projection_y_coordinate':'latitude'}, inplace = True)


# %%
## Plot map with all individual tracks:
import cartopy.crs as ccrs
fig_map,ax_map=plt.subplots(figsize=(8,5),
                    subplot_kw={'projection': ccrs.PlateCarree()})
ax_map=tobac.map_tracks(Track,axis_extent=axis_extent,axes=ax_map)
plt.savefig(f"../plots/fv3/track_map_{thres}.png", transparent=False)
print("saved")
plt.show()
# %%
# Create animation of tracked clouds and outlines with OLR as a background field
animation_test_tobac=tobac.animation_mask_field(Track,Features,iwp_cube,Mask_OLR,
                                          axis_extent=axis_extent,#figsize=figsize,orientation_colorbar='horizontal',pad_colorbar=0.2,
                                          vmin=1e-4,vmax=1,extend="max",
                                          plot_outline=True,plot_marker=True,marker_track='x',plot_number=True,plot_features=True)


# %%
## Display animation:
# from IPython.display import HTML, Image, display
# HTML(animation_test_tobac.to_html5_video())

# %%
# Save animation to file:
savefile_animation=f'../plots/fv3/Animation_{thres}.mp4'
animation_test_tobac.save(savefile_animation,dpi=200)
print(f'animation saved to {savefile_animation}')

# %%
# Lifetimes of tracked clouds:
fig_lifetime,ax_lifetime=plt.subplots()
tobac.plot_lifetime_histogram_bar(Track,axes=ax_lifetime,bin_edges=np.arange(0,1440,60),density=False,width_bar=60)
ax_lifetime.set_xlabel('lifetime (hr)')
ax_lifetime.set_ylabel('counts')
ax_lifetime.set_yscale('log')
ax_lifetime.set_xticks(np.arange(0,1441,60))
ax_lifetime.set_xticklabels(range(0,25), rotation=30)
plt.savefig(f"../plots/fv3/lifetime_histogram_{thres}.png",facecolor='white', transparent=False)
plt.show()
# %%
fig_area,ax_area=plt.subplots(1,1)
tobac.area_histogram(Features,Mask_OLR,bin_edges=np.arange(0,30000,500),
                   density=False,method_area=None,
                   return_values=False,representative_area=False)
ax_area.set_xlabel('Area (m)')
ax_area.set_ylabel('counts')
plt.savefig(f"../plots/fv3/area_histogram_{thres}.png")
print("done")
plt.show()
# %%
track18 = Track[Track.cell==18]
print(track18)
#%%
plt.subplots(1,1, figsize=(5,4))
sc = plt.scatter(track18.longitude, track18.latitude, c=track18.frame, 
            s=5, cmap="cividis")
plt.title("cell=18, track path, life time = "+str(list(track18.time_cell)[-1]))
plt.colorbar(sc, label="timesteps since "+str(list(track18.time)[0])[:-3]+" UTC")
plt.xlim([143,153])
plt.ylim([-5,5])
plt.savefig("../plots/fv3/track18.jpg")
plt.show()
# %%
track4475 = Track[Track.cell==4475]
print(track4475)
# %%
plt.subplots(1,1, figsize=(5,4))
sc = plt.scatter(track4475.longitude, track4475.latitude, c=track4475.frame, 
            s=5, cmap="cividis")
plt.title("cell=4475, track path, life time = "+str(list(track4475.time_cell)[-1]))
plt.colorbar(sc, label="timesteps since "+str(list(track4475.time)[0])[:-3]+" UTC")
plt.xlim([143,153])
plt.ylim([-5,5])
plt.savefig("../plots/fv3/track4475.jpg")
plt.show()
# %%
