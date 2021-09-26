# %%
# practice with iris and tobac
import iris.plot as iplt
import iris.quickplot as qplt
import iris
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import tobac
import pandas as pd

import utility.analysis_parameters as ap
from utility import load, load01deg, util
RUN = False

# %%
olr = load01deg.get_olr("FV3", "TWP")[96*2:960]
olr_cube = olr.to_iris()

# %%
# t = 268
for t in range(448,600,4):
    plt.figure(figsize=(6,6))
    qplt.pcolormesh(olr_cube[t,:,:], cmap="viridis_r",
                vmin=80, vmax=300)
    plt.gca().coastlines()
    plt.xticks(np.arange(143,154,2))
    plt.yticks(np.arange(-5,6,2))
    plt.title(util.tstring(t*0.25))
    plt.savefig(f"../plots/fv3/olr_cube_plot_{t}.png")
    plt.show()

# %%
## olr coord name change
# olr_cube.coord("latitude").standard_name  = "projection_x_coordinate"
# olr_cube.coord("longitude").standard_name = "projection_y_coordinate"

# %%
# Feature Detection
## Determine temporal and spatial sampling:
dxy,dt = tobac.get_spacings(olr_cube, grid_spacing=10)

# %%
# Dictionary containing keyword arguments for feature detection step (Keywords could also be given directly in the function call).
if RUN:
    parameters_features={}
    parameters_features['position_threshold']='weighted_diff'
    parameters_features['sigma_threshold']=0.5
    parameters_features['min_num']=4
    parameters_features['target']='minimum'
    parameters_features['threshold']=np.linspace(125,250,4)
    print("thresholds", parameters_features['threshold'])
# %%
## Perform feature detection:
if RUN:
    print('starting feature detection')
    Features=tobac.feature_detection_multithreshold(olr_cube,dxy,**parameters_features)
    Features.to_hdf('../plots/fv3/Features.h5','table')
    print('feature detection performed and saved')

# %%
# Segmentation step
## Dictionary containing keyword options for the segmentation step:
if RUN:
    parameters_segmentation={}
    parameters_segmentation['target']='minimum'
    parameters_segmentation['method']='watershed'
    parameters_segmentation['threshold']=250

# %%
## Perform segmentation and save results:
if RUN:
    print('Starting segmentation based on OLR.')
    Mask_OLR,Features_OLR=tobac.segmentation_2D(Features,olr_cube,dxy,**parameters_segmentation)
    print('segmentation OLR performed, start saving results to files')
    iris.save([Mask_OLR],'../plots/fv3/Mask_Segmentation_OLR.nc',zlib=True,complevel=4)                
    Features_OLR.to_hdf('../plots/fv3/Features_OLR.h5','table')
    print('segmentation OLR performed and saved')
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
    Track=tobac.linking_trackpy(Features,olr_cube,dt=dt,dxy=dxy,**parameters_linking)
    Track.to_hdf('../plots/fv3/Track.h5','table')

# %%
# Visualization
## Set extent of maps created in the following cells:
axis_extent=[143,153,-5,5]
# %%
# load saved track analysis
if not(RUN):
    Track = pd.read_hdf("../plots/fv3/Track.h5", "table")
    Features = pd.read_hdf("../plots/fv3/Features.h5", "table")
    Mask_OLR = iris.load_cube("../plots/fv3/Mask_Segmentation_OLR.nc")
    Features_OLR = pd.read_hdf("../plots/fv3/Features_OLR.h5", "table")
    # Track.rename(columns = {'projection_x_coordinate':'longitude',
    #             'projection_y_coordinate':'latitude'}, inplace = True)

# %%
## Plot map with all individual tracks:
import cartopy.crs as ccrs
fig_map,ax_map=plt.subplots(figsize=(8,5),
                    subplot_kw={'projection': ccrs.PlateCarree()})
ax_map=tobac.map_tracks(Track,axis_extent=axis_extent,axes=ax_map)
plt.savefig("../plots/fv3/track_map.png", transparent=False)
print("saved")
plt.show()
# %%
# Create animation of tracked clouds and outlines with OLR as a background field
animation_test_tobac=tobac.animation_mask_field(Track,Features,olr_cube,Mask_OLR,
                                          axis_extent=axis_extent,#figsize=figsize,orientation_colorbar='horizontal',pad_colorbar=0.2,
                                          vmin=80,vmax=330,
                                          plot_outline=True,plot_marker=True,marker_track='x',plot_number=True,plot_features=True)


# %%
## Display animation:
# from IPython.display import HTML, Image, display
# HTML(animation_test_tobac.to_html5_video())

# %%
# Save animation to file:
savefile_animation='../plots/fv3/Animation.mp4'
animation_test_tobac.save(savefile_animation,dpi=200)
print(f'animation saved to {savefile_animation}')

# %%
# Lifetimes of tracked clouds:
fig_lifetime,ax_lifetime=plt.subplots()
tobac.plot_lifetime_histogram_bar(Track,axes=ax_lifetime,bin_edges=np.arange(0,200,20),density=False,width_bar=10)
ax_lifetime.set_xlabel('lifetime (min)')
ax_lifetime.set_ylabel('counts')
plt.savefig("../plots/fv3/lifetime_histogram.png")
plt.show()
# %%
