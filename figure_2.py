import numpy as np
from netCDF4 import Dataset
import os
from scipy import stats
from mpl_toolkits.basemap import Basemap 
import matplotlib.pyplot as plt


years = np.arange(1979,2018,1)

# =================================== ERA =============================================
era_dir = ""
os.chdir(era_dir)

# Get all filenames
directory = os.getcwd()

pressures_era = [200,250,300,400,500,600,700,850,925,1000]

# Array to store 6 hourly temps
hourly_polar = np.zeros((len(years),366*4,len(pressures_era)))
hourly_polar[:] = np.nan 

hourly_tropic = np.zeros((len(years),366*4,len(pressures_era)))
hourly_tropic[:] = np.nan 

# Loop over all files
for year in years:
	doy=0
	for filename in sorted(os.listdir(directory)):
		if str(year) in filename[24:28]:
			print filename
			# Open file
			f = Dataset(filename,mode='r')
		
			# Grid
			lats_era = f.variables['g4_lat_2'][:]
			lons_era = f.variables['g4_lon_3'][:]
			levs = f.variables['lv_ISBL1'][:]			

			lats_polar = lats_era[0:29]
			lats_tropic = lats_era[28::]

			levidx1 = np.where(levs==200)[0][0]
			levidx2 = np.where(levs==1000)[0][0]

			weights_polar = np.cos(np.radians(lats_polar))
			weights_tropic = np.cos(np.radians(lats_tropic))

			polarbox = f.variables['T_GDS4_ISBL'][:,levidx1:levidx2+1,0:29,:].squeeze()
			tropicbox = f.variables['T_GDS4_ISBL'][:,levidx1:levidx2+1,28::,:].squeeze()	

			polarbox_zm = np.mean(polarbox,axis=-1)
			tropicbox_zm = np.mean(tropicbox,axis=-1)

			polarbox_mean = np.average(polarbox_zm,axis=-1,weights=weights_polar)	
			tropicbox_mean = np.average(tropicbox_zm,axis=-1,weights=weights_tropic)	

			# Store value
			hourly_tropic[np.where(years==year)[0][0],doy,:] = tropicbox_mean
			hourly_polar[np.where(years==year)[0][0],doy,:] = polarbox_mean

			f.close()
			doy+=1
		
# Calculate annual-mean
ann_mean_polar = np.nanmean(hourly_polar,axis=1)
ann_mean_tropic = np.nanmean(hourly_tropic,axis=1)
ann_mean_gradient = ann_mean_tropic - ann_mean_polar

# Statistical analysis
# Array to store slopes
slopes_era = np.zeros((len(pressures_era)))
# Array to store slope error
confidence_era = np.zeros(np.shape(slopes_era))
# Array to store p values
pvals_era = np.zeros(np.shape(slopes_era))

for p in range(np.shape(ann_mean_gradient)[1]):
	slope,intercept,rvalue,pvalue,stderr = stats.linregress(years,ann_mean_gradient[:,p])
	slopes_era[p] = 10*slope # trend per decade, -1 to flip lats
	confidence_era[p] = 1.96*10*stderr # 95% confidence on slope per decade
	pvals_era[p] = pvalue

# Positioning and labels for figure y axis
labs_era = []
for i in range(len(pressures_era)):
	labs_era.append(str(pressures_era[i]))
pos_era = np.arange(len(pressures_era))

# Colour depending on if positive or negative
cols_era = []
for item in slopes_era:
	if item>0:
		cols_era.append('r')
	else:
		cols_era.append('b')

# =================================== JRA =============================================
# Set working directory
jra_dir = ""
os.chdir(jra_dir)

pressures_jra = [200,250,300,400,500,600,700,850,925,1000]

# List to store files sorted by year
file_list = np.chararray((len(years),12),itemsize=51)

# Get all filenames
directory = os.getcwd()
for filename in sorted(os.listdir(directory)):
	year = int(filename[17:21])
	count = np.where(years==year)[0][0]
	month = int(filename[21:23])-1 
	file_list[count,month] = filename

# Array to store annual-mean zonal winds
ann_mean_gradient = np.zeros((len(years),len(pressures_jra)))

# Loop over files by year & month
count = 0 # Index for year
for i in range(np.shape(file_list)[0]):
	print years[i]
	monthly_temp = np.zeros((np.shape(file_list)[1],len(pressures_jra)))
	for j in range(np.shape(file_list)[1]):
		# Open file
		fname = file_list[i,j]
		print fname
		f = Dataset(fname,mode='r')
		
		# Grid
		lats_jra = f.variables['g0_lat_2'][:]
		lons_jra = f.variables['g0_lon_3'][:]
		levs = f.variables['lv_ISBL1'][:]
		
		levidx1 = np.where(levs==200)[0][0]
		levidx2 = np.where(levs==1000)[0][0]

		# Calculate two-boxes
		# Bounding boxes
		lat_polar1 = np.where(lats_jra==70)[0][0]
		lat_polar2 = np.where(lats_jra==50)[0][0]
		lat_tropic1 = np.where(lats_jra==50)[0][0]
		lat_tropic2 = np.where(lats_jra==30)[0][0]

		lats_polar = lats_jra[lat_polar1:lat_polar2+1]
		lats_tropic = lats_jra[lat_tropic1:lat_tropic2+1]

		# Extract
		polarbox = f.variables['TMP_GDS0_ISBL'][:,levidx1:levidx2+1,lat_polar1:lat_polar2+1,:]
		tropicbox = f.variables['TMP_GDS0_ISBL'][:,levidx1:levidx2+1,lat_tropic1:lat_tropic2+1,:]
		
		# Zonally average
		polar_zm = np.mean(polarbox,axis=-1)
		tropic_zm = np.mean(tropicbox,axis=-1)

		# Weight for latitude
		weights_polar = np.cos(np.radians(lats_polar))
		weights_tropic = np.cos(np.radians(lats_tropic))
		polar_mean = np.average(polar_zm,axis=-1,weights=weights_polar)
		tropic_mean = np.average(tropic_zm,axis=-1,weights=weights_tropic)		

		# Time-mean
		polar_month = np.mean(polar_mean,axis=0)
		tropic_month = np.mean(tropic_mean,axis=0)

		# Difference
		temp_month = tropic_month-polar_month
		
		# Store value
		monthly_temp[j,:] = temp_month		

	# Calculate annual-mean
	ann_mean_temp = np.mean(monthly_temp,axis=0)

	ann_mean_gradient[count,:] = ann_mean_temp
	
	f.close()

	count+=1

# Statistical analysis
# Array to store slopes
slopes_jra = np.zeros((len(pressures_jra)))
# Array to store slope error
confidence_jra = np.zeros(np.shape(slopes_jra))
# Array to store p values
pvals_jra = np.zeros(np.shape(slopes_jra))

for p in range(np.shape(ann_mean_gradient)[1]):
	slope,intercept,rvalue,pvalue,stderr = stats.linregress(years,ann_mean_gradient[:,p])
	slopes_jra[p] = 10*slope # trend per decade, -1 to flip lats
	confidence_jra[p] = 1.96*10*stderr # 95% confidence on slope per decade
	pvals_jra[p] = pvalue

# Positioning and labels for figure y axis
labs_jra = []
for i in range(len(pressures_jra)):
	labs_jra.append(str(pressures_jra[i]))
pos_jra = np.arange(len(pressures_jra))

# Colour depending on if positive or negative
cols_jra = []
for item in slopes_jra:
	if item>0:
		cols_jra.append('r')
	else:
		cols_jra.append('b')

#==================================== NCEP ============================================
ncep_dir = ""
os.chdir("")

pressures_ncep = [ 1000,   925,   850,   700,   600,   500,   400,   300, 250,   200]

# Array to store annual-mean temperature
ann_mean_gradient = np.zeros((len(years),len(pressures_ncep)))

count=0
for year in years:
	fname = "air."+str(year)+".nc"
	print fname

	f = Dataset(fname,mode='r')
	lats = f.variables['lat'][:]
	lons = f.variables['lon'][:]
	levs = f.variables['level'][:]

	# Bounding boxes
	lat_polar1 = np.where(lats==70)[0][0]
	lat_polar2 = np.where(lats==50)[0][0]
	lat_tropic1 = np.where(lats==50)[0][0]
	lat_tropic2 = np.where(lats==30)[0][0]

	lats_polar = lats[lat_polar1:lat_polar2+1]
	lats_tropic = lats[lat_tropic1:lat_tropic2+1]
	
	lonidx1 = np.where(lons==280)[0][0]
	lonidx2 = np.where(lons==350)[0][0]

	# Lats and lons of Atlantic sector
	lons_ncep = lons[lonidx1:lonidx2+1]

	levelidx1 = np.where(levs==1000)[0][0]
	levelidx2 = np.where(levs==200)[0][0]
	
	# Temperature
	polarbox = f.variables['air'][:,levelidx1:levelidx2+1,lat_polar1:lat_polar2+1,lonidx1:lonidx2+1]
	tropicbox = f.variables['air'][:,levelidx1:levelidx2+1,lat_tropic1:lat_tropic2+1,lonidx1:lonidx2+1]

	# Area average
	tropic_zm = np.mean(tropicbox,axis=-1)
	polar_zm = np.mean(polarbox,axis=-1)

	weights_tropic = np.cos(np.radians(lats_tropic))
	weights_polar = np.cos(np.radians(lats_polar))

	tropic_mean = np.average(tropic_zm,weights=weights_tropic,axis=-1)
	polar_mean = np.average(polar_zm,weights=weights_polar,axis=-1)

	diff = tropic_mean - polar_mean

	# Annual mean
	ann_mean_diff = np.mean(diff,axis=0)
	ann_mean_gradient[count,:] = ann_mean_diff
	
	f.close()

	count+=1

# Statistical analysis
# Array to store slopes
slopes_ncep = np.zeros((len(pressures_ncep)))
# Array to store slope error
confidence_ncep = np.zeros(np.shape(slopes_ncep))
# Array to store p values
pvals_ncep = np.zeros(np.shape(slopes_ncep))

for p in range(np.shape(ann_mean_gradient)[1]):
	slope,intercept,rvalue,pvalue,stderr = stats.linregress(years,ann_mean_gradient[:,p])
	slopes_ncep[p] = 10*slope # trend per decade, -1 to flip lats
	confidence_ncep[p] = 1.96*10*stderr # 95% confidence on slope per decade
	pvals_ncep[p] = pvalue

# Positioning and labels for figure y axis
labs_ncep = []
for i in range(len(pressures_ncep)):
	labs_ncep.append(str(pressures_ncep[i]))
pos_ncep = np.arange(len(pressures_ncep))

# Colour depending on if positive or negative
cols_ncep = []
for item in slopes_ncep:
	if item>0:
		cols_ncep.append('r')
	else:
		cols_ncep.append('b')

# ======================================== Figure ==================================================
plt.rcParams.update({'font.size': 7})
fig, (ax1,ax2,ax3) = plt.subplots(3,1,figsize=(3.5,9.7),sharex=True)
cpsize=2
ax1.barh(pos_era,slopes_era,edgecolor ="None",align='center',color=cols_era,xerr = confidence_era,ecolor='black',capsize=cpsize)
ax1.set_yticks(pos_era)
ax1.set_yticklabels(labs_era)
ax1.set_ylim(-1,len(pos_era))
ax1.axvline(x=0,color='k')
ax1.invert_yaxis()
ax1.grid(linestyle='--')
ax1.set_ylabel("pressure (hPa)")
ax1.set_xlim(-0.6,0.4)
ax1.text(0, 1.022, 'a',verticalalignment='bottom', horizontalalignment='right',transform=ax1.transAxes,color='k', weight='bold',fontsize=8)
ax1.text(0.5, 1.022, 'ERA-Interim',verticalalignment='bottom', horizontalalignment='center',transform=ax1.transAxes,color='k', weight='bold')

ax2.barh(pos_ncep,slopes_ncep,edgecolor ="None",align='center',color=cols_ncep,xerr = confidence_ncep,ecolor='black',capsize=cpsize)
ax2.set_yticks(pos_ncep)
ax2.set_yticklabels(labs_ncep)
ax2.set_ylim(-1,len(pos_ncep))
ax2.axvline(x=0,color='k')
ax2.grid(linestyle='--')
ax2.set_ylabel("pressure (hPa)")
ax2.set_xlim(-0.6,0.4)
ax2.text(0, 1.022, 'b',verticalalignment='bottom', horizontalalignment='right',transform=ax2.transAxes,color='k', weight='bold',fontsize=8)
ax2.text(0.5, 1.022, 'NCEP/NCAR',verticalalignment='bottom', horizontalalignment='center',transform=ax2.transAxes,color='k', weight='bold')

ax3.barh(pos_jra,slopes_jra,edgecolor ="None",align='center',color=cols_jra,xerr = confidence_jra,ecolor='black',capsize=cpsize)
ax3.set_yticks(pos_jra)
ax3.set_yticklabels(labs_jra)
ax3.set_ylim(-1,len(pos_jra))
ax3.axvline(x=0,color='k')
ax3.grid(linestyle='--')
ax3.invert_yaxis()
ax3.set_ylabel("pressure (hPa)")
ax3.set_xlabel("trend (K decade$^{-1}$)")
ax3.set_xlim(-0.6,0.4)
ax3.text(0, 1.022, 'c',verticalalignment='bottom', horizontalalignment='right',transform=ax3.transAxes,color='k', weight='bold',fontsize=8)
ax3.text(0.5, 1.022, 'JRA-55',verticalalignment='bottom', horizontalalignment='center',transform=ax3.transAxes,color='k', weight='bold')

plt.savefig("figure_2.png",dpi=400,bbox_inches="tight")
plt.savefig("figure_2.eps",dpi=400,bbox_inches="tight")
