"""
Produces Figure 1 from Lee, Williams, and Frame 2019.
Python 2.
NCEP, ERA, JRA data.
"""
import numpy as np
from netCDF4 import Dataset
import os
from scipy import stats
from mpl_toolkits.basemap import Basemap 
import matplotlib.pyplot as plt

years = np.arange(1979,2018,1)

# Set pressure level to analyse on
level = 250.

#===================================== NCEP ============================================

# working directory for ncep data
ncep_dir = ""
os.chdir(ncep_dir)

# Array to store annual-mean temperature
ann_mean_temps = np.zeros((len(years),17,32))

for count,year in enumerate(years):
	fname = "air."+str(year)+".nc"
	print fname

	f = Dataset(fname,mode='r')
	lats = f.variables['lat'][:]
	lons = f.variables['lon'][:]
	levs = f.variables['level'][:]

	# Bounding box
	latidx1 = np.where(lats==70)[0][0]
	latidx2 = np.where(lats==30)[0][0]
	lonidx1 = np.where(lons==280)[0][0]
	lonidx2 = np.where(lons==357.5)[0][0]

	levelidx = np.where(levs==level)[0][0]

	
	# Temperature
	temperature = f.variables['air'][:,levelidx,latidx1:latidx2+1,lonidx1:lonidx2+1]

	# Annual mean
	ann_mean_temp = np.mean(temperature,axis=0)

	# Store annual mean
	ann_mean_temps[count,:,:]=ann_mean_temp
	
	f.close()

# Lats and lons of Atlantic sector
lons_ncep = lons[lonidx1:lonidx2+1]
lats_ncep = lats[latidx1:latidx2+1]	

# Statistical analysis
# Array to store slopes
slopes_ncep = np.zeros((np.shape(ann_mean_temps)[1],np.shape(ann_mean_temps)[2]))
# Array to store p values
pvals_ncep = np.zeros(np.shape(slopes_ncep))

for i in range(np.shape(ann_mean_temps)[1]):
	for j in range(np.shape(ann_mean_temps)[2]):
		slope,intercept,rvalue,pvalue,stderr = stats.linregress(years,ann_mean_temps[:,i,j])
		slopes_ncep[i,j] = slope
		pvals_ncep[i,j] = pvalue

# Mask non-sig regions
pvals_sig_ncep = np.ma.masked_greater(pvals_ncep,0.05)

# Figure
# Convert slopes to per decade
decade_slopes_ncep = 10*slopes_ncep

#================================================== ERA ===============================================
# working directory for era data
era_dir=""
os.chdir(era_dir)

# List to store files sorted by year
file_list = np.chararray((len(years),366),itemsize=47)

# Get all filenames
directory = os.getcwd()

# Array to store 6 hourly temps
hourly_temps = np.zeros((len(years),366*4,57,99))
hourly_temps[:] = np.nan 

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

			levidx = np.where(levs==level)[0][0]

			# Extract temperature
			temperature = f.variables['T_GDS4_ISBL'][:,levidx,:,:].squeeze()

			# Store value
			hourly_temps[np.where(years==year)[0][0],doy,:,:] = temperature

			f.close()
			doy+=1
		
# Calculate annual-mean
ann_mean_temps = np.nanmean(hourly_temps,axis=1)

# Statistical analysis
# Array to store slopes
slopes_era = np.zeros((np.shape(ann_mean_temps)[1],np.shape(ann_mean_temps)[2]))
# Array to store p values
pvals_era = np.zeros(np.shape(slopes_era))

for i in range(np.shape(ann_mean_temps)[1]):
	for j in range(np.shape(ann_mean_temps)[2]):
		slope,intercept,rvalue,pvalue,stderr = stats.linregress(years,ann_mean_temps[:,i,j])
		slopes_era[i,j] = slope
		pvals_era[i,j] = pvalue

# Mask non-sig regions
pvals_sig_era = np.ma.masked_greater(pvals_era,0.05)

# Figure
# Convert slopes to per decade
decade_slopes_era = 10*slopes_era


#================================================ JRA ==================================================
# Set working directory for jra data
jra_dir = ""
os.chdir(jra_dir)

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
ann_mean_temps = np.zeros((len(years),33,57))

# Loop over files by year & month
for count,i in enumerate(range(np.shape(file_list)[0])):
	print years[i]
	monthly_temp = np.zeros((np.shape(file_list)[1],33,57))
	for j in range(np.shape(file_list)[1]):
		# Open file
		fname = file_list[i,j]
		print fname
		f = Dataset(fname,mode='r')
		
		# Grid
		lats_jra = f.variables['g0_lat_2'][:]
		lons_jra = f.variables['g0_lon_3'][:]
		levs = f.variables['lv_ISBL1'][:]

		levidx = np.where(levs==level)[0][0]

		# Extract temperature
		temperature = f.variables['TMP_GDS0_ISBL'][:,levidx,:,:]

		# Time-mean
		temp_month = np.mean(temperature,axis=0)
		
		# Store value
		monthly_temp[j,:,:] = temp_month

		f.close()

	# Calculate annual-mean
	annual_temp = np.mean(monthly_temp,axis=0)
	# Store value
	ann_mean_temps[count,:,:] = annual_temp

# Statistical analysis
# Array to store slopes
slopes_jra = np.zeros((np.shape(ann_mean_temps)[1],np.shape(ann_mean_temps)[2]))
# Array to store p values
pvals_jra = np.zeros(np.shape(slopes_jra))

for i in range(np.shape(ann_mean_temps)[1]):
	for j in range(np.shape(ann_mean_temps)[2]):
		slope,intercept,rvalue,pvalue,stderr = stats.linregress(years,ann_mean_temps[:,i,j])
		slopes_jra[i,j] = slope
		pvals_jra[i,j] = pvalue

# Mask non-sig regions
pvals_sig_jra = np.ma.masked_greater(pvals_jra,0.05)

# Convert slopes to per decade
decade_slopes_jra = 10*slopes_jra

#=================================== FIGURE ==================================================================
# Set contouring levels
plt.rcParams.update({'font.size': 7})
temp_levs = np.arange(-0.3,0.32,0.02)

fig, (ax1,ax2,ax3) = plt.subplots(3,1,figsize=(3.5,9.7))

m = Basemap(projection='mill',llcrnrlat=30,urcrnrlat=70,llcrnrlon=280,urcrnrlon=350,resolution='l',ax=ax1)        
m.drawcoastlines(linewidth=1.25, color='#444444')     
m.drawmeridians(np.arange(0, 360, 10), labels=[0,0,0,1],linewidth=1,color='gray')
m.drawparallels(np.arange(0, 90, 10), labels=[1,0,0,0], linewidth=1,color='gray')   
x,y = np.meshgrid(lons_era,lats_era)
x1, y1 = m(x, y)
m.contourf(x1,y1,decade_slopes_era,temp_levs,cmap="bwr",extend="both")
m.pcolor(x1, y1, pvals_sig_era, hatch='.',alpha=0.,zorder=10)
ax1.text(0, 1.04, 'a',verticalalignment='bottom', horizontalalignment='right',transform=ax1.transAxes,color='k', weight='bold')
ax1.set_title('ERA-Interim',weight="bold")

m = Basemap(projection='mill',llcrnrlat=30,urcrnrlat=70,llcrnrlon=280,urcrnrlon=350,resolution='l',ax=ax2)        
m.drawcoastlines(linewidth=1.25, color='#444444')     
m.drawmeridians(np.arange(0, 360, 10), labels=[0,0,0,1],linewidth=1,color='gray')
m.drawparallels(np.arange(0, 90, 10), labels=[1,0,0,0], linewidth=1,color='gray')   
x,y = np.meshgrid(lons_ncep,lats_ncep)
x1, y1 = m(x, y)
m.contourf(x1,y1,decade_slopes_ncep,temp_levs,cmap="bwr",extend="both")
m.pcolor(x1, y1, pvals_sig_ncep, hatch='.',alpha=0.,zorder=10)
ax2.text(0, 1.04, 'b',verticalalignment='bottom', horizontalalignment='right',transform=ax2.transAxes,color='k', weight='bold')
ax2.set_title('NCEP/NCAR',weight="bold")

m = Basemap(projection='mill',llcrnrlat=30,urcrnrlat=70,llcrnrlon=280,urcrnrlon=350,resolution='l',ax=ax3)        
m.drawcoastlines(linewidth=1.25, color='#444444')     
m.drawmeridians(np.arange(0, 360, 10), labels=[0,0,0,1],linewidth=1,color='gray')
m.drawparallels(np.arange(0, 90, 10), labels=[1,0,0,0], linewidth=1,color='gray')   
x,y = np.meshgrid(lons_jra,lats_jra)
x1, y1 = m(x, y)
cf = m.contourf(x1,y1,decade_slopes_jra,temp_levs,cmap="bwr",extend="both")
m.pcolor(x1, y1, pvals_sig_jra,hatch='.',alpha=0.,zorder=10)
ax3.text(0, 1.04, 'c',verticalalignment='bottom', horizontalalignment='right',transform=ax3.transAxes,color='k', weight='bold')
ax3.set_title('JRA-55',weight="bold")

cb_ax = fig.add_axes([0.1,0.05,0.8,0.02])
cbar = fig.colorbar(cf, cax=cb_ax,ticks = np.arange(min(temp_levs),max(temp_levs)+0.1,0.1),orientation="horizontal")
cbar.set_label("trend (K decade$^{-1}$)")

plt.savefig("figure_1.png",dpi=400,bbox_inches="tight")
plt.savefig("figure_1.eps",dpi=400,bbox_inches="tight")