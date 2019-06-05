import numpy as np
from netCDF4 import Dataset,num2date
import os
import matplotlib.pyplot as plt
from scipy import stats

# Years to analyse over
years = np.arange(1979,2018,1)

# ================================================ ERA ===========================================================
# Set working directory
os.chdir(era_dir)

# Pressure levels
pressure = [ 150, 200, 250,  300,  350,  400,  450,  500,  550,  600, 650,  700,  750,  800,  850,  900,  950, 1000]

# Number of pressure levels
npres = len(pressure)

# Get all filenames
directory = os.getcwd()

# Array to store winds
hourly_winds = np.zeros((len(years),366*4,npres))
hourly_winds[:] = np.nan 

# Loop over all files
for year in years:
	doy=0
	for filename in sorted(os.listdir(directory)):
		if str(year) in filename[24:28]:
			print filename
			# Open file
			f = Dataset(filename,mode='r')
		
			# Grid
			lats = f.variables['g4_lat_2'][:]
			levs = f.variables['lv_ISBL1'][:]

			# Weights for area-mean
			weights = np.cos(np.radians(lats))

			# Extract u-wind
			uwnd = f.variables['U_GDS4_ISBL'][:,1::,:,:]
			# Zonally average
			uwnd_zm = np.mean(uwnd,axis=-1)
			# Area-weighted mean
			uwnd_mean = np.average(uwnd_zm,axis=-1,weights=weights)

			# Store value
			hourly_winds[np.where(years==year)[0][0],doy,:] = uwnd_mean

			f.close()
			doy+=1

			
# Calculate annual-mean
list_winds_era = np.nanmean(hourly_winds,axis=1)

# Calculate shear by finite centred differences
list_shear_era = np.zeros((np.shape(list_winds_era)[0],len(pressure)-2))
for lev in range(len(pressure)-2):
	shear = (list_winds_era[:,lev]-list_winds_era[:,lev+2])/(pressure[lev+2]-pressure[lev])
	list_shear_era[:,lev] = shear

# Trends in winds
slopes_winds_era = []
confidence_winds_era = []
pvalues_winds_era = []
for lev in range(len(pressure)):
	slope,intercept,rvalue,pvalue,stderr = stats.linregress(years,list_winds_era[:,lev])
	slopes_winds_era.append(slope*10)
	confidence_winds_era.append(stderr*10*1.96)
	pvalues_winds_era.append(pvalue)

# Trends in shear
slopes_shear_era = []
confidence_shear_era = []
pvalues_shear_era = []
for lev in range(len(pressure)-2):
	slope,intercept,rvalue,pvalue,stderr = stats.linregress(years,list_shear_era[:,lev])
	slopes_shear_era.append(slope*10*100)
	confidence_shear_era.append(stderr*10*1.96*100)
	pvalues_shear_era.append(pvalue)

# ============================================= NCEP/NCAR ============================================
os.chdir(ncep_dir)

# Pressure levels
pressure = [ 1000,   925,   850,   700,   600,   500,   400,   300, 250,   200,   150]

# Lists to store annual-mean values
list_winds_ncep = np.zeros((len(years),len(pressure)))

# Loop over years
for count,year in enumerate(years):	
	fname = "uwnd."+str(year)+".nc"
	print fname

	f = Dataset(fname,mode='r')
	lats = f.variables['lat'][:]
	lons = f.variables['lon'][:]
	levs = f.variables['level'][:]

	# Required levels
	levidx1 = np.where(levs==1000)[0][0]
	levidx2 = np.where(levs==150)[0][0]

	# Bounding box
	latidx1 = np.where(lats==70)[0][0]
	latidx2 = np.where(lats==30)[0][0]
	lonidx1 = np.where(lons==280)[0][0]
	lonidx2 = np.where(lons==350)[0][0]

	# Weights for area-mean
	weights = np.cos(np.radians(lats[latidx1:latidx2+1]))

	# Extract u-wind
	uwnd = f.variables['uwnd'][:,levidx1:levidx2+1,latidx1:latidx2+1,lonidx1:lonidx2+1]
	uwnd_zm = np.mean(uwnd,axis=-1)
	uwnd_mean = np.average(uwnd_zm,axis=-1,weights=weights)

	# Annual-mean
	uwnd_annual_mean = np.mean(uwnd_mean,axis=0)
	
	# Store in array
	list_winds_ncep[count,:] = uwnd_annual_mean

	f.close()

# Calculate shear
list_shear_ncep = np.zeros((np.shape(list_winds_ncep)[0],np.shape(list_winds_ncep)[1]-2))
for level in range(len(pressure)-2):
	list_shear_ncep[:,level] = (list_winds_ncep[:,level+2] - list_winds_ncep[:,level])/(pressure[level]-pressure[level+2])

# Trends in winds
slopes_winds_ncep = []
confidence_winds_ncep = []
pvalues_winds_ncep = []
for lev in range(len(pressure)):
	slope,intercept,rvalue,pvalue,stderr = stats.linregress(years,list_winds_ncep[:,lev])
	slopes_winds_ncep.append(slope*10)
	confidence_winds_ncep.append(stderr*10*1.96)
	pvalues_winds_ncep.append(pvalue)

# Trends in shear
slopes_shear_ncep = []
confidence_shear_ncep = []
pvalues_shear_ncep = []
for lev in range(len(pressure)-2):
	slope,intercept,rvalue,pvalue,stderr = stats.linregress(years,list_shear_ncep[:,lev])
	slopes_shear_ncep.append(slope*10*100)
	confidence_shear_ncep.append(stderr*100*10*1.96)
	pvalues_shear_ncep.append(pvalue)

# ================================================ JRA =================================================================
# Set working directory
os.chdir(jra_dir)

# Pressure levels
pressure_jra = [ 150, 200, 250,  300,  350,  400,  450,  500,  550,  600, 650,  700,  750,  800,  850,  900,  950, 1000]

# Number of pressure levels
npres = len(pressure_jra)

# List to store files sorted by year
file_list = np.chararray((len(years),12),itemsize=52)

# Get all filenames
directory = os.getcwd()
for filename in sorted(os.listdir(directory)):
	year = int(filename[18:22])
	count = np.where(years==year)[0][0]
	month = int(filename[22:24])-1 
	file_list[count,month] = filename

# List to store annual-mean zonal winds
list_winds_jra = np.zeros((len(years),npres))

# Loop over files by year & month
for count,i in enumerate(range(np.shape(file_list)[0])):
	print years[i]
	monthly_wind = np.zeros((np.shape(file_list)[1],npres))
	for j in range(np.shape(file_list)[1]):
		# Open file
		fname = file_list[i,j]
		print fname
		f = Dataset(fname,mode='r')
		
		# Grid
		lats = f.variables['g0_lat_2'][:]
		levs = f.variables['lv_ISBL1'][:]

		# Weights for area-mean
		weights = np.cos(np.radians(lats))

		# Extract u-wind
		uwnd = f.variables['UGRD_GDS0_ISBL'][:,1::,:,:]
		# Zonally average
		uwnd_zm = np.mean(uwnd,axis=-1)
		# Area-weighted mean
		uwnd_mean = np.average(uwnd_zm,axis=-1,weights=weights)
		# Time-mean
		uwnd_month = np.mean(uwnd_mean,axis=0)
		
		# Store value
		monthly_wind[j,:] = uwnd_month

		f.close()

	# Calculate annual-mean
	annual_wind = np.mean(monthly_wind,axis=0)
	# Store value
	list_winds_jra[count,:] = annual_wind

# Calculate shear by finite centred differences
list_shear_jra = np.zeros((np.shape(list_winds_jra)[0],len(pressure_jra)-2))
for lev in range(len(pressure_jra)-2):
	shear = (list_winds_jra[:,lev]-list_winds_jra[:,lev+2])/(pressure_jra[lev+2]-pressure_jra[lev])
	list_shear_jra[:,lev] = shear

# Trends in winds
slopes_winds_jra = []
confidence_winds_jra = []
for lev in range(len(pressure_jra)):
	slope,intercept,rvalue,pvalue,stderr = stats.linregress(years,list_winds_jra[:,lev])
	slopes_winds_jra.append(slope*10)
	confidence_winds_jra.append(stderr*10*1.96)

# Trends in shear
slopes_shear_jra = []
confidence_shear_jra = []
for lev in range(len(pressure_jra)-2):
	slope,intercept,rvalue,pvalue,stderr = stats.linregress(years,list_shear_jra[:,lev])
	slopes_shear_jra.append(slope*10*100)
	confidence_shear_jra.append(stderr*10*1.96*100)

# ================================================ FIGURE ===============================================
# Pressure levels in JRA & ERA
pressure = [ 150, 200, 250,  300,  350,  400,  450,  500,  550,  600, 650,  700,  750,  800,  850,  900,  950, 1000]
shear_labs = []
for i in range(1,len(pressure)-1):
	shear_labs.append(str(pressure[i]))
shearpos = np.arange(len(pressure)-2)
wind_labs = []
for item in pressure:
	wind_labs.append(str(item))
windpos = np.arange(len(pressure))

# Pressure levels in NCEP
pressure_ncep = [ 1000,   925,   850,   700,   600,   500,   400,   300, 250,   200,   150]
shear_labs_ncep = []
for i in range(1,len(pressure_ncep)-1):
	shear_labs_ncep.append(str(pressure_ncep[i]))
shearpos_ncep = np.arange(len(pressure_ncep)-2)
wind_labs_ncep = []
for item in pressure_ncep:
	wind_labs_ncep.append(str(item))
windpos_ncep = np.arange(len(pressure_ncep))

# Colour depending on pos or neg
era_cols_wind = []
jra_cols_wind = []
ncep_cols_wind = []
for item in slopes_winds_era:
	if item>0:
		era_cols_wind.append('r')
	else:
		era_cols_wind.append('b')
for item in slopes_winds_jra:
	if item>0:
		jra_cols_wind.append('r')
	else:
		jra_cols_wind.append('b')
for item in slopes_winds_ncep:
	if item>0:
		ncep_cols_wind.append('r')
	else:
		ncep_cols_wind.append('b')
era_cols_shear = []
jra_cols_shear = []
ncep_cols_shear = []
for item in slopes_shear_era:
	if item>0:
		era_cols_shear.append('r')
	else:
		era_cols_shear.append('b')
for item in slopes_shear_jra:
	if item>0:
		jra_cols_shear.append('r')
	else:
		jra_cols_shear.append('b')
for item in slopes_shear_ncep:
	if item>0:
		ncep_cols_shear.append('r')
	else:
		ncep_cols_shear.append('b')

wind_labs = wind_labs[::2]
shear_labs = shear_labs[::2]

shearpos_old = shearpos

shearpos = np.arange(1,17,1)
shearpos_ncep = np.arange(1,10,1)

# Make figure
plt.rcParams.update({'font.size': 7})
fig, ((ax1, ax2), (ax3, ax4), (ax5,ax6)) = plt.subplots(3, 2, sharex='col', sharey=False,figsize=(7.2,9.7))
cpsize=2
# ERA
ax1.barh(shearpos,slopes_shear_era,edgecolor ="None",align='center',color=era_cols_shear,xerr = confidence_shear_era,ecolor='black',capsize=cpsize)
ax1.set_yticks(shearpos[::2])
ax1.set_yticklabels(shear_labs)
ax1.set_ylim(-1,len(shearpos)+2)
ax1.invert_yaxis()
ax1.axvline(x=0,color='k')
ax1.grid(linestyle='--')
ax1.set_xlim(-0.1,0.2)
ax1.set_ylabel("pressure (hPa)")
ax1.text(0, 1.022, 'a',verticalalignment='bottom', horizontalalignment='right',transform=ax1.transAxes,color='k', weight='bold')
ax1.text(0.5, 1.022, 'ERA-Interim',verticalalignment='bottom', horizontalalignment='center',transform=ax1.transAxes,color='k', weight='bold')

ax2.barh(windpos,slopes_winds_era,edgecolor ="None",align='center',color=era_cols_wind,xerr = confidence_winds_era,ecolor='black',capsize=cpsize)
ax2.set_yticks(windpos[::2])
ax2.set_yticklabels(wind_labs)
ax2.set_ylim(-1,len(windpos))
ax2.invert_yaxis()
ax2.grid(linestyle='--')
ax2.axvline(x=0,color='k')
ax2.set_xlim(-0.25,0.25)
ax2.text(0, 1.022, 'd',verticalalignment='bottom', horizontalalignment='right',transform=ax2.transAxes,color='k', weight='bold')
ax2.text(0.5, 1.022, 'ERA-Interim',verticalalignment='bottom', horizontalalignment='center',transform=ax2.transAxes,color='k', weight='bold')

# JRA
ax5.barh(shearpos,slopes_shear_jra,edgecolor ="None",align='center',color=jra_cols_shear,xerr = confidence_shear_jra,ecolor='black',capsize=cpsize)
ax5.set_yticks(shearpos[::2])
ax5.set_yticklabels(shear_labs)
ax5.set_ylim(-1,len(shearpos)+2)
ax5.invert_yaxis()
ax5.axvline(x=0,color='k')
ax5.grid(linestyle='--')
ax5.set_xlim(-0.1,0.2)
ax5.set_ylabel("pressure (hPa)")
ax5.text(0, 1.022, 'c',verticalalignment='bottom', horizontalalignment='right',transform=ax5.transAxes,color='k', weight='bold')
ax5.text(0.5, 1.022, 'JRA-55',verticalalignment='bottom', horizontalalignment='center',transform=ax5.transAxes,color='k', weight='bold')
ax5.set_xlabel("trend (m s$^{-1}$ 100 hPa$^{-1}$ decade$^{-1}$)")

ax6.barh(windpos,slopes_winds_jra,edgecolor ="None",align='center',color=jra_cols_wind,xerr = confidence_winds_jra,ecolor='black',capsize=cpsize)
ax6.set_yticks(windpos[::2])
ax6.set_yticklabels(wind_labs)
ax6.set_ylim(-1,len(windpos))
ax6.invert_yaxis()
ax6.grid(linestyle='--')
ax6.axvline(x=0,color='k')
ax6.set_xlim(-0.25,0.25)
ax6.text(0, 1.022, 'f',verticalalignment='bottom', horizontalalignment='right',transform=ax6.transAxes,color='k', weight='bold')
ax6.text(0.5, 1.022, 'JRA-55',verticalalignment='bottom', horizontalalignment='center',transform=ax6.transAxes,color='k', weight='bold')
ax6.set_xlabel("trend (m s$^{-1}$ decade$^{-1}$)")

# NCEP
ax3.barh(shearpos_ncep,slopes_shear_ncep,edgecolor ="None",align='center',color=ncep_cols_shear,xerr = confidence_shear_ncep,ecolor='black',capsize=cpsize)
ax3.set_yticks(shearpos_ncep)
ax3.set_yticklabels(shear_labs_ncep)
ax3.set_ylim(-1,len(shearpos_ncep)+2)
ax3.axvline(x=0,color='k')
ax3.grid(linestyle='--')
ax3.set_xlim(-0.1,0.2)
ax3.set_ylabel("pressure (hPa)")
ax3.text(0, 1.022, 'b',verticalalignment='bottom', horizontalalignment='right',transform=ax3.transAxes,color='k', weight='bold')
ax3.text(0.5, 1.022, 'NCEP/NCAR',verticalalignment='bottom', horizontalalignment='center',transform=ax3.transAxes,color='k', weight='bold')

ax4.barh(windpos_ncep,slopes_winds_ncep,edgecolor ="None",align='center',color=ncep_cols_wind,xerr = confidence_winds_ncep,ecolor='black',capsize=cpsize)
ax4.set_yticks(windpos_ncep)
ax4.set_yticklabels(wind_labs_ncep)
ax4.set_ylim(-1,len(windpos_ncep))
ax4.grid(linestyle='--')
ax4.axvline(x=0,color='k')
ax4.set_xlim(-0.25,0.25)
ax4.text(0, 1.022, 'e',verticalalignment='bottom', horizontalalignment='right',transform=ax4.transAxes,color='k', weight='bold')
ax4.text(0.5, 1.022, 'NCEP/NCAR',verticalalignment='bottom', horizontalalignment='center',transform=ax4.transAxes,color='k', weight='bold')

plt.savefig("ed_fig1.png",bbox_inches="tight",dpi=400)
plt.savefig("ed_fig1.eps",bbox_inches="tight",dpi=400)
