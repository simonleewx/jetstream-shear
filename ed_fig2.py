import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from netCDF4 import Dataset
import os
from scipy import stats

# List of years
years = np.arange(1979,2018,1)

# ============================ NCEP/NCAR ==============================================
os.chdir(ncep_dir)

# Array to store annual-mean values
shearmax_annual_ncep = np.zeros((len(years)))

# Loop over years
for count,year in enumerate(years):	
	fname = "uwnd."+str(year)+".nc"
	print fname

	f = Dataset(fname,mode='r')
	lats = f.variables['lat'][:]
	lons = f.variables['lon'][:]
	levs = f.variables['level'][:]

	# Required levels
	levidx1 = np.where(levs==300)[0][0]
	levidx2 = np.where(levs==200)[0][0]

	# Bounding box
	latidx1 = np.where(lats==70)[0][0]
	latidx2 = np.where(lats==30)[0][0]
	lonidx1 = np.where(lons==280)[0][0]
	lonidx2 = np.where(lons==357.5)[0][0]

	latbox = lats[latidx1:latidx2+1]

	# Extract u-wind
	uwnd = f.variables['uwnd'][:,levidx1:levidx2+1,latidx1:latidx2+1,lonidx1:lonidx2+1]

	# Calculate 250 hPa shear in m/s/100 hPa
	shear = 100*((uwnd[:,2,:,:] - uwnd[:,0,:,:])/100.)

	# Find maximum value per 6 hours
	day_shearmax=[]
	for i in range(np.shape(shear)[0]):
		# The maximum shear
		shearmax = np.max(shear[i,:,:])
		# Store values
		day_shearmax.append(shearmax)

	# Take annual mean
	shearmax_annual = np.mean(day_shearmax)

	# Store 
	shearmax_annual_ncep[count] = shearmax_annual

	f.close()

# Statistical analysis
slope_ncep,intercept_ncep,rvalue_ncep,pvalue_ncep,stderr_ncep = stats.linregress(years,shearmax_annual_ncep)

# Percent change in shearmax
percent_shear_ncep = 100*(((slope_ncep*years+intercept_ncep)[-1]-(slope_ncep*years+intercept_ncep)[0])/(slope_ncep*years+intercept_ncep)[0])


#================================= ERA-Interim ======================================
# Set working directory
os.chdir(era_dir)

# Get all filenames
directory = os.getcwd()

# Array to store winds
hourly_winds = np.zeros((len(years),366*4))
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

			levidx1 = np.where(levs==200)[0][0]
			levidx2 = np.where(levs==300)[0][0]

			# Extract u-wind
			uwnd = f.variables['U_GDS4_ISBL'][:,levidx1:levidx2+1,:,:].squeeze()

			# Calculate 250 hPa shear
			shear = 100*((uwnd[0,:,:]-uwnd[2,:,:])/100.)

			# Find maximum value
			shearmax = np.max(shear)

			# Store value
			hourly_winds[np.where(years==year)[0][0],doy] = shearmax

			f.close()
			doy+=1

# Calculate annual mean
shearmax_annual_era = np.nanmean(hourly_winds,axis=1)

# Statistical analysis
slope_era,intercept_era,rvalue_era,pvalue_era,stderr_era = stats.linregress(years,shearmax_annual_era)
# Percent change
percent_shear_era = 100*(((slope_era*years+intercept_era)[-1]-(slope_era*years+intercept_era)[0])/(slope_era*years+intercept_era)[0])


#======================================= JRA =========================================================
# Set working directory
os.chdir(jra_dir)

# List to store files sorted by year
file_list = np.chararray((len(years),12),itemsize=52)

# Get all filenames
directory = os.getcwd()
for filename in sorted(os.listdir(directory)):
	year = int(filename[18:22])
	count = np.where(years==year)[0][0]
	month = int(filename[22:24])-1 
	file_list[count,month] = filename

# Array to store annual-means
shearmax_annual_jra = np.zeros((len(years)))

# Loop over files by year & month
for count,i in enumerate(range(np.shape(file_list)[0])):
	print years[i]
	monthly_wind = np.zeros((np.shape(file_list)[1]))
	for j in range(np.shape(file_list)[1]):
		# Open file
		fname = file_list[i,j]
		print fname
		f = Dataset(fname,mode='r')
		
		# Grid
		lats = f.variables['g0_lat_2'][:]
		levs = f.variables['lv_ISBL1'][:]

		levidx1 = np.where(levs==200)[0][0]
		levidx2 = np.where(levs==300)[0][0]

		# Extract u-wind
		uwnd = f.variables['UGRD_GDS0_ISBL'][:,levidx1:levidx2+1,:,:]

		# Calculate shear
		shear = 100*((uwnd[:,0,:,:]-uwnd[:,2,:,:])/100.)
		
		# Find maximum value per 6 hours
		day_shearmax=[]
		for k in range(np.shape(shear)[0]):
			# The maximum shear
			shearmax = np.max(shear[k,:,:])
			# Store values
			day_shearmax.append(shearmax)

		# Monthly mean		
		monthly_wind[j] = np.mean(day_shearmax)

		f.close()

	# Calculate annual-mean
	annual_wind = np.mean(monthly_wind)

	# Store value
	shearmax_annual_jra[count] = annual_wind

# Statistical analysis
slope_jra,intercept_jra,rvalue_jra,pvalue_jra,stderr_jra = stats.linregress(years,shearmax_annual_jra)
# Percent change
percent_shear_jra = 100*(((slope_jra*years+intercept_jra)[-1]-(slope_jra*years+intercept_jra)[0])/(slope_jra*years+intercept_jra)[0])

# ========================== Ensemble mean ==============================
ensemble_shearmax = (shearmax_annual_jra+shearmax_annual_ncep+shearmax_annual_era)/3.
slope_ens,intercept_ens,rvalue_ens,pvalue_ens,stderr_ens = stats.linregress(years,ensemble_shearmax)
percent_shear_ens = 100*(((slope_ens*years+intercept_ens)[-1]-(slope_ens*years+intercept_ens)[0])/(slope_ens*years+intercept_ens)[0])

# ========================== Figure =====================================
plt.rcParams.update({'font.size': 7})
plt.figure(figsize=(4.7,3.5))
plt.plot(years,shearmax_annual_era,label="ERA-Interim",color='blue',alpha=0.7)
plt.plot(years,shearmax_annual_ncep,label="NCEP/NCAR",color='orange',alpha=0.7)
plt.plot(years,shearmax_annual_jra,label="JRA-55",color='green',alpha=0.7)
plt.plot(years,ensemble_shearmax,label="mean",marker='o',linewidth=2,color='k')
plt.plot(years,slope_ens*years+intercept_ens,linestyle='--',color='k',linewidth=2,label="mean trend")
plt.xlim(1978,2018)
plt.ylabel("vertical zonal shear (m s$^{-1}$ (100 hPa)$^{-1}$)")
plt.xlabel("year")
plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15),frameon=False, ncol=5)
plt.savefig("ed_fig2.png",bbox_inches="tight",dpi=400)
plt.savefig("ed_fig2.eps",bbox_inches="tight",dpi=400)
