import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from netCDF4 import Dataset
import os
from scipy import stats

years = np.arange(1979,2018,1)

# ============================ NCEP/NCAR ==============================================
os.chdir(ncep_dir)

# Array to store annual-mean values
shear_annual_ncep = np.zeros((len(years)))
speed_annual_ncep = np.zeros((len(years)))

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

	# Weights for area mean
	weights = np.cos(np.radians(latbox))

	# Extract u-wind
	uwnd = f.variables['uwnd'][:,levidx1:levidx2+1,latidx1:latidx2+1,lonidx1:lonidx2+1]

	# Calculate 250 hPa shear in m/s/100 hPa
	shear = 100*((uwnd[:,2,:,:] - uwnd[:,0,:,:])/100.)

	# 250 hPa wind
	speed = uwnd[:,1,:,:]

	# Take zonal-mean
	shear_zm = np.mean(shear,axis=-1)
	speed_zm = np.mean(speed,axis=-1)

	# Take area-mean
	shear_mean = np.average(shear_zm,weights=weights,axis=-1)
	speed_mean = np.average(speed_zm,weights=weights,axis=-1)

	# Take annual mean
	shear_annual = np.mean(shear_mean)
	speed_annual = np.mean(speed_mean)

	# Store 
	shear_annual_ncep[count] = shear_annual
	speed_annual_ncep[count] = speed_annual

	f.close()

# Statistical analysis
slope_ncep,intercept_ncep,rvalue_ncep,pvalue_ncep,stderr_ncep = stats.linregress(years,shear_annual_ncep)
slope_ncep_speed,intercept_ncep_speed,rvalue_ncep_speed,pvalue_ncep_speed,stderr_ncep_speed = stats.linregress(years,speed_annual_ncep)

# Percent change in shearmax
percent_shear_ncep = 100*(((slope_ncep*years+intercept_ncep)[-1]-(slope_ncep*years+intercept_ncep)[0])/(slope_ncep*years+intercept_ncep)[0])


#================================= ERA-Interim ======================================
# Set working directory
os.chdir(era_dir)

# List to store files sorted by year
file_list = np.chararray((len(years),366),itemsize=47)

# Get all filenames
directory = os.getcwd()

# Array to store winds
hourly_winds = np.zeros((len(years),366*4))
hourly_winds[:] = np.nan 

hourly_speed = np.zeros((len(years),366*4))
hourly_speed[:] = np.nan 

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

			# Weights for area mean
			weights = np.cos(np.radians(lats))

			# Extract u-wind
			uwnd = f.variables['U_GDS4_ISBL'][:,levidx1:levidx2+1,:,:].squeeze()

			# Calculate 250 hPa shear
			shear = 100*((uwnd[0,:,:]-uwnd[2,:,:])/100.)

			# 250 hPa wind
			speed = uwnd[1,:,:]

			# Take zonal-mean
			shear_zm = np.mean(shear,axis=-1)
			speed_zm = np.mean(speed,axis=-1)

			# Take regional mean
			shear_mean = np.average(shear_zm,weights=weights,axis=-1)
			speed_mean = np.average(speed_zm,weights=weights,axis=-1)

			# Store value
			hourly_winds[np.where(years==year)[0][0],doy] = shear_mean
			hourly_speed[np.where(years==year)[0][0],doy] = speed_mean

			f.close()
			doy+=1

# Calculate annual mean
shear_annual_era = np.nanmean(hourly_winds,axis=1)
speed_annual_era = np.nanmean(hourly_speed,axis=1)

# Statistical analysis
slope_era,intercept_era,rvalue_era,pvalue_era,stderr_era = stats.linregress(years,shear_annual_era)
slope_era_speed,intercept_era_speed,rvalue_era_speed,pvalue_era_speed,stderr_era_speed = stats.linregress(years,speed_annual_era)

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
shear_annual_jra = np.zeros((len(years)))
speed_annual_jra = np.zeros((len(years)))

# Loop over files by year & month
for count,i in enumerate(range(np.shape(file_list)[0])):
	print years[i]
	monthly_shear = np.zeros((np.shape(file_list)[1]))
	monthly_speed = np.zeros((np.shape(file_list)[1]))
	for j in range(np.shape(file_list)[1]):
		# Open file
		fname = file_list[i,j]
		print fname
		f = Dataset(fname,mode='r')
		
		# Grid
		lats = f.variables['g0_lat_2'][:]
		levs = f.variables['lv_ISBL1'][:]

		weights=np.cos(np.radians(lats))

		levidx1 = np.where(levs==200)[0][0]
		levidx2 = np.where(levs==300)[0][0]

		# Extract u-wind
		uwnd = f.variables['UGRD_GDS0_ISBL'][:,levidx1:levidx2+1,:,:]

		# Calculate shear
		shear = 100*((uwnd[:,0,:,:]-uwnd[:,2,:,:])/100.)

		# 250 hPa wind
		speed = uwnd[:,1,:,:]

		# Zonal mean
		shear_zm = np.mean(shear,axis=-1)
		speed_zm = np.mean(speed,axis=-1)

		# Area mean weighted
		shear_mean = np.average(shear_zm,axis=-1,weights=weights)
		speed_mean = np.average(speed_zm,axis=-1,weights=weights)
		
		# Monthly mean		
		monthly_shear[j] = np.mean(shear_mean)
		monthly_speed[j] = np.mean(speed_mean)

		f.close()

	# Calculate annual-mean
	annual_shear = np.mean(monthly_shear)
	annual_speed = np.mean(monthly_speed)

	# Store value
	shear_annual_jra[count] = annual_shear
	speed_annual_jra[count] = annual_speed

# Statistical analysis
slope_jra,intercept_jra,rvalue_jra,pvalue_jra,stderr_jra = stats.linregress(years,shear_annual_jra)
slope_jra_speed,intercept_jra_speed,rvalue_jra_speed,pvalue_jra_speed,stderr_jra_speed = stats.linregress(years,speed_annual_jra)

# Percent change
percent_shear_jra = 100*(((slope_jra*years+intercept_jra)[-1]-(slope_jra*years+intercept_jra)[0])/(slope_jra*years+intercept_jra)[0])

# ========================== Ensemble mean ==============================
ensemble_shear = (shear_annual_jra+shear_annual_ncep+shear_annual_era)/3.

slope_ens,intercept_ens,rvalue_ens,pvalue_ens,stderr_ens = stats.linregress(years,ensemble_shear)
percent_shear_ens = 100*(((slope_ens*years+intercept_ens)[-1]-(slope_ens*years+intercept_ens)[0])/(slope_ens*years+intercept_ens)[0])

ensemble_speed = (speed_annual_jra+speed_annual_ncep+speed_annual_era)/3.

slope_ens_speed,intercept_ens_speed,rvalue_ens_speed,pvalue_ens_speed,stderr_ens_speed = stats.linregress(years,ensemble_speed)

# ========================== Figure =====================================
plt.rcParams.update({'font.size': 7})

fig,(ax1,ax2) = plt.subplots(2,1,figsize=(4.7,7),sharex=True)

ax1.plot(years,shear_annual_era,label="ERA-Interim",color='blue',alpha=0.7)
ax1.plot(years,shear_annual_ncep,label="NCEP/NCAR",color='orange',alpha=0.7)
ax1.plot(years,shear_annual_jra,label="JRA-55",color='green',alpha=0.7)
ax1.plot(years,ensemble_shear,label="mean",marker='o',linewidth=2,color='k')
ax1.plot(years,slope_ens*years+intercept_ens,linestyle='--',color='k',linewidth=2,label="mean trend")
ax1.set_xlim(1978,2018)
ax1.set_ylabel("vertical zonal shear (m s$^{-1}$ (100 hPa)$^{-1}$)")
ax1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15),frameon=False, ncol=5)
ax1.text(-0.1,1.02, 'a', weight="bold",horizontalalignment='left',\
         verticalalignment='center', transform=ax1.transAxes)

ax2.plot(years,speed_annual_era,label="ERA-Interim",color='blue',alpha=0.7)
ax2.plot(years,speed_annual_ncep,label="NCEP/NCAR",color='orange',alpha=0.7)
ax2.plot(years,speed_annual_jra,label="JRA-55",color='green',alpha=0.7)
ax2.plot(years,ensemble_speed,label="Mean",marker='o',linewidth=2,color='k')
ax2.plot(years,slope_ens_speed*years+intercept_ens_speed,linestyle='--',color='k',linewidth=2,label="Mean trend")
ax2.set_xlim(1978,2018)
ax2.set_ylabel("zonal speed (m s$^{-1}$)")
ax2.set_xlabel("year")
ax2.text(-0.1, 1.02, 'b', weight="bold",horizontalalignment='left',\
         verticalalignment='center', transform=ax2.transAxes)
fig.tight_layout()
plt.savefig("figure_3.png",bbox_inches="tight",dpi=400)
plt.savefig("figure_3.eps",bbox_inches="tight",dpi=400)
