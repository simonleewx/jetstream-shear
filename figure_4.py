import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from netCDF4 import Dataset
import os
from scipy import stats
from mpl_toolkits.basemap import Basemap 

# List of years
years = np.arange(1979,2018,1)

# ============================ NCEP/NCAR ==============================================
ncep_dir = ""
os.chdir(ncep_dir)

# Array to store annual-mean values
shear_annual_ncep = np.zeros((len(years),17,32))
speed_annual_ncep = np.zeros((len(years),17,32))

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
	lonidx2 = np.where(lons==350)[0][0]

	lats_ncep = lats[latidx1:latidx2+1]
	lons_ncep = lons[lonidx1:lonidx2+1]

	# Extract u-wind
	uwnd = f.variables['uwnd'][:,levidx1:levidx2+1,latidx1:latidx2+1,lonidx1:lonidx2+1]

	# Calculate 250 hPa shear in m/s/100 hPa
	shear = 100*((uwnd[:,2,:,:] - uwnd[:,0,:,:])/100.)

	# 250 hPa wind
	speed = uwnd[:,1,:,:]

	# Take annual mean
	shear_annual = np.mean(shear,axis=0)
	speed_annual = np.mean(speed,axis=0)

	# Store 
	shear_annual_ncep[count] = shear_annual
	speed_annual_ncep[count] = speed_annual

	f.close()


# Statistical analysis
# Array to store slopes
slopes_ncep = np.zeros((np.shape(shear_annual_ncep)[1],np.shape(shear_annual_ncep)[2]))
# Array to store p values
pvals_ncep = np.zeros(np.shape(slopes_ncep))
for i in range(np.shape(shear_annual_ncep)[1]):
	for j in range(np.shape(shear_annual_ncep)[2]):
		slope,intercept,rvalue,pvalue,stderr = stats.linregress(years,shear_annual_ncep[:,i,j])
		slopes_ncep[i,j] = slope*10 
		pvals_ncep[i,j] = pvalue

# Mask non-sig regions
pvals_sig_ncep = np.ma.masked_greater(pvals_ncep,0.05)

#================================= ERA-Interim ======================================
# Set working directory
era_dir = ""
os.chdir(era_dir)

# Get all filenames
directory = os.getcwd()

# Array to store winds
hourly_winds = np.zeros((len(years),366*4,57,99))
hourly_winds[:] = np.nan 

hourly_speed = np.zeros((len(years),366*4,57,99))
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
			lats_era = f.variables['g4_lat_2'][:]
			lons_era = f.variables['g4_lon_3'][:]
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

			# Store value
			hourly_winds[np.where(years==year)[0][0],doy] = shear
			hourly_speed[np.where(years==year)[0][0],doy] = speed

			f.close()
			doy+=1

# Calculate annual mean
shear_annual_era = np.nanmean(hourly_winds,axis=1)
speed_annual_era = np.nanmean(hourly_speed,axis=1)

# Statistical analysis
# Array to store slopes
slopes_era = np.zeros((np.shape(shear_annual_era)[1],np.shape(shear_annual_era)[2]))
# Array to store p values
pvals_era = np.zeros(np.shape(slopes_era))
for i in range(np.shape(shear_annual_era)[1]):
	for j in range(np.shape(shear_annual_era)[2]):
		slope,intercept,rvalue,pvalue,stderr = stats.linregress(years,shear_annual_era[:,i,j])
		slopes_era[i,j] = slope*10
		pvals_era[i,j] = pvalue

# Mask non-sig regions
pvals_sig_era = np.ma.masked_greater(pvals_era,0.05)

#======================================= JRA =========================================================
# Set working directory
jra_dir = ""
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
shear_annual_jra = np.zeros((len(years),33,57))
speed_annual_jra = np.zeros((len(years),33,57))

# Loop over files by year & month
count = 0 # Index for year
for i in range(np.shape(file_list)[0]):
	print years[i]
	monthly_shear = np.zeros((np.shape(file_list)[1],33,57))
	monthly_speed = np.zeros((np.shape(file_list)[1],33,57))
	for j in range(np.shape(file_list)[1]):
		# Open file
		fname = file_list[i,j]
		print fname
		f = Dataset(fname,mode='r')
		
		# Grid
		lats_jra = f.variables['g0_lat_2'][:]
		lons_jra = f.variables['g0_lon_3'][:]
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
	
		# Monthly mean		
		monthly_shear[j] = np.mean(shear,axis=0)
		monthly_speed[j] = np.mean(speed,axis=0)

		f.close()

	# Calculate annual-mean
	annual_shear = np.mean(monthly_shear,axis=0)
	annual_speed = np.mean(monthly_speed,axis=0)

	# Store value
	shear_annual_jra[count] = annual_shear
	speed_annual_jra[count] = annual_speed
	count+=1

# Statistical analysis
# Array to store slopes
slopes_jra = np.zeros((np.shape(shear_annual_jra)[1],np.shape(shear_annual_jra)[2]))
# Array to store p values
pvals_jra = np.zeros(np.shape(slopes_jra))
for i in range(np.shape(shear_annual_jra)[1]):
	for j in range(np.shape(shear_annual_jra)[2]):
		slope,intercept,rvalue,pvalue,stderr = stats.linregress(years,shear_annual_jra[:,i,j])
		slopes_jra[i,j] = slope*10
		pvals_jra[i,j] = pvalue

# Mask non-sig regions
pvals_sig_jra = np.ma.masked_greater(pvals_jra,0.05)

# Climatological zonal-mean wind for 1980-2017
climo_ncep = np.mean(speed_annual_ncep,axis=0)
climo_era = np.mean(speed_annual_era,axis=0)
climo_jra = np.mean(speed_annual_jra,axis=0)

# ============================= CALCULATED FROM THERMAL WIND BALANCE ====================
# Set pressure level to analyse on
level = 250.

#===================================== NCEP ============================================
os.chdir(ncep_dir)

# Array to store annual-mean temperature
ann_mean_gradient = np.zeros((len(years),17,32))

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

	# Annual mean gradient
	annual_gradient = np.gradient(ann_mean_temp)[0]

	# Store annual mean
	ann_mean_gradient[count,:,:] = annual_gradient
	
	f.close()


# Lats and lons of Atlantic sector
lons_ncep = lons[lonidx1:lonidx2+1]
lats_ncep = lats[latidx1:latidx2+1]	

R=287.
# Coriolis parameter
coriolis = np.zeros((np.shape(ann_mean_gradient)[1],np.shape(ann_mean_gradient)[2]))
for i in range(len(lats_ncep)):
	for j in range(len(lons_ncep)):
		coriolis[i,j] = 2*7.29e-5*np.sin(np.radians(lats_ncep[i]))

ann_mean_gradient = (R/(250.*coriolis))*(ann_mean_gradient/(6400e3*2.5*np.pi/180.))

# Statistical analysis
# Array to store slopes
slopes_ncep_twb = np.zeros((np.shape(ann_mean_gradient)[1],np.shape(ann_mean_gradient)[2]))
# Array to store p values
pvals_ncep_twb = np.zeros(np.shape(slopes_ncep_twb))

for i in range(np.shape(ann_mean_gradient)[1]):
	for j in range(np.shape(ann_mean_gradient)[2]):
		slope,intercept,rvalue,pvalue,stderr = stats.linregress(years,ann_mean_gradient[:,i,j])
		slopes_ncep_twb[i,j] = slope*10*100
		pvals_ncep_twb[i,j] = pvalue

# Mask non-sig regions
pvals_sig_ncep_twb = np.ma.masked_greater(pvals_ncep_twb,0.05)

#================================================== ERA ===============================================
os.chdir(era_dir)

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
ann_mean_gradient = np.gradient(ann_mean_temps)[1]

# Calculate thermal wind shear
R=287.
# Coriolis parameter
coriolis = np.zeros((np.shape(ann_mean_gradient)[1],np.shape(ann_mean_gradient)[2]))
for i in range(len(lats_era)):
	for j in range(len(lons_era)):
		coriolis[i,j] = 2*7.29e-5*np.sin(np.radians(lats_era[i]))

ann_mean_gradient = (R/(250.*coriolis))*(ann_mean_gradient/(6400e3*0.75*np.pi/180.))

# Statistical analysis
# Array to store slopes
slopes_era_twb = np.zeros((np.shape(ann_mean_gradient)[1],np.shape(ann_mean_gradient)[2]))
# Array to store p values
pvals_era_twb = np.zeros(np.shape(slopes_era_twb))

for i in range(np.shape(ann_mean_gradient)[1]):
	for j in range(np.shape(ann_mean_gradient)[2]):
		slope,intercept,rvalue,pvalue,stderr = stats.linregress(years,ann_mean_gradient[:,i,j])
		slopes_era_twb[i,j] = slope*10*100
		pvals_era_twb[i,j] = pvalue

# Mask non-sig regions
pvals_sig_era_twb = np.ma.masked_greater(pvals_era_twb,0.05)

#================================================ JRA ==================================================
# Set working directory
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
ann_mean_gradient = np.zeros((len(years),33,57))

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
	annual_gradient = np.gradient(annual_temp)[0]
	# Store value
	ann_mean_gradient[count,:,:] = annual_gradient

# Calculate thermal wind shear
R=287.
# Coriolis parameter
coriolis = np.zeros((np.shape(ann_mean_gradient)[1],np.shape(ann_mean_gradient)[2]))
for i in range(len(lats_jra)):
	for j in range(len(lons_jra)):
		coriolis[i,j] = 2*7.29e-5*np.sin(np.radians(lats_jra[i]))

ann_mean_gradient = (R/(250.*coriolis))*(ann_mean_gradient/(6400e3*1.25*np.pi/180.))

# Statistical analysis
slopes_jra_twb = np.zeros(np.shape(coriolis))
pvals_jra_twb = np.zeros(np.shape(slopes_jra_twb))

for i in range(np.shape(ann_mean_gradient)[1]):
	for j in range(np.shape(ann_mean_gradient)[2]):
		slope,intercept,rvalue,pvalue,stderr = stats.linregress(years,ann_mean_gradient[:,i,j])
		slopes_jra_twb[i,j] = slope*10*100
		pvals_jra_twb[i,j] = pvalue

# Mask non-sig regions
pvals_sig_jra_twb = np.ma.masked_greater(pvals_jra_twb,0.05)

#=================================== FIGURE ==================================================================
# Set contouring levels
shear_levs = np.arange(-0.3,0.32,0.02)
wind_levs = np.arange(5,30,5)
plt.rcParams.update({'font.size': 7})
fig, ((ax1,ax4),(ax2,ax5),(ax3,ax6)) = plt.subplots(3,2,figsize=(7.2,9.7))

m = Basemap(projection='mill',llcrnrlat=30,urcrnrlat=70,llcrnrlon=280,urcrnrlon=350,resolution='l',ax=ax1)        
m.drawcoastlines(linewidth=1.25, color='#444444')    
m.drawmeridians(np.arange(0, 360, 20), labels=[0,0,0,1],linewidth=1,color='gray')
m.drawparallels(np.arange(0, 90, 20), labels=[1,0,0,0], linewidth=1,color='gray')    
x,y = np.meshgrid(lons_era,lats_era)
x1, y1 = m(x, y)
m.contourf(x1,y1,slopes_era,shear_levs,cmap="bwr",extend="both")
m.pcolor(x1, y1, pvals_sig_era, hatch='.',alpha=0.,zorder=10)
CR = m.contour(x1,y1,climo_era,wind_levs,colors='k')
plt.clabel(CR, inline=1,fmt = '%1.0f')
ax1.text(0, 1.04, 'a',verticalalignment='bottom', horizontalalignment='right',transform=ax1.transAxes,color='k', weight='bold')
ax1.set_title('ERA-Interim',weight="bold")

m = Basemap(projection='mill',llcrnrlat=30,urcrnrlat=70,llcrnrlon=280,urcrnrlon=350,resolution='l',ax=ax2)        
m.drawcoastlines(linewidth=1.25, color='#444444')     
m.drawmeridians(np.arange(0, 360, 20), labels=[0,0,0,1],linewidth=1,color='gray')
m.drawparallels(np.arange(0, 90, 20), labels=[1,0,0,0], linewidth=1,color='gray')    
x,y = np.meshgrid(lons_ncep,lats_ncep)
x1, y1 = m(x, y)
m.contourf(x1,y1,slopes_ncep,shear_levs,cmap="bwr",extend="both")
m.pcolor(x1, y1, pvals_sig_ncep, hatch='.',alpha=0.,zorder=10)
CR = m.contour(x1,y1,climo_ncep,wind_levs,colors='k')
plt.clabel(CR, inline=1,fmt = '%1.0f')
ax2.text(0, 1.04, 'b',verticalalignment='bottom', horizontalalignment='right',transform=ax2.transAxes,color='k', weight='bold')
ax2.set_title('NCEP/NCAR',weight="bold")

m = Basemap(projection='mill',llcrnrlat=30,urcrnrlat=70,llcrnrlon=280,urcrnrlon=350,resolution='l',ax=ax3)        
m.drawcoastlines(linewidth=1.25, color='#444444')    
m.drawmeridians(np.arange(0, 360, 20), labels=[0,0,0,1],linewidth=1,color='gray')
m.drawparallels(np.arange(0, 90, 20), labels=[1,0,0,0], linewidth=1,color='gray')     
x,y = np.meshgrid(lons_jra,lats_jra)
x1, y1 = m(x, y)
cf = m.contourf(x1,y1,slopes_jra,shear_levs,cmap="bwr",extend="both")
m.pcolor(x1, y1, pvals_sig_jra, hatch='.',alpha=0.,zorder=10)
CR = m.contour(x1,y1,climo_jra,wind_levs,colors='k')
plt.clabel(CR, inline=1,fmt = '%1.0f')
ax3.text(0, 1.04, 'c',verticalalignment='bottom', horizontalalignment='right',transform=ax3.transAxes,color='k', weight='bold')
ax3.set_title('JRA-55',weight="bold")

cb_ax = fig.add_axes([0.1,0.05,0.8,0.02])
cbar = fig.colorbar(cf, cax=cb_ax,ticks = np.arange(min(shear_levs),max(shear_levs)+0.1,0.1),orientation="horizontal")
cbar.set_label("trend (m s$^{-1}$ (100 hPa)$^{-1}$ decade$^{-1}$)")

m = Basemap(projection='mill',llcrnrlat=30,urcrnrlat=70,llcrnrlon=280,urcrnrlon=350,resolution='l',ax=ax4)        
m.drawcoastlines(linewidth=1.25, color='#444444')   
m.drawmeridians(np.arange(0, 360, 20), labels=[0,0,0,1],linewidth=1,color='gray')
m.drawparallels(np.arange(0, 90, 20), labels=[1,0,0,0], linewidth=1,color='gray')      
x,y = np.meshgrid(lons_era,lats_era)
x1, y1 = m(x, y)
m.contourf(x1,y1,slopes_era_twb,shear_levs,cmap="bwr",extend="both")
m.pcolor(x1, y1, pvals_sig_era_twb, hatch='.',alpha=0.,zorder=10)
CR = m.contour(x1,y1,climo_era,wind_levs,colors='k')
plt.clabel(CR, inline=1,fmt = '%1.0f')
ax4.text(0, 1.04, 'd',verticalalignment='bottom', horizontalalignment='right',transform=ax4.transAxes,color='k', weight='bold')
ax4.set_title('ERA-Interim',weight="bold")

m = Basemap(projection='mill',llcrnrlat=30,urcrnrlat=70,llcrnrlon=280,urcrnrlon=350,resolution='l',ax=ax5)        
m.drawcoastlines(linewidth=1.25, color='#444444')    
m.drawmeridians(np.arange(0, 360, 20), labels=[0,0,0,1],linewidth=1,color='gray')
m.drawparallels(np.arange(0, 90, 20), labels=[1,0,0,0], linewidth=1,color='gray')     
x,y = np.meshgrid(lons_ncep,lats_ncep)
x1, y1 = m(x, y)
m.contourf(x1,y1,slopes_ncep_twb,shear_levs,cmap="bwr",extend="both")
m.pcolor(x1, y1, pvals_sig_ncep_twb, hatch='.',alpha=0.,zorder=10)
CR = m.contour(x1,y1,climo_ncep,wind_levs,colors='k')
plt.clabel(CR, inline=1,fmt = '%1.0f')
ax5.text(0, 1.04, 'e',verticalalignment='bottom', horizontalalignment='right',transform=ax5.transAxes,color='k', weight='bold')
ax5.set_title('NCEP/NCAR',weight="bold")

m = Basemap(projection='mill',llcrnrlat=30,urcrnrlat=70,llcrnrlon=280,urcrnrlon=350,resolution='l',ax=ax6)        
m.drawcoastlines(linewidth=1.25, color='#444444')    
m.drawmeridians(np.arange(0, 360, 20), labels=[0,0,0,1],linewidth=1,color='gray')
m.drawparallels(np.arange(0, 90, 20), labels=[1,0,0,0], linewidth=1,color='gray')     
x,y = np.meshgrid(lons_jra,lats_jra)
x1, y1 = m(x, y)
m.contourf(x1,y1,slopes_jra_twb,shear_levs,cmap="bwr",extend="both")
m.pcolor(x1, y1, pvals_sig_jra_twb, hatch='.',alpha=0.,zorder=10)
CR = m.contour(x1,y1,climo_jra,wind_levs,colors='k')
plt.clabel(CR, inline=1,fmt = '%1.0f')
ax6.text(0, 1.04, 'f',verticalalignment='bottom', horizontalalignment='right',transform=ax6.transAxes,color='k', weight='bold')
ax6.set_title("JRA-55",weight="bold")

plt.savefig("figure_4.png",dpi=400,bbox_inches="tight")
plt.savefig("figure_4.eps",dpi=400,bbox_inches="tight")
