import mhkit
from mhkit.wave import resource, contours, graphics
import matplotlib.pyplot as plt
from mhkit.wave.io import ndbc
from scipy import stats
import pandas as pd
import numpy as np
from array import array
import seaborn as sns

# Specify the parameter as spectral wave density and the buoy number to be 42022
parameter = 'swden'
buoy_number = '42002'
ndbc_available_data= ndbc.available_data(parameter, buoy_number)
ndbc_available_data.head()

# Slice the available data to only include through year 2012
years_of_interest = ndbc_available_data[ndbc_available_data.year < 2013]
years_of_interest.head()

# Get dictionary of parameter data by year
filenames= years_of_interest['filename']
ndbc_requested_data = ndbc.request_data(parameter, filenames)

# Lastly we will convert a DateTime Index
ndbc_data={}

# Create a Datetime Index and remove NOAA date columns for each year
for year in ndbc_requested_data:
    year_data = ndbc_requested_data[year]
    ndbc_data[year] = ndbc.to_datetime_index(parameter, year_data)
    
# Display DataFrame of 46022 data from 1996
ndbc_data['1996'].head()

# Intialize empty lists to store the results from each year
Hm0_list=[]
Te_list=[]

# Iterate over each year and save the result in the initalized dictionary
for year in ndbc_data:
    year_data = ndbc_data[year]
    Hm0_list.append(resource.significant_wave_height(year_data.T))
    Te_list.append(resource.energy_period(year_data.T))
# Concatenate list of Series into a single DataFrame
Te = pd.concat(Te_list ,axis=0)
Hm0 = pd.concat(Hm0_list ,axis=0)
Hm0_Te = pd.concat([Hm0,Te],axis=1)
# Drop any NaNs created from the calculation of Hm0 or Te

Hm0_Te.dropna(inplace=True)
# Sort the DateTime index
# prints matrices of energy period and wave height data from chosen buoy 
Hm0_Te.sort_index(inplace=True)
Hm0_Te
print(Te)
print(Hm0)

## This part is unfinished - the matrix=mhkit... line has an error due to the use of the mhkit function

# Return period (years) of interest
period = 100

# Remove Hm0 Outliers
Hm0_Te_clean = Hm0_Te[Hm0_Te.Hm0 < 20]

# Get only the values from the DataFrame
Hm0 = Hm0_Te_clean.Hm0.values
Te  = Hm0_Te_clean.Te.values

df=pd.DataFrame({'Te': Te, 'Hm0': Hm0})
sns.jointplot(data=df,x='Te', y='Hm0', kind='kde')
matrix1 = df.to_numpy()

Hm0_bins=np.array([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15])
Te_bins=np.array([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15])
L= np.ones((15, 15), dtype=int)
matrix = mhkit.wave.performance.capture_length_matrix(Hm0, Te, L, 'frequency', Hm0_bins, Te_bins)
print(matrix)

## this section creates contour of height and period values that can be used to see most common values

# Delta time of sea-states
dt = (Hm0_Te_clean.index[2]-Hm0_Te_clean.index[1]).seconds

# Get the contour values
copula = contours.environmental_contours(Hm0, Te, dt, period, 'PCA', return_PCA=True)
Hm0_contour=copula['PCA_x1']
Te_contour=copula['PCA_x2']

fig,ax=plt.subplots(figsize=(8,4))
#%matplotlib inline
ax=graphics.plot_environmental_contour(Te, Hm0,
                                      Te_contour, Hm0_contour,
                                      data_label='NDBC 46022',
                                      contour_label='100 Year Contour',
                                      x_label = 'Energy Period, $Te$ [s]',
                                      y_label = 'Sig. wave height, $Hm0$ [m]',
                                      ax=ax)