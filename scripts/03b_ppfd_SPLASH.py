# Install pyrealm if not already installed
# !pip install pyrealm

import numpy as np
import pandas as pd
from pyrealm import splash
from datetime import datetime
import os


# Read csv file for SW_IN_F from file
try:
    df_rad = pd.read_csv('../data/FLX_CH-Oe2_FLUXNET2015_FULLSET_2004-2023_1-3/FLX_CH-Oe2_FLUXNET2015_FULLSET_DD_2004-2023_1-3.csv')
except FileNotFoundError:
    print("CSV file not found. Please check the file path.")
    raise

# Only take the following columns from the data with TIMESTAMP as date with format YYYY-MM-DD
df_rad = df_rad[['TIMESTAMP', 'SW_IN_F']]
# Convert TIMESTAMP to datetime and set as index
df_rad['TIMESTAMP'] = pd.to_datetime(df_rad['TIMESTAMP'], format='%Y%m%d')
df_rad.set_index('TIMESTAMP', inplace=True)

# Define the location and elevation
latitude = 47.286417
longitude = 7.733750
elevation = 452  # meters

# Date range
start_date = "2004-01-01"
end_date = "2023-12-31"
dates = pd.date_range(start=start_date, end=end_date, freq='D')

# Prepare input DataFrame for SPLASH
df = pd.DataFrame({
    'date': dates,
    'year': dates.year,
    'doy': dates.dayofyear,
    'latitude': latitude,
    'longitude': longitude,
    'elevation': elevation
})

# Merge SW_IN_F (W m-2) into df by date
df = df.merge(df_rad, left_on='date', right_index=True, how='left')

# Convert SW_IN_F from W m-2 to MJ m-2 d-1: 1 W m-2 = 0.0864 MJ m-2 d-1 (for daily mean)
df['SW_IN_F_MJ'] = df['SW_IN_F'] * 0.0864

# Calculate PPFD from SW_IN_F if available, otherwise use SPLASH
# Conversion: 1 W m-2 ≈ 2.04 μmol m-2 s-1 (for PAR, assuming 45% of SW is PAR, and 4.57 μmol/J)
# For daily sum: PPFD (mol m-2 d-1) = SW_IN_F (W m-2) * 0.45 * 4.57 * 86400 / 1e6
def sw_to_ppfd(sw):
    # sw: W m-2 (daily mean)
    if pd.isna(sw):
        return np.nan
    return sw * 0.45 * 4.57 * 86400 / 1e6

df['ppfd_obs'] = df['SW_IN_F'].apply(sw_to_ppfd)

# Calculate PPFD: use observed if available, otherwise SPLASH
ppfd = []
for i, row in df.iterrows():
    if not pd.isna(row['ppfd_obs']):
        ppfd.append(row['ppfd_obs'])
    else:
        try:
            result = splash.run_one_day(lat=latitude, lon=longitude, elv=elevation, 
                                       n=row['doy'], year=row['year'], sf=1.0)
            ppfd.append(result['ppfd'])
        except Exception as e:
            print(f"Error calculating SPLASH for {row['date']}: {e}")
            ppfd.append(np.nan)

df['ppfd'] = ppfd

# Show first few rows
print(df[['date', 'SW_IN_F', 'ppfd_obs', 'ppfd']].head())

# Create output directory if it doesn't exist
output_dir = '../data/01_data_prep'
os.makedirs(output_dir, exist_ok=True)

# save the DataFrame to a CSV file
df.to_csv(os.path.join(output_dir, '04_splash_ppfd.csv'), index=False)
print(f"Data saved to {output_dir}/04_splash_ppfd.csv")