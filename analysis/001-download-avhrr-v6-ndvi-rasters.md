to download version 6 of the AVHRR NDVI data save your token as a bash variable called "MY_TOKEN" using:

MY_TOKEN="yourTokenHere"

then run the code similar to the lines below to download all NDVI data from each satellite.

`wget -e robots=off -m -np -R .html,.tmp -nH --cut-dirs=3 "https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/PATH_TO_DATA_DIRECTORY" --header "Authorization: Bearer $MY_TOKEN" -P TARGET_DIRECTORY_ON_YOUR_FILE_SYSTEM`

code used to download all AVHRR data up to the end of 2013 (VIIRS data started in early 2014):

wget -e robots=off -m -np -A .nc -R .html,.tmp -nH -nd -q "https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/466/N07_AVH13C1/" --header "Authorization: Bearer $MY_TOKEN" -P H/GitHub/ndvi-stochasticity/data/avhrr-viirs-ndvi/raster-files/AVHRR-Land_v006/N07_AVH13C1
wget -e robots=off -m -np -A .nc -R .html,.tmp -nH -nd -q "https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/466/N09_AVH13C1/" --header "Authorization: Bearer $MY_TOKEN" -P H/GitHub/ndvi-stochasticity/data/avhrr-viirs-ndvi/raster-files/AVHRR-Land_v006/N09_AVH13C1
wget -e robots=off -m -np -A .nc -R .html,.tmp -nH -nd -q "https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/466/N11_AVH13C1/" --header "Authorization: Bearer $MY_TOKEN" -P H/GitHub/ndvi-stochasticity/data/avhrr-viirs-ndvi/raster-files/AVHRR-Land_v006/N11_AVH13C1
wget -e robots=off -m -np -A .nc -R .html,.tmp -nH -nd -q "https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/466/N14_AVH13C1/" --header "Authorization: Bearer $MY_TOKEN" -P H/GitHub/ndvi-stochasticity/data/avhrr-viirs-ndvi/raster-files/AVHRR-Land_v006/N14_AVH13C1
wget -e robots=off -m -np -A .nc -R .html,.tmp -nH -nd -q "https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/466/N16_AVH13C1/" --header "Authorization: Bearer $MY_TOKEN" -P H/GitHub/ndvi-stochasticity/data/avhrr-viirs-ndvi/raster-files/AVHRR-Land_v006/N16_AVH13C1
wget -e robots=off -m -np -A .nc -R .html,.tmp -nH -nd -q "https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/466/N18_AVH13C1/" --header "Authorization: Bearer $MY_TOKEN" -P H/GitHub/ndvi-stochasticity/data/avhrr-viirs-ndvi/raster-files/AVHRR-Land_v006/N18_AVH13C1
wget -e robots=off -m -np -A .nc -R .html,.tmp -nH -nd -q "https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/466/N19_AVH13C1/2009/" --header "Authorization: Bearer $MY_TOKEN" -P H/GitHub/ndvi-stochasticity/data/avhrr-viirs-ndvi/raster-files/AVHRR-Land_v006/N19_AVH13C1
wget -e robots=off -m -np -A .nc -R .html,.tmp -nH -nd -q "https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/466/N19_AVH13C1/2010/" --header "Authorization: Bearer $MY_TOKEN" -P H/GitHub/ndvi-stochasticity/data/avhrr-viirs-ndvi/raster-files/AVHRR-Land_v006/N19_AVH13C1
wget -e robots=off -m -np -A .nc -R .html,.tmp -nH -nd -q "https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/466/N19_AVH13C1/2011/" --header "Authorization: Bearer $MY_TOKEN" -P H/GitHub/ndvi-stochasticity/data/avhrr-viirs-ndvi/raster-files/AVHRR-Land_v006/N19_AVH13C1
wget -e robots=off -m -np -A .nc -R .html,.tmp -nH -nd -q "https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/466/N19_AVH13C1/2012/" --header "Authorization: Bearer $MY_TOKEN" -P H/GitHub/ndvi-stochasticity/data/avhrr-viirs-ndvi/raster-files/AVHRR-Land_v006/N19_AVH13C1
wget -e robots=off -m -np -A .nc -R .html,.tmp -nH -nd -q "https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/466/N19_AVH13C1/2013/" --header "Authorization: Bearer $MY_TOKEN" -P H/GitHub/ndvi-stochasticity/data/avhrr-viirs-ndvi/raster-files/AVHRR-Land_v006/N19_AVH13C1

*NOTE:* initially, it may seem that nothing is downloading because `wget` is downloading the index.html file from each directory and deleting it

`wget` arguments:

- `-e cmnd` or `--execute cmnd`: execute command `cmnd`
- `-o=logfile` or `--output-file=logfile`: log all messages to `logfile`
- `-q` or `-quiet` turn off the output printing
- `-m` or `--mirror` turn on options suitable for mirroring:
  - turn on recursion and time-stamping,
  - set infinite recursion depth
  - keep FTP directory listings
  - currently equivalent to `-r -N -l inf --no-remove-listing`.  
- `-np` never ascend to the parent directory when retrieving recursively
- `-A .nc` only accept nc files (avoids downloading files ending in "details=1")
- `-R .html,.tmp` reject html and tmp files
- `-nH` or `--no-host-directories` disable generation of host-prefixed directories (i.e., the higher-level folders)
- `-nd` ior `--no-directories` do not create a hierarchy of directories when retrieving recursively
- `--cut-dirs=3` ignore 3 directory components (i.e., subfolders) from the given location
- `--header=header-line` send additional header info along with the request
- `-P prefix’` or `--directory-prefix=prefix` will save all files to the `prefix` folder (default is `.`, the current directory)
