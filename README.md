# Regnie

- regnie is outdated:
> The data set "REGNIE" does not exist anymore. Is has been replayed by the "HYRAS-PRE-DE" data set in January 2022. Please use this new data set as a replacement. You can find it on our OpenData server here:
```
Daily data: 	https://opendata.dwd.de/climate_environment/CDC/grids_germany/daily/hyras_de/precipitation/
Monthly data: 	https://opendata.dwd.de/climate_environment/CDC/grids_germany/monthly/hyras_de/precipitation/
Multi-annual means: 	https://opendata.dwd.de/climate_environment/CDC/grids_germany/multi_annual/hyras_de/precipitation/
```

- Abteilung Hydrometeorologie: REGNIE (REGionalisierte NIEederschläge): Verfahrensbeschreibung & Nutzeranleitung, interner Bericht im DWD, Offenbach 2017.
- hydromet@dwd.de 
```
Rauthe, M., H. Steiner, U. Riediger, A. Mazurkiewicz and A. Gratzki, 2013: A Central
European precipitation climatology – Part I: Generation and validation of a high-
resolution gridded daily data set (HYRAS), Vol. 22(3), p 235–256, DOI:10.1127/0941-
2948/2013/0436.
```

1. Die Anzahl der Stationen schwankt von Jahr zu Jahr und beeinflusst die Qualität der Rasterdaten, da Niederschlagsfelder häufig, vor allem im Sommer, räumlich sehr inhomogen sind.
2. Ein Knackpunkt der Herstellung der Rasterdaten ist die Berechnung der Hintergrundfelder. Von deren Qualität hängt die des Endproduktes im Wesentlichen ab. Dabei ist vor allem in den Randbereichen die Bestimmung der Regressionskoeffizienten durch die nicht vorhandenen Daten außerhalb Deutschlands unsicher.

- filenames:
```
   monthly:
     rasaYYMM # absolute monthly precipitaion sum of month MM of year YY (in mm) 
     rasaYY13 # absolute annual precipitation sum of year YY (in mm)
     rasrYY13 # relative annual precipiation sum relative to 1961-1990 peroiod (in mm)
   daily:
     raYYMMDD # absolute daily precipitation sum of day DD of month MM of year YY (in mm)  
 attention:
   first file of 20th century: 
     monthly: rasa3101 --> January 1931
     daily:   ra310101 --> January 1st, 1931
   first file of 21th century: 
     monthly: rasa0001 --> January 2000
     daily:   ra000101 --> January 1st, 2000
   !!! this file naming convention needs to be changed from January 2031 on !!! 
```

# download

- download monthly data:
```
cd dwd/data/CDC/grids_germany/daily/regnie
wget ftp://ftp-cdc.dwd.de/pub/CDC/grids_germany/monthly/regnie/*.tar .
ls *.tar | xargs -i tar xf {}
#rm *.tar # if not needed anymore
rename RASA rasa *RASA* # repair flawed filenames
#rm *RASR* # if not needed: precip relative to 1961-1990
```
- download daily data:
```
cd dwd/data/CDC/grids_germany/daily/regnie
wget ftp://ftp-cdc.dwd.de/pub/CDC/grids_germany/daily/regnie/*.tar .
ls *.tar | xargs -i tar xf {}
#rm *.tar # if not needed anymore
```

# storage

- monthly files from 193101-201902 

| Type                         	| Size on disk 	| Rel Size 	|
|------------------------------	|--------------	|----------	|
| tar (90 files)               	| 307 Mb       	| 0.12     	|
| ascii (1146 files)           	| 2.6 Gb       	| 1        	|
| netcdf (1 file)              	| 2.4 Gb       	| 0.92     	|
| fst (1 file; compression 50) 	| 871 Mb       	| 0.34     	

- daily files from 19310101-20190316

| Type                         	| Size on disk 	| Rel Size 	|
|------------------------------	|--------------	|----------	|
| tar (90 files)               	| 4.3 Gb       	|          	|
| ascii (32218 files)          	| 72 Gb        	|          	|
| netcdf (1 file)              	|              	|          	|
| fst (1 file; compression 50) 	|              	|          	|


