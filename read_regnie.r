## R

## DWD monthly and daily precipitation sum REGNIE (in mm)
## ftp://ftp-cdc.dwd.de/pub/CDC/grids_germany/monthly/regnie/
## ftp://ftp-cdc.dwd.de/pub/CDC/grids_germany/daily/regnie/
## hydromet@dwd.de 

# Rauthe, M., H. Steiner, U. Riediger, A. Mazurkiewicz and A. Gratzki, 2013: A Central
# European precipitation climatology – Part I: Generation and validation of a high-
# resolution gridded daily data set (HYRAS), Vol. 22(3), p 235–256, DOI:10.1127/0941-
# 2948/2013/0436.

# Abteilung Hydrometeorologie: REGNIE (REGionalisierte NIEederschläge): Verfah-
# rensbeschreibung & Nutzeranleitung, interner Bericht im DWD, Offenbach 2017.

# 1. Die Anzahl der Stationen schwankt von Jahr zu Jahr und beeinflusst die Quali-
# tät der Rasterdaten, da Niederschlagsfelder häufig, vor allem im Sommer,
# räumlich sehr inhomogen sind.

# 2. Ein Knackpunkt der Herstellung der Rasterdaten ist die Berechnung der Hin-
# tergrundfelder. Von deren Qualität hängt die des Endproduktes im Wesentli-
# chen ab. Dabei ist vor allem in den Randbereichen die Bestimmung der Re-
# gressionskoeffizienten durch die nicht vorhandenen Daten außerhalb
# Deutschlands unsicher.

# filenames:
#   monhtly:
#     rasaYYMM # absolute monthly precipitaion sum of month MM of year YY (in mm) 
#     rasaYY13 # absolute annual precipitation sum of year YY (in mm)
#     rasrYY13 # relative annual precipiation sum relative to 1961-1990 peroiod (in mm)
#   daily:
#     raYYMMDD # absolute daily precipitation sum of day DD of month MM of year YY (in mm)  
# attention:
#   first file of 20th century: 
#     monthly: rasa3101 --> January 1931
#     daily:   ra310101 --> January 1st, 1931
#   first file of 21th century: 
#     monthly: rasa0001 --> January 2000
#     daily:   ra000101 --> January 1st, 2000
#   !!! this file naming convention needs to be changed from January 2031 on !!! 

# before running this script:
# 1: cd tar
# 2: if monthly: wget ftp://ftp-cdc.dwd.de/pub/CDC/grids_germany/monthly/regnie/*.tar .
# 2: if daily:   wget ftp://ftp-cdc.dwd.de/pub/CDC/grids_germany/daily/regnie/*.tar .
# 3: ls *.tar | xargs -i tar xf {}
# 4: mv *.gz ../ascii/.
# 5: cd ../bin
# 6: gzip -d *.gz
# 7: if monthly: rename RASA rasa *RASA* # repair flawed filenames
# 8: if monthly: if not needed: rm *RASR* # precip relative to 1961-1990
# 9: this script

graphics.off()
rm(list=ls())

################### user input ###################

if (T) {
    timetype <- "monthly" # "monthly" or "daily"
    from <- 193101 # YYYYMM if timetype = "monthly"
    to <- 201902 # YYYYMMDD if timetype = "daily"
    from <- 193101
    to <- 193112
} else if (F) {
    timetype <- "daily" # "monthly" or "daily"
    from <- 19310101 # YYYYMM if timetype = "monthly"
    to <- 20190316 # YYYYMMDD if timetype = "daily"
}
regniepath <- getwd()
inpath <- paste0(regniepath, "/data/", timetype, "/ascii/")
outpath <- paste0(regniepath, "/post/")
outname_nc <- paste0("regnie_", timetype, "_", from, "-", to, ".nc")
outname_fst <- paste0("regnie_", timetype, "_", from, "-", to, ".fst")

################# user input end #################

# my check
if (timetype == "daily" && Sys.info()[4] == "K") stop("no.")

# check if binaries already exist
if (file.access(paste0(outpath, outname_nc), mode=0) == 0 &&
    file.access(paste0(outpath, outname_fst), mode=0) == 0) {
    message("outfiles")
    message(paste0("   ", outpath, c(outname_nc, outname_fst), collapse="\n"))
    message("already exist. nothing to do.")

} else {

    # get ascii files between 'from' and 'to'
    fs <- list.files(inpath)
    if (timetype == "monthly") {
        # remove annual sum files in if monthly data
        fs <- fs[-which(substr(fs, 7, 8) == "13")]
        fs <- gsub("rasa", "", fs) # 0001 --> YYMM
        if (!all(nchar(fs) == 4)) {
            stop(paste0("unknown file format for timetype = ", timetype))
        }
    } else if (timetype == "daily") {
        fs <- gsub("ra", "", fs) # 190310 --> YYMMDD
        if (!all(nchar(fs) == 6)) {
            stop(paste0("unknown file format for timetype = ", timetype))
        }
    }
    yy <- as.numeric(substr(fs, 1, 2))
    mm <- as.numeric(substr(fs, 3, 4))
    timevec <- rep(2000, t=length(fs))
    timevec[yy >= 31] <- 1900
    timevec <- timevec + yy
    timevec <- paste0(timevec, sprintf("%.2i", mm))
    if (timetype == "daily") {
        timevec <- paste0(timevec, sprintf("%.2i", dd))
    }
    timevec <- as.numeric(timevec)

    # sort files chronologically
    inds <- sort(timevec, index.return=T)$ix
    fs <- fs[inds]
    timevec <- timevec[inds]

    nf <- length(fs)
    if (nf == 0) stop(paste0("no data found in ", inpath))

    # make numeric time
    time <- as.integer(substr(timevec, 1, 4)) # YYYY
    months <- as.numeric(substr(timevec, 5, 6)) # MM
    time <- time + 1/12*months - 1/12
    stop("asd")
    if (timetype == "daily") {
	    days <-  as.numeric(substr(timevec, 7, 8)) # DD
        #time <- time + 
    }
    if (any(diff(time) > diff(time)[1])) {
        message("warning: dt of 'time' is not constant")
    }

    # load packages
    library(ff)
    library(stringi)
    library(ncdf4)
    library(fst)

    # Das Gitter besteht aus 611 Rasterpunkten in West/Ost-Richtung und 971 Raster-
    # punkten in Nord/Süd-Richtung.
    # Die Auflösung beträgt 60 geogr. Sekunden längenparallel und 30 geogr. Sekunden
    # breitenparallel, d.h. 1 Längengrad setzt sich aus 60 Gitterpunkten und 1 Breitengrad
    # aus 120 Gitterpunkten zusammen.
    # Die geographischen Koordinaten der Rasterpunkte lassen sich wie folgt berechnen:
    # Schrittweite längenparallel (xdelta) : 60 geographische Sekunden
    #   xdelta = 1.0 * (1./60.) (Grad)
    # Schrittweite breitenparallel (ydelta) : 30 geographische Sekunden
    #   ydelta = 1.0 * (1./120.) (Grad)
    # Die geographischen Koordinaten (x,y) eines beliebigen Rasterpunktes R(m,n) :
    #   x(m,n) = ( 6. - 10.*xdelta) + (n - 1)*xdelta
    #   y(m,n) = (55. + 10.*ydelta) - (m - 1)*ydelta,
    # wobei x in West/Ost-Richtung (n von 1 bis 611) und y in Nord/Süd-Richtung (m von 1
    # bis 971) verläuft.
    nlon <- 611
    nlat <- 971
    xdelta <- 1.0 * (1/60) # deg
    ydelta <- 1.0 * (1/120) # deg
    xmat <- array(NA, c(nlat, nlon))
    ymat <- xmat
    for (m in 1:nlat) {
        for (n in 1:nlon) {
            xmat[m,n] = ( 6 - 10*xdelta) + (n - 1)*xdelta 
            ymat[m,n] = (55 + 10*ydelta) - (m - 1)*ydelta
        }
    }
    lon <- xmat[1,]
    lat <- ymat[,1]
    
    # check if lon and lat are increasing
    lon_rev <- lat_rev <- F
    if (any(diff(lon) < 0)) {
        lon <- rev(lon)
        lon_rev <- T
    }
    if (any(diff(lat) < 0)) {
        lat <- rev(lat)
        lat_rev <- T
    }

    # In den Daten sind die Rasterfelder zeilenweise von Nord nach Süd und West nach
    # Ost abgespeichert so dass sie sich wie folgt mit FORTRAN auslesen:
    #   INTEGER NDH(971,611)
    #   DO I=1,971
    #     WRITE (2,'(611i4)') (NDH(I,J),J=1,611)
    #   ENDDO
    # Der Bezugspunkt der Daten ist immer der Rastermittelpunkt.

    # Nichtbesetzte Rasterpunkte sind mit Fehlwerten -999 besetzt. 
    
    precision <- "integer"
    ninteger <- 4 # --> so for this data format, the maximum monthly precip = 9999 mm?!??!?!
    ascii_mv <- -999
    data_all <- ff::ff(0, vmode="integer", dim=c(nlon, nlat, nf))
    #data_all <- array(NA, c(nlon, nlat, nf)) # i like dimension order (lon,lat) more than (lat,lon)
    timevec <- rep(NA, t=nf)

    for (fi in 1:nf) {
        
        message(fi, "/", nf, ": ", fs[fi], " --> ", timevec[fi])

        # data
        ascii_lines <- readLines(paste0(inpath, fs[fi]))
        mat <- array(NA, c(nlon, nlat)) 
        for (i in 1:nlat) {
            # split long string every 'ninteger'
            tmp <- as.numeric(stringi::stri_sub(str=ascii_lines[i], 
                                                from=seq(1, nchar(ascii_lines[i]), b=ninteger), 
                                                l=ninteger))
            tmp[tmp == ascii_mv] <- NA
            mat[,i] <- tmp
        }

        # Die Dimension der monatlichen und jährlichen Niederschlagshöhen 
        # beträgt mm, die der täglichen mm/10.
        if (timetype == "daily") mat <- mat*10

        # check if a dimension needs to get flipped
        if (lon_rev) mat <- mat[nlon:1,]
        if (lat_rev) mat <- mat[,nlat:1]
        if (F) {
            oce::imagep(lon, lat, mat, decimate=F, useRaster=T)
            #fields::image.plot(lon, lat, mat, useRaster=T)
        }
        data_all[,,fi] <- mat
    } # for fi nf

    stop("asd")

    # make numeric time
    if (timetype == "monthly") {
        time <- as.integer(substr(timevec, 1, 4)) # YYYY
        months <- as.numeric(substr(timevec, 5, 6)) # MM
        time <- time + 1/12*months - 1/12
    } else if (timetype == "daily") {
        stop("asd")
    }
    if (any(diff(time) > diff(time)[1])) {
        message("warning: dt of 'time' is not constant")
    }

    # Messzeitraum: 07:30 GZ bis 07:30 GZ (gesetzliche Zeit) Folgetag, d.h. jeweils
    # MEZ bzw. MSEZ im Winter bzw. Sommer

    # Ab dem 01.01.1971 ist die Niederschlagshöhe auf den Tag bezogen, an dem
    # der Niederschlag gefallen ist (ein Tag bevor dieser gemessen wurde). Vor
    # 1971 ist der Tag Bezugstag, an dem gemessen wurde. Um eine kontinuierli-
    # che Zeitreihe zu bekommen, müssen die Raster vor 1971 um einen Tag zu-
    # rückdatiert werden.
    if (timetype == "daily") {
        stop("asd")
    }

    # Es wird routinemäßig täglich ein vorläufiges Produkt erstellt (dieses hat den
    # Zeitbezug 5:50 UTC bis 5:50 UTC), das alle automatisch meldenden Statio-
    # nen enthält. Die übrigen Daten werden nach Eingang und Qualitätskontrolle
    # nachträglich eingepflegt und damit die endgültigen Rasterdaten mit dem Zeit-
    # bezug (07:30-07:30 GZ) erstellt.

    # netcdf output
    if (file.access(paste0(outpath, outname_nc), mode=0) == -1) {
        if (any(precision == c("double", "float"))) {
            nc_mv <- NA
        } else {
            nc_mv <- ascii_mv
        }
        message("save ", outpath, outname_nc, " ...")
        time_dim <- ncdim_def("time", "", time)
        lon_dim <- ncdim_def("longitude", "", lon)
        lat_dim <- ncdim_def("latitude", "", lat)
        data_var <- ncvar_def(paste0("regni_", timetype), "mm",
                              list(lon_dim, lat_dim, time_dim),
                              missval=nc_mv, prec=precision)
        timevec_var <- ncvar_def("timevec", "#", time_dim, prec="integer")
        outnc <- nc_create(paste0(outpath, outname_nc), 
                           list(data_var, timevec_var), force_v4=T)
        ncvar_put(outnc, data_var, data_all[,,])
        ncvar_put(outnc, timevec_var, timevec)
        ncatt_put(outnc, 0, "ref", 
                  paste0("ftp://ftp-cdc.dwd.de/pub/CDC/grids_germany/", timetype, "/regnie/"))
        nc_close(outnc)
    }

    # fst output
    if (file.access(paste0(outpath, outname_fst), mode=0) == -1) {
        # fst saves 2-dimensional data frames
        data_all_df <- array(data_all, dim=c(nlon*nlat, nf))
        data_all_df <- as.data.frame(data_all_df)
        write_fst(data_all_df, paste0(outpath, outname_fst))
        time_df <- data.frame(time=time, timevec=timevec)
        write_fst(time_df, paste0(outpath, outname_fst, ".time"))
        lon_df <- data.frame(lon=lon)
        write_fst(lon_df, paste0(outpath, outname_fst, ".lon"))
        lat_df <- data.frame(lat=lat)
        write_fst(lat_df, paste0(outpath, outname_fst, ".lat"))
    }

} # if output files do not exist

message("finish")
