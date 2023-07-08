# r

graphics.off()
rm(list=ls())

################### user input ###################

if (F) {
    regnie_temporal_resolution <- "monthly" # "monthly" or "daily"
    from <- 193101 # YYYYMM if regnie_temporal_resolution = "monthly"
    to <- 201902 # YYYYMMDD if regnie_temporal_resolution = "daily"
    from <- 193101
    to <- 193112
} else if (T) {
    regnie_temporal_resolution <- "daily" # "monthly" or "daily"
    # first data: 19301231
    # last data: 20190316 
    from <- 19301231 # YYYYMM if regnie_temporal_resolution = "monthly" 
    to <- 20190316 # YYYYMMDD if regnie_temporal_resolution = "daily"
    #from <- 19310101
    #to <- 19310131
    #to <- 19310101
    #from <- 19701201
    #to <- 19710131
}
regniepath <- getwd()
inpath <- paste0(regniepath, "/data/", regnie_temporal_resolution, "/gz/")
outpath <- paste0(regniepath, "/post/")
save_as <- "annual_mean" # "annual_mean" "array" or "mat"
outname <- paste0("regnie_", regnie_temporal_resolution, "_", from, "-", to, "_", save_as)
outname_nc <- paste0(outname, ".nc")
fst_compress <- 100
outname_fst <- paste0(outname, "_compress_", fst_compress, ".fst")
outname_rdata <- paste0(outname, ".rdata")

################# user input end #################

# my check
if (regnie_temporal_resolution == "daily" && Sys.info()[4] == "K") stop("no.")

# check if binaries already exist
if (file.access(paste0(outpath, outname_nc), mode=0) == 0 &&
    file.access(paste0(outpath, outname_fst), mode=0) == 0) {
    message("outfiles")
    message(paste0("   ", outpath, c(outname_nc, outname_fst), collapse="\n"))
    message("already exist. nothing to do.")

} else {

    # get ascii files between 'from' and 'to'
    message("read directory ", inpath, " ... ", appendLF=F)
    fs <- list.files(inpath)
    message(length(fs), " files found.")

    if (length(fs) == 0) {
        message("could not find any data in ", inpath)
        message("note that you need to fulfill the following steps, before running this script:")
        message("   1: $ cd /where/to/save/the/regnie/tar/files")
        if (regnie_temporal_resolution == "monthly") {
            message("   2: $ wget ftp://ftp-cdc.dwd.de/pub/CDC/grids_germany/monthly/regnie/*.tar .")
        } else if (regnie_temporal_resolution == "daily") {
            message("   2: $ wget ftp://ftp-cdc.dwd.de/pub/CDC/grids_germany/daily/regnie/*.tar .")
        }
        message("   3: $ ls *.tar | xargs -i tar xf {}")
        message("   --> .gz files needed for this script are created (you can delete the *.tar files if you want)")
        if (regnie_temporal_resolution == "monthly") {
            message("   4: $ rename RASA rasa *RASA* # force all low case filenames")
            message("   5 (optional): $ rm *RASR* # precip relative to 1961-1990 data are not needed for this script")
        } else if (regnie_temporal_resolution == "daily") {
            message("   4: take care for special date RA701231.ZIP (if this data from Dec 31 1970 is needed):")
            message("      $ mv RA701231.ZIP /where/all/other/*.gz/files/are")
            message("      --> you maybe will have two files of the same date, ra701231.gz and RA701231.ZIP. thats ok.") 
        }
        stop("abort")
    }

    message("prepare ...")
    if (regnie_temporal_resolution == "monthly") {
        
        # remove annual sum files in if monthly data
        fs <- fs[-which(substr(fs, 7, 8) == "13")]
        
        # rasa0001.gz --> YYMM
        if (!all(nchar(fs) == 11)) {
            inds <- which(nchar(fs) != 11)
            stop(paste0("files ", paste0(fs[inds], collapse=","), " are of unknown file format"))
        }
        yy <- as.numeric(substr(fs, 5, 6))
        mm <- as.numeric(substr(fs, 7, 8))
    
    } else if (regnie_temporal_resolution == "daily") {
        
        # ra190310.gz --> YYMMDD
        if (!all(any(nchar(fs) == 11)) ||
            !all(any(nchar(fs) == 12))) {
            inds <- which(nchar(fs) != 11 | nchar(fs) != 12)
            stop(paste0("files ", paste0(fs[inds], collapse=","), " are of unknown file format"))
        }
        yy <- as.numeric(substr(fs, 3, 4))
        mm <- as.numeric(substr(fs, 5, 6))
        dd <- as.numeric(substr(fs, 7, 8))
    }

    # YYYYMMDD format
    time <- rep(2000, t=length(fs))
    time[yy >= 31] <- 1900
    time <- time + yy
    time <- paste0(time, sprintf("%.2i", mm))
    if (regnie_temporal_resolution == "daily") {
        time <- paste0(time, sprintf("%.2i", dd))
    }

    # sort files chronologically 
    # (use this YYYYMMDD format since not possible to get the sort index from date-class object below?)
    inds <- sort(as.numeric(time), index.return=T)$ix
    time <- time[inds]
    fs <- fs[inds]

    # convert to date class
    tz <- "CET" # time zone necessary for converting from and to data and numeric classes
    origin <- "1970-01-01" # default date necessary for later nc output
    oldtime <- time
    time <- as.POSIXct(as.character(time), format="%Y%m%d", tz=tz)
    
    # for testing 
    if (F) {
        inds <- which(time >= as.POSIXct("19701201", format="%Y%m%d", tz=tz) & 
                      time <= as.POSIXct("19710131", format="%Y%m%d", tz=tz))
        time <- time[inds]
        fs <- fs[inds]
    }

    # Messzeitraum: 07:30 GZ bis 07:30 GZ (gesetzliche Zeit) Folgetag, d.h. jeweils
    # MEZ bzw. MSEZ im Winter bzw. Sommer

    # Ab dem 01.01.1971 ist die Niederschlagshöhe auf den Tag bezogen, an dem
    # der Niederschlag gefallen ist (ein Tag bevor dieser gemessen wurde). Vor
    # 1971 ist der Tag Bezugstag, an dem gemessen wurde. Um eine kontinuierli-
    # che Zeitreihe zu bekommen, müssen die Raster vor 1971 um einen Tag zu-
    # rückdatiert werden.

    # after my mail to DWD:
    # Der 31.12.1970 (Beobachtungszeitraum 31.12.1970 - 01.01.1971) liegt als
    # Extra-Datei vor : RA701231
    if (regnie_temporal_resolution == "daily") {
        inds_before710101 <- time < as.POSIXct("19710101", format="%Y%m%d", tz=tz)
        if (any(inds_before710101)) {
            message("date all times before 1971-01-01 back one day ...")
            if (any(fs[inds_before710101] == "RA701231.ZIP")) {
                message("take care of RA701231.ZIP ...")
                inds_before710101[which(fs == "RA701231.ZIP")] <- F
            }
            time[inds_before710101] <- time[inds_before710101] - 86400 # minus one day
        } # if data before January 1 1971
    } # if regnie_temporal_resolution == "daily"

    # Es wird routinemäßig täglich ein vorläufiges Produkt erstellt (dieses hat den
    # Zeitbezug 5:50 UTC bis 5:50 UTC), das alle automatisch meldenden Statio-
    # nen enthält. Die übrigen Daten werden nach Eingang und Qualitätskontrolle
    # nachträglich eingepflegt und damit die endgültigen Rasterdaten mit dem Zeit-
    # bezug (07:30-07:30 GZ) erstellt.

    # select time 'from' 'to'
    fromind <- which(time == as.POSIXct(as.character(from), format="%Y%m%d", tz=tz))
    toind <- which(time == as.POSIXct(as.character(to), format="%Y%m%d", tz=tz))
    if (length(fromind) == 0 || length(toind) == 0) {
        stop("could not find any files between ", from, "-", to)
    }
    fs <- fs[fromind:toind]
    time <- time[fromind:toind]
    nf <- length(fs)
    
    # all data is too big:
    if (save_as == "daily" && nlon*nlat*nf > .Machine$integer.max) {
        stop(paste0("nlon*nlat*nf = ", nlon*nlat*nf, " > .Machine$integer.max (",
                    .Machine$integer.max, "). change 'save_as' to e.g. 'monthly'"))
    }

    # make numeric time
    if (F) {
        time <- as.integer(substr(timevec, 1, 4)) # YYYY
        months <- as.numeric(substr(timevec, 5, 6)) # MM
        time <- time + months/12 - 1/12
        if (regnie_temporal_resolution == "daily") {
            days <-  as.numeric(substr(timevec, 7, 8)) # DD
            # check if data takes leap years into account
            inds_feb <- months == 2
            inds_day29 <- days == 29
            inds_feb29 <- which(inds_feb & inds_day29)
            if (length(inds_feb29) > 0) {
                source("~/scripts/r/functions/leap_function.r")
                if (any(!is.leap(floor(time[inds_feb29])))) {
                    message("warn: these februaries contain day 29 although its not a leap year:")
                    print(data.frame(fs=fs[inds_feb29], 
                                     timevec=timevec[inds_feb29], 
                                     is.leap=is.leap(floor(time[inds_feb29]))))
                }
            } # if any februrary 29 exist
            ndays_per_month <- rep(NA, t=nf)
            inds <- which(months == 1 | months == 3 | months == 5 |
                          months == 7 | months == 8 | months == 10 | months == 12)
            ndays_per_month[inds] <- 31
            inds <- which(months == 2 | months == 4 | months == 6 |
                          months == 9 | months == 11)
            ndays_per_month[inds] <- 30
            time <- time + (days - 1)/ndays_per_month*1/12
        } # if daily
    } # old

    # load packages
    message("load packges ...")
    if (T) {
        library(readr)
        library(ncdf4)
        library(fst)
    }
    if (F) library(ff)

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
    message("get lon,lat coordinates ...")
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
    
    if (save_as == "annual_mean") {
        #years_all <- floor(time)
        years_all <- format(time, "%Y") # date class --> character class
        years_all <- as.numeric(years_all) # needed later
        years <- unique(years_all)  
        nyears <- length(years)
        if (nlon*nlat*nyears > .Machine$integer.max) {
            message("nlon*nlat*nyears = ", nlon*nlat*nyears, " > .Machine$integer.max (", .Machine$integer.max, ")")
            stop("o_O this should not happen")
        }
        data_all <- array(NA, c(nlon, nlat, nyears))
        ndata_all <- array(0, c(nlon, nlat, nyears))
        ndays_mm_01 <- ndata_all
        yeari <- 0 # counter
    } else {
        if (save_as == "array") {
            #data_all <- array(NA, c(nlon, nlat, nf)) # i like dimension order (lon,lat) more than (lat,lon)
            data_all <- ff::ff(0, vmode="integer", dim=c(nlon, nlat, nf))
        } else if (save_as == "mat") {
            #data_all <- matrix(NA, nrow=nlon*nlat, ncol=nf)
            data_all <- ff::ffdf(ff::ff(0, vmode="integer", dim=c(nlon*nlat, nf)))
        }
    }

    # for all monthly of daily files
    for (fi in 1:nf) {

        message(fi, "/", nf, ": ", fs[fi], " --> ", time[fi], " --> year ", years_all[fi])

        # read data
        # in each archive is an ascii file with ~971 = nlat rows; sometimes one more due to a meta line at the end)
        # note: for whatever reason, readr::read_fwf is faster than LaF::laf_open_fwf
        # read_fwf can read the .gz / .zip archives! 
        # stolen from https://github.com/EDiLD/esmisc
        mat <- read_fwf(paste0(inpath, fs[fi]),  
                        col_positions=fwf_widths(rep(ninteger, t=nlon)), 
                        col_types=cols(.default=col_integer()),
                        n_max=nlat, na=as.character(ascii_mv))
        mat <- t(as.matrix(mat))
        
        # check if a dimension needs to get flipped
        if (lon_rev) mat <- mat[nlon:1,]
        if (lat_rev) mat <- mat[,nlat:1]
      
        # Die Dimension der monatlichen und jährlichen Niederschlagshöhen 
        # beträgt mm, die der täglichen mm/10.
        # --> 150 mm/10 --> 15 mm
        # --> original data is integer --> smallest value = 1 mm/10 --> 0.1 mm 
        if (regnie_temporal_resolution == "daily") mat <- mat/10 # mm/10 --> mm

        # check 
        if (F) {
            oce::imagep(lon, lat, mat, decimate=F, useRaster=T)
            #fields::image.plot(lon, lat, mat, useRaster=T)
        }

        # save annual chunks 
        if (save_as == "annual_mean") {
            
            # tmp arrays for summing all days of a year
            if (fi == 1 || # initial
                (years_all[fi] != years_all[fi-1])) { # or this files year is different than one file before
                data_tmp <- array(0, c(nlon, nlat))
                ndata_tmp <- array(0, c(nlon, nlat))
                ndays_mm_01_tmp <- ndata_tmp
            }
            data_inds <- !is.na(mat)
            if (any(data_inds)) {
                data_tmp[data_inds] <- data_tmp[data_inds] + mat[data_inds]
                ndata_tmp[data_inds] <- ndata_tmp[data_inds] + 1
                mm_01_inds <- mat > 0 # mm 
                if (any(is.na(mm_01_inds))) {
                    mm_01_inds[is.na(mm_01_inds)] <- F
                }
                ndays_mm_01_tmp[mm_01_inds] <- ndays_mm_01_tmp[mm_01_inds] + 1
            } else {
                message("note: not a single data point ...")
            }

            # put to annual arrays
            if (fi == nf || # last
                (years_all[fi+1] != years_all[fi])) { # or next files year is different than current
                yeari <- yeari + 1
                message("put to year ", yeari, "/", nyears, ": ", years[yeari])
                data_all[,,yeari] <- data_tmp
                ndata_all[,,yeari] <- ndata_tmp
                ndays_mm_01[,,yeari] <- ndays_mm_01_tmp
            }

        } else {
            if (save_as == "array") {
                data_all[,,fi] <- mat
            } else if (save_as == "mat") {
                data_all[,fi] <- as.vector(mat)
            }
        }
    } # for fi nf

    # netcdf output
    if (any(save_as == c("annual_mean", "array"))) {

        # if output does not exist already
        if (file.access(paste0(outpath, outname_nc), mode=0) == -1) {
            
            if (save_as == "array") {
                estimated_filesize <- nlon*nlat*nf*4/1024^3 # Gb; rough estimate
            } else if (save_as == "annual_mean") {
                estimated_filesize <- nlon*nlat*nyears*4/1024^3
            }
            if (estimated_filesize > 10) { 
                message("outfile ", outname_nc, " would need ~", estimated_filesize, " Gb. skip.")
            } else {
                if (any(precision == c("double", "float"))) {
                    nc_mv <- NA
                } else {
                    nc_mv <- ascii_mv
                }
                
                message("save ", outpath, outname_nc, " ...")
                lon_dim <- ncdim_def("longitude", "degrees_east", lon)
                lat_dim <- ncdim_def("latitude", "degrees_north", lat)
                if (save_as == "annual_mean") {
                    year_dim <- ncdim_def("year", "", years)
                    data_var <- ncvar_def("sum", "mm",
                                          list(lon_dim, lat_dim, year_dim),
                                          missval=nc_mv, prec=precision)
                    ndata_var <- ncvar_def("ndata", "#",
                                           list(lon_dim, lat_dim, year_dim),
                                           missval=nc_mv, prec="integer")
                    ndays_mm_01_var <- ncvar_def("ndays_mm_01", "#",
                                                 list(lon_dim, lat_dim, year_dim),
                                                 missval=nc_mv, prec="integer")
                    outnc <- nc_create(paste0(outpath, outname_nc), 
                                       list(data_var, ndata_var, ndays_mm_01_var), force_v4=T)
                } else {
                    time_dim <- ncdim_def(name="time", units=paste0("seconds since ", origin), 
                                          vals=as.numeric(time), 
                                          calendar=paste0("tz='", tz, "'"))
                    data_var <- ncvar_def(paste0("regni_", regnie_temporal_resolution), "mm",
                                          list(lon_dim, lat_dim, time_dim),
                                          missval=nc_mv, prec=precision)
                    outnc <- nc_create(paste0(outpath, outname_nc), 
                                       list(data_var, timevec_var), force_v4=T)
                }
                ncvar_put(outnc, data_var, data_all[,,])
                if (save_as == "annual_mean") {
                    ncvar_put(outnc, ndata_var, ndata_all[,,])
                    ncvar_put(outnc, ndays_mm_01_var, ndays_mm_01[,,])
                }
                ncatt_put(outnc, 0, "ref", 
                          paste0("ftp://ftp-cdc.dwd.de/pub/CDC/grids_germany/", regnie_temporal_resolution, "/regnie/"))
                nc_close(outnc)
            }
        } # if outname_nc does not exist
    } # if save_as == "annual_mean" "array"

    # fst output
    if (save_as != "annual_mean") {

        f <- "/work/ba0941/a270073/data/dwd/Regnie/post/regnie_daily_19301231-20190316_annual_mean.nc"
        nc <- nc_open(f)
        years <- ncvar_get(nc, "year")
        ndays <- ncvar_get(nc, "sum")
        ndays_df <- array(ndays, c(prod(dim(ndays)[1:2]), dim(ndays)[3]))
        ndays_df <- as.data.frame(ndays_df)
        colnames(ndays_df) <- years
        write_fst(ndays_df, "post/regnie_daily_19301231-20190316_annual_sum_compress100.fst", compress=100)
        ndays <- ncvar_get(nc, "ndays_mm_01")
        ndays_df <- array(ndays, c(prod(dim(ndays)[1:2]), dim(ndays)[3]))
        ndays_df <- as.data.frame(ndays_df)
        colnames(ndays_df) <- years
        write_fst(ndays_df, "post/regnie_daily_19301231-20190316_annual_ndays_mm_01_compress100.fst", compress=100)
        ndays <- ncvar_get(nc, "ndata")
        ndays_df <- array(ndays, c(prod(dim(ndays)[1:2]), dim(ndays)[3]))
        ndays_df <- as.data.frame(ndays_df)
        colnames(ndays_df) <- years
        write_fst(ndays_df, "post/regnie_daily_19301231-20190316_annual_ndata_compress100.fst", compress=100)
        lon <- ncvar_get(nc, "longitude")
        lat <- ncvar_get(nc, "latitude")
        lmax <- max(length(lon), length(lat), length(years))
        lon2 <- c(lon, rep(NA, t=lmax-length(lon)))
        lat2 <- c(lat, rep(NA, t=lmax-length(lat)))
        years2 <- c(years, rep(NA, t=lmax-length(years)))
        df < data.frame(lon=lon2, lat=lat2, years=years2)
        write_fst(ndays_df, "post/regnie_daily_19301231-20190316_coords_time_compress100.fst", compress=100)

        if (file.access(paste0(outpath, outname_fst), mode=0) == -1) {
            # fst saves 2-dimensional data frames
            if (save_as == "array") {
                stop("asd")
                data_all <- array(data_all, dim=c(nlon*nlat, nf))
                data_all <- as.data.frame(data_all)
            }
            message("save ", outpath, outname_fst, " ...")
            write_fst(data_all[,], paste0(outpath, outname_fst), compress=fst_compress)
            
            # save time, lon, lat in normal Rdata
            if (file.access(paste0(outpath, outname_rdata), mode=0) == -1) {
                message("save ", outpath, outname_rdata, " ...")
                save(time, lon, lat, file=paste0(outpath, outname_rdata))
            }
        }
    } # if save_as != "annual_mean"

} # if output files do not exist

message("finish")
