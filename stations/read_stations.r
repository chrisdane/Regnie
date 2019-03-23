## R

read_regnie_stations <- function(path="stations/", 
                                 f="RR_Tageswerte_Beschreibung_Stationen.txt") {

    # ftp://ftp-cdc.dwd.de/pub/CDC/observations_germany/climate/daily/more_precip/recent/RR_Tageswerte_Beschreibung_Stationen.txt
    if (F) {
        path <- "stations/"
        f <- "RR_Tageswerte_Beschreibung_Stationen.txt"
    }
    message("read ", path, f, " ...")
    con <- file(paste0(path, f), encoding="latin1") # ä ö ü
    widths <- c(5, 9, 9, 15, 12, 11, 41)
    widths <- c(widths, 200-sum(widths))
    header <- strsplit(readLines(con, n=1), " ")[[1]]
    stations <- read.fwf(con, skip=2, widths=widths, col.names=header)
    #close(con)
    #stations

}
