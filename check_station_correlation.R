setwd("/home/leonard/Documents/Praktikum/")

# 
# prepare atmospheric pressure data for whenuapai
#
station <- "whenuapai"
files <- list.files("./climate_data/Luftdruck/hourly/", station, full.names = T)
frames <- lapply(files, function(f) read.csv(f, sep = ",", dec = ".", skip = 9))
dataPressure <- Reduce(function(...) rbind(...), frames)
dataPressure <- dataPressure[dataPressure[,1] == "1410" | dataPressure[,1] == "23976", ]
dataPressure <- dataPressure[!duplicated(as.character(dataPressure[,2])), ]
datetimePressure_whenuapai <- strptime(dataPressure[,2], "%Y%m%d:%H%M", tz="GMT")
pressure_whenuapai <- dataPressure[,3]

dates <- as.Date(datetimePressure_whenuapai)
dates_whenuapai <- unique(dates)
daily_whenuapai <- sapply(dates_whenuapai, function(dt) mean(pressure_whenuapai[dates == dt], na.rm = T))
if(any(is.na(datetimePressure_whenuapai))){
  stop()
}

# 
# prepare atmospheric pressure data for 1400
#
station <- "1400"
files <- list.files("./climate_data/Luftdruck/hourly/", station, full.names = T)
frames <- lapply(files, function(f) read.csv(f, sep = ",", dec = ".", skip = 8))
dataPressure <- Reduce(function(...) rbind(...), frames)
dataPressure <- dataPressure[dataPressure[,1] == station, ]
dataPressure <- dataPressure[!duplicated(as.character(dataPressure[,2])), ]
datetimePressure_1400 <- strptime(dataPressure[,2], "%Y%m%d:%H%M", tz="GMT")
pressure_1400 <- dataPressure[,3]

dates <- as.Date(datetimePressure_1400)
dates_1400 <- unique(dates)
daily_1400 <- sapply(dates_1400, function(dt) mean(pressure_1400[dates == dt], na.rm = T))
if(any(is.na(datetimePressure_1400))){
  stop()
}

# 
# prepare atmospheric pressure data for 1340
#
station <- "1340"
files <- list.files("./climate_data/Luftdruck/hourly/", station, full.names = T)
frames <- lapply(files, function(f) read.csv(f, sep = ",", dec = ".", skip = 8))
dataPressure <- Reduce(function(...) rbind(...), frames)
dataPressure <- dataPressure[dataPressure[,1] == station, ]
dataPressure <- dataPressure[!duplicated(as.character(dataPressure[,2])), ]
datetimePressure_1340 <- strptime(dataPressure[,2], "%Y%m%d:%H%M", tz="GMT")
pressure_1340 <- dataPressure[,3]

dates <- as.Date(datetimePressure_1340)
dates_1340 <- unique(dates)
daily_1340 <- sapply(dates_1340, function(dt) mean(pressure_1340[dates == dt], na.rm = T))
if(any(is.na(datetimePressure_1340))){
  stop()
}

#
# correlate 1400 and whenuapai
#
if(is.unsorted(datetimePressure_1400) | is.unsorted(datetimePressure_whenuapai)){
  stop("Sth went wrong with the dates")
}
datesCommon <- dates_1400[dates_1400 %in% dates_whenuapai]
if(!all.equal(datesCommon, dates_whenuapai[dates_whenuapai %in% dates_1400])){
  stop("Merging went wrong (very bad)")
}

common1400 <- daily_1400[dates_1400 %in% dates_whenuapai]
commonWhenuapai <- daily_whenuapai[dates_whenuapai %in% dates_1400]
plot(common1400, commonWhenuapai)


#
# correlate 1400 and 1340
#
if(is.unsorted(datetimePressure_1400) | is.unsorted(datetimePressure_1340)){
  stop("Sth went wrong with the dates")
}
datesCommon <- dates_1400[dates_1400 %in% dates_1340]
if(!all.equal(datesCommon, dates_1340[dates_1340 %in% dates_1400])){
  stop("Merging went wrong (very bad)")
}

common1400 <- daily_1400[dates_1400 %in% dates_1340]
common1340 <- daily_1340[dates_1340 %in% dates_1400]
plot(common1400, common1340)


#
# correlate whenuapai and 1340
#
datesCommon <- dates_whenuapai[dates_whenuapai %in% dates_1340]
if(!all.equal(datesCommon, dates_1340[dates_1340 %in% dates_whenuapai])){
  stop("Merging went wrong (very bad)")
}

commonWhenuapai <- daily_whenuapai[dates_whenuapai %in% dates_1340]
common1340 <- daily_1340[dates_1340 %in% dates_whenuapai]
plot(commonWhenuapai, common1340)
