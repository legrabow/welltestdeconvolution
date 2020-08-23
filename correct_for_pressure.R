library(zoo)
setwd("/home/leonard/Documents/Praktikum/")

# barometric efficiency:
be <- 0.58

sanity_check_dates <- function(dates){
  if(any(is.na(dates))){
    stop("Conversion to datetime object went wrong!")
  }
  return()
}

sanity_check_merging <- function(){
  timesWhenuapai <- time(mergeWhenuapai)
  times1340 <- time(merge1340)
  #check that times are correct
  dateCheck <- tail(c(timesWhenuapai, times1340), -1)
  if(is.character(all.equal(dateCheck, dateWL))){
    stop("Pressure data set is not complete!")
  }
  
  #check that values are correct
  for(dt in datePressureWhenuapai){
    if(dt %in% timesWhenuapai){
      p1 <- pressureWhenuapai[datePressureWhenuapai == dt]
      p2 <- as.numeric(mergeWhenuapai[timesWhenuapai == dt])
      if(p1 != p2){
        stop("Merging went wrong!")
      }
    }
  }
  for(dt in datePressure1340){
    if(dt %in% times1340){
      p1 <- pressure1340[datePressure1340 == dt]
      p2 <- as.numeric(merge1340[times1340 == dt])
      if(p1 != p2){
        stop("Merging went wrong!")
      }
    }
  }
  return("Everything O.K.")
}

plot_averaging <- function(starttime = "01/01/2000", endtime = "01/01/2001"){
  start = strptime(starttime, "%d/%m/%Y", tz="GMT")
  end = strptime(endtime, "%d/%m/%Y", tz="GMT")
  #compares start and end have 00:00 as %H:%M by default
  toTakeHourly <-  datetimePressure >= start & datetimePressure <= end
  if(sum(toTakeHourly) < 2){
    return("Given time range is to narrow")
  }
  print(paste("Plot for station:", station))
  title <- paste(station, "\n", format(start, "%d/%m/%y-%H:%M:%S"), "until", format(end, "%d/%m/%y-%H:%M:%S"))
  plot(datetimePressure[toTakeHourly], pressure[toTakeHourly], xlab = "Time", ylab = "Pressure", type = "l", main = title, cex.main=0.8)
  halfADay <- 12*60*60
  if(station == "1340"){
    lines(as.POSIXlt.date(datePressure1340, tz="GMT") + halfADay, pressure1340, col="red")
  }else{
    lines(as.POSIXlt.date(datePressureWhenuapai, tz="GMT") + halfADay, pressureWhenuapai, col="red")
  }
  return()
}

plot_merging <- function(){
  minDate <- min(c(datePressureWhenuapai, datePressure1340))
  maxDate <- max(c(datePressureWhenuapai, datePressure1340))
  plot(datePressureWhenuapai, pressureWhenuapai, type="l", col = "green", xlim = c(minDate, maxDate))
  lines(datePressure1340, pressure1340, col="red" )
  lines(c(time(mergeWhenuapai), time(merge1340)),as.numeric(c(mergeWhenuapai, merge1340)), col = "brown")
  return()
}

plot_total_pressure <- function(starttime = "01/01/2005", endtime = "01/01/2021"){
  start = as.Date(starttime, "%d/%m/%Y", tz="GMT")
  end = as.Date(endtime, "%d/%m/%Y", tz="GMT")
  toTake <-  dateWL >= start & dateWL <= end
  if(sum(toTake) < 2){
    return("Given time range is to narrow")
  }
  title <- paste(format(start, "%d/%m/%y-%H:%M:%S"), "until", format(end, "%d/%m/%y-%H:%M:%S"))
  plot(dateWL[toTake], waterlevel[toTake], type = "l", col="red", xlab = "Time", ylab = "Waterlevel")
  shift <- mean(waterlevel[toTake], na.rm = T) - mean(totalPressure[toTake], na.rm = T)
  scaling <- (max(waterlevel[toTake], na.rm = T) - shift) / max(totalPressure[toTake], na.rm = T)
  lines(dateWL[toTake], (totalPressure[toTake] * scaling + shift), col="green")
  abline(h=shift, col="green")
  return()
}

plot_outcome <- function(starttime = "01/01/2005", endtime = "01/01/2021"){
  start = as.Date(starttime, "%d/%m/%Y", tz="GMT")
  end = as.Date(endtime, "%d/%m/%Y", tz="GMT")
  toTake <-  dateWL >= start & dateWL <= end
  if(sum(toTake) < 2){
    return("Given time range is to narrow")
  }
  title <- paste(format(start, "%d/%m/%y-%H:%M:%S"), "until", format(end, "%d/%m/%y-%H:%M:%S"))
  plot(dateWL[toTake], waterlevel[toTake], type = "l", col="red", xlab = "Time", ylab = "Waterlevel")
  lines(dateWL[toTake], waterlevelCorrected[toTake], col="green")
  return()
}

#
# prepare waterlevel data .. waterlevel in m
#
startDate <- as.Date("30/11/1993","%d/%m/%Y")
dataWL <- read.csv("./hydraulic_data/waterlevel_daily.csv", sep = ",")
dateWL <- as.Date(dataWL[,1], "%d.%m.%Y", tz="GMT")
endDate <- tail(dateWL,1)
sanity_check_dates(dateWL)
toTake <- which(dateWL > startDate)
dateWL <- dateWL[toTake]
waterlevel <- dataWL[toTake,2]

# 
# prepare atmospheric pressure data .. pressure in m
#

# read in 1340
station <- "1340"
files <- list.files("./climate_data/Luftdruck/hourly/", station, full.names = T)
frames <- lapply(files, function(f) read.csv(f, sep = ",", dec = ".", skip = 8))
dataPressure <- Reduce(function(...) rbind(...), frames)
dataPressure <- dataPressure[dataPressure[,1] == station, ]
dataPressure <- dataPressure[!duplicated(as.character(dataPressure[,2])), ]
datetimePressure <- strptime(dataPressure[,2], "%Y%m%d:%H%M", tz="GMT")
sanity_check_dates(datetimePressure)
pressure <- dataPressure[,3] / (9.81 * 10)
# average to daily values
dates <- as.Date(datetimePressure)
datePressure1340 <- unique(dates)
pressure1340 <- sapply(datePressure1340, function(dt) mean(pressure[dates == dt], na.rm = T))


#read in whenuapai air field
station <- "whenuapai"
files <- list.files("./climate_data/Luftdruck/hourly/", station, full.names = T)
frames <- lapply(files, function(f) read.csv(f, sep = ",", dec = ".", skip = 9))
dataPressure <- Reduce(function(...) rbind(...), frames)
dataPressure <- dataPressure[dataPressure[,1] == "1410" | dataPressure[,1] == "23976", ]
dataPressure <- dataPressure[!duplicated(as.character(dataPressure[,2])), ]
datetimePressure <- strptime(dataPressure[,2], "%Y%m%d:%H%M", tz="GMT")
sanity_check_dates(datetimePressure)
pressure <- dataPressure[,3]  / (9.81 * 10)
# average to daily values
dates <- as.Date(datetimePressure)
datePressureWhenuapai <- unique(dates)
pressureWhenuapai <- sapply(datePressureWhenuapai, function(dt) mean(pressure[dates == dt], na.rm = T))

#
# apply barometric changes
#

# clip pressure sets and interpolate linearly between missing values
startDate1340 <- as.Date("16/03/2010", "%d/%m/%Y")
mergeWhenuapai <-  na.approx(zoo(pressureWhenuapai, order.by = datePressureWhenuapai), xout = seq(startDate, startDate1340 - 1, "day"))
merge1340 <- na.approx(zoo(pressure1340, order.by = datePressure1340), xout = seq(startDate1340, endDate, "day"))
sanity_check_merging()

# calculate pressure changes
diffWhenuapai <- diff(as.numeric(mergeWhenuapai))
diff1340 <- diff(as.numeric(merge1340))
# calculate subsequent head change
totalPressure <- c(diffWhenuapai, (tail(diffWhenuapai, 1) + diff1340[1]) / 2 , diff1340) * be

#correct waterlevel
waterlevelCorrected <- waterlevel + c(0,totalPressure[1:(length(totalPressure) - 1)])

#
# save output
#

finalFrame <- data.frame(format(dateWL, "%d.%m.%Y"), waterlevelCorrected)
names(finalFrame) <- names(dataWL)
write.csv(finalFrame, file = "./hydraulic_data/waterlevel_daily_CORRECTED.csv", row.names = F, sep = ",")
