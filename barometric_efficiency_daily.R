setwd("/home/leonard/Documents/Praktikum/")

calculate_be <- function(dates){
  print("calculate the barometric efficiency...")
  datesCommon <- dates[dates %in% as.Date(datetimePressure) & dates %in% as.Date(datetimeWL)]
  
  if(length(datesCommon) < 2){
    print("Cannot find enough common data points")
    return()
  }
  
  subWaterlevel <- sapply(datesCommon, function(dt) mean(waterlevel[as.Date(datetimeWL) == dt]))
  subPressure <- sapply(datesCommon, function(dt) mean(pressure[as.Date(datetimePressure) == dt]))
  
  if(any(is.na(subWaterlevel)) || any(is.na(subPressure))){
    stop()
  }
  
  subWaterlevelDiffs <- -diff(subWaterlevel)
  subPressureDiffs <- diff(subPressure)
  fit <- lm(subWaterlevelDiffs ~ subPressureDiffs)
  n <- fit$coefficients[1]
  m <- fit$coefficients[2]
  corCoef <- cor(subWaterlevelDiffs, subPressureDiffs)
  amountDatapoints <- length(subWaterlevelDiffs)
  toReturn <- c(n, m , corCoef, amountDatapoints)
  
  png(paste(datesCommon[1], "png", sep = "."))
  titleName <- paste(datesCommon[1], "bis", tail(datesCommon, 1))
  titleName <- paste(titleName, corCoef, sep = " ,CorCoef:")
  titleName <- paste(titleName, m, sep = " ,Slope:")
  
  plot(subPressureDiffs, subWaterlevelDiffs, main=titleName, ylab = "Änderung Wasserstand [cm]", xlab = "Änderung atmosphärischer Druck [cm]")
  abline(fit)
  grid()
  dev.off()
  
  return(toReturn)
}


#
# prepare precipitation data 
#
data1400 <- read.csv("./climate_data/niederschlag_daily/niwa/1400.txt", skip = 8, sep = ",", dec = ".",)
dataWarkworth <- read.csv("./climate_data/niederschlag_daily/niwa/warkworth.txt", skip = 9, sep = ",", dec = ".")

date1400 <- as.Date(data1400[,2], "%Y%m%d:%H%M")
dateWarkworth <- as.Date(dataWarkworth[,2], "%Y%m%d:%H%M")
if(any(is.na(date1400)) || any(is.na(dateWarkworth))){
  stop()
}
prec1400 <- data1400[,3]
precWarkworth <- dataWarkworth[,3]

#cut both to the same times
if(is.unsorted(date1400) | is.unsorted(dateWarkworth)){
  stop("Dates are not sorted!")
}
dateUnified <- date1400[date1400 %in% dateWarkworth]
prec1400Unified <- prec1400[date1400 %in% dateWarkworth]
precWarkworthUnified <- precWarkworth[dateWarkworth %in% date1400]

#
# prepare waterlevel data .. waterlevel in cm
#
dataWL <- read.csv("./hydraulic_data/Wasserstand-tief_Hourly.csv", sep = ";", dec = ",", skip = 3)
waterlevelFACTOR <- dataWL[2:nrow(dataWL),5]
waterlevel <- as.numeric(gsub(",", ".", as.character(waterlevelFACTOR)))
datetimeWL <- strptime(dataWL[2:nrow(dataWL),1], "%d.%m.%Y %H:%M", tz="GMT")
if(any(is.na(datetimeWL))){
  stop()
}
toRemove <- is.na(waterlevel)
waterlevel <- waterlevel[!toRemove] * 100
datetimeWL <- datetimeWL[!toRemove]

# 
# prepare atmospheric pressure data .. pressure in cm
#
station <- "1340"
#station <- "whenuapai"
files <- list.files("./climate_data/Luftdruck/hourly/", station, full.names = T)
frames <- lapply(files, function(f) read.csv(f, sep = ",", dec = ".", skip = 8))
#frames <- lapply(files, function(f) read.csv(f, sep = ",", dec = ".", skip = 9))
dataPressure <- Reduce(function(...) rbind(...), frames)
dataPressure <- dataPressure[dataPressure[,1] == station, ]
#dataPressure <- dataPressure[dataPressure[,1] == "1410" | dataPressure[,1] == "23976", ]
dataPressure <- dataPressure[!duplicated(as.character(dataPressure[,2])), ]
datetimePressure <- strptime(dataPressure[,2], "%Y%m%d:%H%M", tz="GMT")
pressure <- dataPressure[,3] * 10 / 9.81
if(any(is.na(datetimePressure))){
  stop()
}

#
# calculate barometric efficiency
#
threshold <- 0.5
barometricEfficiency <- data.frame(c(NA, NA, NA, NA))
noRainDates <- dateUnified[prec1400Unified <= threshold & precWarkworthUnified <= threshold]
diffs <- diff(noRainDates)
consistentDays <- c(noRainDates[1])
for (idx in seq_along(diffs)){
  difference <- diffs[idx]
  if(difference == 1){
    consistentDays <- c(consistentDays, noRainDates[idx+1])
  }else{
    print("Reached end of one time period without rain: ")
    print(consistentDays)
          be <- calculate_be(consistentDays)
    if(length(be) != 0){
      timeRange <- paste(consistentDays[1], tail(consistentDays, 1), sep = "...")
      subBarometricEfficiency <- data.frame(be,row.names = c("n", "m", "corrCoef", "AmountDatapoints"))
      names(subBarometricEfficiency) <- c(timeRange)
      barometricEfficiency <- cbind(barometricEfficiency, subBarometricEfficiency)
    }
    consistentDays <- c(noRainDates[idx+1])
  }
}
barometricEfficiency <- barometricEfficiency[,2:ncol(barometricEfficiency)]
#
#plotting
#
dates <- strsplit(names(barometricEfficiency), "\\.\\.\\.")
dates <- sapply(dates, function(dt) (dt[1]))
dates <- as.Date(dates, "%Y-%m-%d")
be <- barometricEfficiency[2,]
corrCoef <- barometricEfficiency[3,]
library(ggplot2)
d = data.frame(x=dates,y=as.numeric(be),Korrelationskoeffizient=as.numeric(corrCoef))
out <- ggplot(data = d, aes(x = x, y = y)) +  geom_point(aes(colour = Korrelationskoeffizient))
out+ scale_colour_gradientn(colours=topo.colors(5),
                            breaks=c(-1,0,1),labels=c("-1","0","1"),
                            limits=c(-1,1)) + labs(title= "Mindestlänge eines Zeitraums: 12 Tage",
                                                   y="Barometrische Effizienz", x = "Zeit")

ggsave("lala",device="png",height = 4, width = 6, limitsize = F)
