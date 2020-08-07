import pandas as pd

def get_waterlevel(starttime, endtime):
	### read in water level data
	### return a pandas data frame
	path = dataDir + "/waterlevel_daily.csv"
	dataFrame = pd.read_csv(path)
	dataFrame["Time"] = pd.to_datetime(dataFrame["Time"], format="%d.%m.%Y")
	toTake = (dataFrame["Time"] >= starttime) & (dataFrame["Time"] <= endtime)
	dataFrame = dataFrame[toTake]
	return(dataFrame)


def get_rates(starttime, endtime):
	### read in rates data
	### return a pandas data frame
	path = dataDir + "/Produktionsraten_Tagesmittel.csv"
	dataFrame = pd.read_csv(path)
	dataFrame["Datum"] = pd.to_datetime(dataFrame["Datum"], format="%d/%m/%Y")
	toTake = (dataFrame["Datum"] >= starttime) & (dataFrame["Datum"] <= endtime)
	dataFrame = dataFrame[toTake]
	return(dataFrame)

def process_data(waterlevelFrame, ratesFrame, maxGaps):
	### do some data checking and cleaning
	### return the water level and rates series as numpy objects, the time series and 
	### a boolean vector telling if entries have been interpolated
	## check if gaps exist in the time series of both frames (should not)
 	check_time_gaps(waterlevelFrame["Time"])
	check_time_gaps(ratesFrame["Datum"])
	## cut both data frames to the same length in case one starts/ends with Na values
	waterlevelFrame, ratesFrame = cut_time_range(waterlevelFrame, ratesFrame)
	## get the boolean list which values are Na for future weighing in the VP-algorithm
	isInterWl = waterlevelFrame["Waterlevel"].isna()
	isInterRates = ratesFrame["Gesamt"].isna()
	## check if there exist very long data gaps that might effect the deconvolution outcome
	check_gap_length(waterlevelFrame["Time"], maxGaps, isInterWl)
	check_gap_length(ratesFrame["Datum"], maxGaps, isInterRates)
	## Interpolate all gaps linearly
	waterlevelFrame = waterlevelFrame.interpolate()
	ratesFrame = ratesFrame.interpolate()
	## check if everything went fine
	if any(ratesFrame["Datum"] != waterlevelFrame["Time"]):
		raise Exception("Something went wrong during the data treatment!")
	## convert data to expected data types
	timeseries = list(waterlevel["Time"])
	waterlevel = np.array(waterlevelFrame["Waterlevel"])
	rates = np.array(ratesFrame["Gesamt"])
	return(timeseries, waterlevel, rates, isInterWl, isInterRates)

def cut_time_range(waterlevelFrame, ratesFrame):
	### in case Na's occour at the start or end of waterlevelFrame or ratesFrame, remove these and cut both
	### frames to the same start- and/or end time
	### return the cutted pandas frames
	waterlevelFrameReturn = waterlevelFrame
	ratesFrameReturn = ratesFrame
	notNaWL = [b for b in waterlevelFrame["Waterlevel"].isna()]
	notNaRates = [b for b in ratesFrame["Gesamt"].isna()]
	if not notNaWL[0] & notNaRates[0]:
		warnings.warn("\nFound Na values at the beginning. Rate and waterlevel data will be cutted to the first common date without Na")
		commonStart = np.maximum(waterLevel["Time"][notNaWL], ratesFrame["Datum"][notNaRates])
		ratesFrameReturn = ratesFrameReturn[ratesFrameReturn["Datum"] >= commonStart]
		waterlevelFrameReturn = waterlevelFrameReturn[waterlevelFrameReturn["Time"] >= commonStart]
	if not notNaWL[-1] & notNaRates[-1]:
		warnings.warn("\nFound Na values at the end. Rate and waterlevel data will be cutted to the last common date without Na")
		commonEnd = np.minimum(waterLevel["Time"][notNaWL], ratesFrame["Datum"][notNaRates])
		ratesFrameReturn = ratesFrameReturn[ratesFrameReturn["Datum"] <= commonEnd]
		waterlevelFrameReturn = waterlevelFrameReturn[waterlevelFrameReturn["Time"] <= commonEnd]
	return(waterlevelFrameReturn, ratesFrameReturn)

def check_time_gaps(timeseries):
	### check if there exist any gaps in the time series and issue a warning if so
	### void function
	if any(timeseries.diff().dt.days > 1):
		warnings.warn("\nTime gap found. This should not exist in the original data set!")
	return()

def check_gap_length(timeseries, maxGaps, isna):
	### check if there exist Na-gaps longer than maxGaps (provided by the user) in both data sets
	### and issue a warning if so
	### void function
	counter = 1
	priorIdx = -1
	for idx, a in enumerate(isna):
		if a:
			if (priorIdx + 1) == idx:
				counter += 1
			else:
				if counter >= maxGaps:
					message = "\nFound a long gap:" \
					"\nStart: " + str(timeseries[priorIdx - counter + 1]) \
					"\nEnd: " + str(timeseries[priorIdx]) \
					"\nGap will be interpolated linearly. Expect poor deconvolution results."
				warnings.warn(message)
				counter = 1
				priorIdx = idx
	return()

