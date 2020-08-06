import pandas as pd

def get_waterlevel(starttime, endtime):
	path = dataDir + "/waterlevel_daily.csv"
	dataFrame = pd.read_csv(path)
	dataFrame["Time"] = pd.to_datetime(dataFrame["Time"], format="%d.%m.%Y")
	toTake = (dataFrame["Time"] >= starttime) & (dataFrame["Time"] <= endtime)
	dataFrame = dataFrame[toTake]
	return(dataFrame)


def get_rates(starttime, endtime):
	path = dataDir + "/Produktionsraten_Tagesmittel.csv"
	dataFrame = pd.read_csv(path)
	dataFrame["Datum"] = pd.to_datetime(dataFrame["Datum"], format="%d/%m/%Y")
	toTake = (dataFrame["Datum"] >= starttime) & (dataFrame["Datum"] <= endtime)
	dataFrame = dataFrame[toTake]
	return(dataFrame)

def merge_sets(waterlevelFrame, ratesFrame, maxGaps):
	times
	isInterWl = waterlevelFrame["Waterlevel"].isna()
	isInterRates = ratesFrame["Gesamt"].isna()
	check_time_gaps(waterlevelFrame["Time"])
	check_time_gaps(ratesFrame["Datum"])
	check_gap_length(waterlevelFrame["Time"], maxGaps, isInterWl)
	check_gap_length(ratesFrame["Datum"], maxGaps, isInterRates)

def cut_time_range(waterlevelFrame, ratesFrame):
	waterlevelFrameReturn = waterlevelFrame
	ratesFrameReturn = ratesFrame
	notNaWL = [b for b in waterlevelFrame["Waterlevel"].isna()]
	notNaRates = [b for b in ratesFrame["Gesamt"].isna()]
	if not notNaWL[0] & notNaRates[0]:
		warnings.warn("Found Na values at the beginning. Rate and waterlevel data will be cutted to the first common date without Na")
		commonStart = np.maximum(waterLevel["Time"][notNaWL], ratesFrame["Datum"][notNaRates])
		ratesFrameReturn = ratesFrameReturn[ratesFrameReturn["Datum"] >= commonStart]
		waterlevelFrameReturn = waterlevelFrameReturn[waterlevelFrameReturn["Time"] >= commonStart]
	if not notNaWL[-1] & notNaRates[-1]:
		warnings.warn("Found Na values at the end. Rate and waterlevel data will be cutted to the last common date without Na")
		commonEnd = np.minimum(waterLevel["Time"][notNaWL], ratesFrame["Datum"][notNaRates])
		ratesFrameReturn = ratesFrameReturn[ratesFrameReturn["Datum"] <= commonEnd]
		waterlevelFrameReturn = waterlevelFrameReturn[waterlevelFrameReturn["Time"] <= commonEnd]
	return(waterlevelFrameReturn, ratesFrameReturn)

def check_time_gaps(timeseries):
	if any(timeseries.diff().dt.days > 1):
		warnings.warn("Time gap found. This should not exist in the original data set!")
	return()

def check_gap_length(timeseries, maxGaps, isna):
	counter = 1
	priorIdx = -1
	for idx, a in enumerate(isna):
		if a:
			if (priorIdx + 1) == idx:
				counter += 1
			else:
				if counter >= maxGaps:
					message = "\nFound a long gap:" + \
					"\nStart: " + str(timeseries[priorIdx - counter + 1]) + \
					"\nEnd: " + str(timeseries[priorIdx]) + \
					"\nGap will be interpolated linearly. Expect poor deconvolution results."
				warnings.warn(message)
				counter = 1
				priorIdx = idx
	return()

