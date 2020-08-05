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

def merge_sets(waterlevelFrame, ratesFrame):
	



	isInterWl = waterlevelFrame["Time"].diff().dt.days != 1
	isInterWl = np.array(isInterWl)
	isInterWl[0] = False

	isInterWl = waterlevelFrame["Time"].diff().dt.days != 1
	isInterWl = np.array(isInterWl)
	isInterWl[0] = False

def check_time_gaps(timeseries):
	if any(timeseries.diff().dt.days > 1):
		warnings.warn("Time gap found. This should not exist in the original data set!")
	return()

def check_gap_length(timeseries, maxGaps):
	timeseries.diff().dt.days !=

def is_na(timeseries):
	isNA = index = dataFrame['Waterlevel'].index[dataFrame['Waterlevel'].apply(np.isnan)]
