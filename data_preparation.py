from datetime import datetime
import warnings
import pandas as pd
import numpy as np

def get_waterlevel(starttime, endtime, dataDir):
    ### read in water level data
    ### return a pandas data frame
    path = dataDir + "/waterlevel_daily.csv"
    dataFrame = pd.read_csv(path)
    dataFrame["Time"] = pd.to_datetime(dataFrame["Time"], format="%d.%m.%Y")
    toTake = (dataFrame["Time"] >= starttime) & (dataFrame["Time"] <= endtime)
    dataFrame = dataFrame[toTake]
    return dataFrame


def get_rates(starttime, endtime, dataDir):
    ### read in rates data
    ### return a pandas data frame
    path = dataDir + "/Produktionsraten_Tagesmittel.csv"
    dataFrame = pd.read_csv(path)
    dataFrame["Datum"] = pd.to_datetime(dataFrame["Datum"], format="%d/%m/%Y")
    toTake = (dataFrame["Datum"] >= starttime) & (dataFrame["Datum"] <= endtime)
    dataFrame = dataFrame[toTake]
    return dataFrame

def process_data(waterlevelFrame, ratesFrame, maxGaps):
    ### do some data checking and cleaning
    ### return the water level and rates series as numpy objects, the RELATIVE time series in  
    ### unit of days and a boolean vector telling if entries have been interpolated
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
    if any(np.array(ratesFrame["Datum"]) != np.array(waterlevelFrame["Time"])):
        raise Exception("Something went wrong during the data treatment!")
    ## convert data to expected data types
    timeseries = waterlevelFrame["Time"] - waterlevelFrame["Time"].iloc[0]
    timeseries = np.array(timeseries) / pd.Timedelta(1, unit='d')
    timeRange = [waterlevelFrame["Time"].iloc[0], waterlevelFrame["Time"].iloc[-1]]
    waterlevel = np.array(waterlevelFrame["Waterlevel"], dtype = "float", ndmin = 2).T
    rates = np.array(ratesFrame["Gesamt"], dtype = "float", ndmin = 2).T
    return timeseries, waterlevel, rates, isInterWl, isInterRates, timeRange

def cut_time_range(waterlevelFrame, ratesFrame):
    ### cut both frames to the same start- and end time and therby neglect Na's at head and tail 
    ### return the cutted pandas frames
    print("Cut the rate- and the water level series to the same end- and start time if necessary")
    ## check if both have the same start time. If not, cut them
    if waterlevelFrame["Time"].iloc[0] != ratesFrame["Datum"].iloc[0]:
        print("\nThe start time of the rate- and the water level sets are not the same")
        commonStart = np.maximum(waterlevelFrame["Time"].iloc[0], ratesFrame["Datum"].iloc[0])
        ratesFrame = ratesFrame[ratesFrame["Datum"] >= commonStart]
        waterlevelFrame = waterlevelFrame[waterlevelFrame["Time"] >= commonStart]
    ## check if both have the same end time. If not, cut them
    if waterlevelFrame["Time"].iloc[-1] != ratesFrame["Datum"].iloc[-1]:
        print("\nThe end time of the rate- and the water level sets are not the same")
        commonEnd = np.minimum(waterlevelFrame["Time"].iloc[-1], ratesFrame["Datum"].iloc[-1])
        ratesFrame = ratesFrame[ratesFrame["Datum"] <= commonEnd]
        waterlevelFrame = waterlevelFrame[waterlevelFrame["Time"] <= commonEnd]
    ## get boolean which values are not Na's
    notNaWL = [not b for b in waterlevelFrame["Waterlevel"].isna()]
    notNaRates = [not b for b in ratesFrame["Gesamt"].isna()]
    ## check if both start without Na's. If not, cut them to the first non-Na start
    if not notNaWL[0] & notNaRates[0]:
        print("\nFound Na values at the beginning")
        commonStart = np.maximum(waterlevelFrame["Time"][notNaWL].iloc[0], ratesFrame["Datum"][notNaRates].iloc[0])
        ratesFrame = ratesFrame[ratesFrame["Datum"] >= commonStart]
        waterlevelFrame = waterlevelFrame[waterlevelFrame["Time"] >= commonStart]
    ## check if both end without Na's. If not, cut them to the latest non-Na end
    if not notNaWL[-1] & notNaRates[-1]:
        print("\nFound Na values at the end")
        commonEnd = np.minimum(waterlevelFrame["Time"][notNaWL].iloc[-1], ratesFrame["Datum"][notNaRates].iloc[-1])
        ratesFrame = ratesFrame[ratesFrame["Datum"] <= commonEnd]
        waterlevelFrame = waterlevelFrame[waterlevelFrame["Time"] <= commonEnd]
    return waterlevelFrame, ratesFrame

def check_time_gaps(timeseries):
    ### check if there exist any gaps in the time series and issue a warning if so
    ### void function
    print("Check for time gaps")
    if any(timeseries.diff().dt.days > 1):
        raise Exception("\nTime gap found. This should not exist in the original data set!")

def check_gap_length(timeseries, maxGaps, isna):
    ### check if there exist Na-gaps longer than maxGaps (provided by the user) in both data sets
    ### and issue a warning if so
    ### void function
    counter = 1
    priorIdx = -2
    for idx, a in enumerate(isna):
        if a:
            if (priorIdx + 1) == idx:
                counter += 1
            else:
                if counter >= maxGaps:
                    message = "\nFound a long gap:" \
                    "\nStart: " + str(timeseries.iloc[priorIdx - counter + 1]) + \
                    "\nEnd: " + str(timeseries.iloc[priorIdx]) + \
                    "\nGap will be interpolated linearly. Expect poor deconvolution results."
                    warnings.warn(message)
                counter = 1
            priorIdx = idx
    if counter >= maxGaps:
        message = "\nFound a long gap:" \
        "\nStart: " + str(timeseries.iloc[priorIdx - counter + 1]) + \
        "\nEnd: " + str(timeseries.iloc[priorIdx]) + \
        "\nGap will be interpolated linearly. Expect poor deconvolution results."
        warnings.warn(message)
