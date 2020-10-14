### main script for deconvolution

## import build-in
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
import numpy as np
import pickle
import matplotlib.pyplot as plt

## import other functions
from functions import *
from data_preparation import *
from variableProjection import *
from weightFunctions import *

## important parameters

zValues = dict()
nodeValues = dict()
wlNatValues = dict()
zInValues = dict()
timeRangeValues = dict()
errorRatioValues = dict()


# start time of time series to look at
starttime = datetime.strptime("01/01/2011","%d/%m/%Y")
# end time of time series to look at
endtime = datetime.strptime("10/01/2011","%d/%m/%Y")
# amount of time nodes of the response function
# increasing the amount increases its resolution but makes the whole TLS-problem "less overdetermined"
# if no amount is given by the user the amount will be equal to the amount of rates
amountNodes = 20

## refining parameters

# directory where to find the rate- and water level data
dataDir = "/home/leonard/Documents/Praktikum/hydraulic_data"
# the critical ratio of decrease in error after one VPA-cycle and the original error telling VPA when to stop
stoppingCriterion = 10**-8
# time of the first node in percentage of one day (so 1 > startNode > 0)
startNode = 10**(-3) / 24
# optional weights to consider
# "consider" defines whether the case should be considered
# "Envelope" defines if the enveloping function of the tidal signal should be taken
# "Variance" defines the sharpness of the negative probability density function 
# "Sharpness" defines the sharpness of the arctan function for the beginning
# "Shift" defines the shift of the arctan function for a later start
# "unconstrainedTime" defnes the length in days after which the 2nd derivative is constrained to zero
# "relWeight" defines the relative weight between the contrained and unconstrained set
weighingFunctions = {
        "tidals":{"Consider":False, "Envelope":True},
        "holidays":{"Consider":False, "Variance":1},
        "beginning":{"Consider":False, "Sharpness":10000,"Shift":1.5},
        "pumping":{"Consider":False},
        "derivativeContrain":{"Consider":False,"unconstrainedTime":30,"relWeight":2000}
        }
# maximal gaps that are acceptable without giving out a warning for bad deconvolution results
maxGaps = 10
# provides the assumed natural water level. If "None" the maximum water level value of that period will be taken
wlNatIn = 2.8

### prepare time series

## get barometrically corrected data
waterlevelRaw = get_waterlevel(starttime = starttime, endtime = endtime, dataDir = dataDir)
ratesRaw = get_rates(starttime = starttime, endtime = endtime, dataDir = dataDir)

## clean and merge data
timeseries, waterlevelTot, ratesTot, isInterWl, isInterRates, timeRange = process_data(waterlevelFrame = waterlevelRaw, ratesFrame = ratesRaw, maxGaps = maxGaps)

waterlevel = waterlevelTot[:-1,]
rates = ratesTot[:-1,]


### set up nodes
nodes = get_nodes(amountNodes = amountNodes, interpolation="linear - time domain", startNode = startNode, timeseries = timeseries)

### set initial values
if not wlNatIn:
    wlNatIn = get_initial_wlNat(waterlevel = waterlevelTot)
yIn = rates
xIn = np.insert(yIn, 0,wlNatIn, axis=0)

### get weighting functions
weights = dict()

## own weights
weights["totalWeightMatrix"] = get_total_weight(weighingFunctions, nodes, rates, len(timeseries), timeRange)

## general weights
weights["rew"] = 0#get_rateError_weight(wlNat = wlNatIn, waterlevel = waterlevel, rates = rates)
weights["dw"] = get_derivate_weight(wlNat = wlNatIn, waterlevel = waterlevel)

### get initial response values
zIn = get_initial_responses(nodes, waterlevel, yIn, wlNatIn, timeseries, weights["totalWeightMatrix"])


### set up monitoring objects
fMatOut = dict()
entriesConvMat = dict()


### solve the non-linear LTS-Problem
y, z, wlNat, finalError = variable_projection(nodes, waterlevel, xIn, rates, zIn, stoppingCriterion, weights, timeseries)

entry = str(i)
zValues[entry] = z.copy()
nodeValues[entry] = np.array(nodes)
wlNatValues[entry] = wlNat.copy()
zInValues[entry] = zIn.copy()
errorRatioValues[entry] = finalError
timeRangeValues[entry] = str(timeRange[0]) + ".." + str(timeRange[1])

### save result
output = {
    "starttime":starttime,
    "endtime":endtime,
    "zIn":zIn,
    "xIn":xIn,
    "stoppingCriterion":stoppingCriterion,
    "weights":weights,
    "nodes":nodes,
    "yOut":y,
    "zOut":z,
    "wlNatOut":wlNat,
    "OtherNotes":""
}

fileName = str(starttime) + "..." + str(endtime)
with open(fileName, 'wb') as outputFile:
    pickle.dump(output, outputFile)
