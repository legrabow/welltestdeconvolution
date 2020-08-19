### main script for deconvolution

## import build-in
from datetime import datetime
import numpy as np
import pickle
import matplotlib.pyplot as plt

## import other functions
#from functions import *
#from data_preparation import *
#from variableProjection import *
#from weightFunctions import *

## important parameters

# start time of time series to look at
starttime = datetime.strptime("01/01/2001","%d/%m/%Y")
# end time of time series to look at (standard: 31/12/2005)
endtime = datetime.strptime("01/02/2006","%d/%m/%Y")
# amount of time nodes of the response function
# increasing the amount increases its resolution but makes the whole TLS-problem "less overdetermined"
# if no amount is given by the user the amount will be equal to the amount of rates
amountNodes = None

## refining parameters

# directory where to find the rate- and water level data
dataDir = "/home/leonard/Documents/Praktikum/hydraulic_data"
# the critical ratio of decrease in error after one VPA-cycle and the original error telling VPA when to stop
stoppingCriterion = 10**-8
# time of the first node in percentage of one day (so 1 > startNode > 0)
startNode = 0.2
# optional weights to consider
# "consider" defines whether the case should be considered
# "Envelope" defines if the enveloping function of the tidal signal should be taken
# "Variance" defines the sharpness of the negative probability density function 
# "Sharpness" defines the sharpness of the arctan function for the beginning
# "Shift" defines the shift of the arctan function for a later start
weighingFunctions = {
 "tidals":{"Consider":False, "Envelope":True},
 "holidays":{"Consider":False, "Variance":1},
 "beginning":{"Consider":False, "Sharpness":0.5,"Shift":0.1},
 "pumping":{"Consider":False}
}
# maximal gaps that are acceptable without giving out a warning about bad deconvolution results
maxGaps = 10
# provide your own assumed natural water level. Very important if the time series does not exhibit recovery periods!
wlNatIn = None
# known transmissivity in m**2/d (pumping test 1979)
transmissivity = 320
# relative size to scale down the gradient vector found in Gauss-Newton
#stepsize = 0.05

### prepare time series

## get barometrically corrected data

waterlevelRaw = get_waterlevel(starttime = starttime, endtime = endtime, dataDir = dataDir)
ratesRaw = get_rates(starttime = starttime, endtime = endtime, dataDir = dataDir)

## clean and merge data

timeseries, waterlevelTot, ratesTot, isInterWl, isInterRates, timeRange = process_data(waterlevelFrame = waterlevelRaw, ratesFrame = ratesRaw, maxGaps = maxGaps)
waterlevel = waterlevelTot[1:,]
rates = ratesTot[1:,]


### set up nodes

nodes = get_nodes(amountNodes = amountNodes, interpolation="linear", startNode = startNode, timeseries = timeseries)

### calculate initial values

if not wlNatIn:
    wlNatIn = get_initial_wlNat(waterlevel = waterlevelTot)
yIn = rates
zIn = get_initial_responses(nodes = nodes, waterlevel = waterlevel, rates = yIn, wlNat = wlNatIn, timeseries = timeseries)
xIn = np.insert(yIn, 0,wlNatIn, axis=0)

### get weighting functions

weights = dict()

## own weights

weights["totalWeightMatrix"] = get_total_weight(weighingFunctions, len(nodes), rates, len(timeseries), timeRange)

## general weights

weights["rew"] = 0#get_rateError_weight(wlNat = wlNatIn, waterlevel = waterlevel, rates = rates)
weights["dw"] = get_derivate_weight(wlNat = wlNatIn, waterlevel = waterlevel)

### set up monitoring objects
fMatOut = dict()
entriesConvMat = dict()


### solve the non-linear LTS-Problem

y, z, wlNat = variable_projection(nodes, waterlevel, xIn, rates, zIn, stoppingCriterion, weights, timeseries)

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

if False:
    fileName = str(starttime) + "..." + str(endtime)
    with open(fileName, 'wb') as outputFile:
        pickle.dump(output, outputFile)

### show result

plt.plot(nodes, z)

