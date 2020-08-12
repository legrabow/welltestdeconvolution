### main script for deconvolution

from datetime import datetime
import warnings
import pandas as pd
import numpy as np
from scipy.linalg import svd, pinv

## important parameters

# start time of time series to look at
starttime = datetime.strptime("01/01/2001","%d/%m/%Y")
# end time of time series to look at
endtime = datetime.strptime("31/12/2005","%d/%m/%Y")
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
startNode = 0.3
# optional weights to consider
weightingFunctions = {
 "tidals":False,
 "holidays":False,
 "beginning":False,
 "pumping":False
}
# specifing the sharpness of the arctan function as a weight function for the beginning
beginningWeightPar = None
# maximal gaps that are acceptable without giving out a warning about bad deconvolution results
maxGaps = 10
# provide your own assumed natural water level. Very important if the time series does not exhibit recovery periods!
wlNatIn = None
# known transmissivity in m**2/d (pumping test 1979)
transmissivity = 320

### prepare time series

## get barometrically corrected data

waterlevelRaw = get_waterlevel(starttime = starttime, endtime = endtime)
ratesRaw = get_rates(starttime = starttime, endtime = endtime)

## clean and merge data

timeseries, waterlevelTot, ratesTot, isInterWl, isInterRates = process_data(waterlevelFrame = waterlevelRaw, ratesFrame = ratesRaw, maxGaps = maxGaps)
waterlevel = waterlevelTot[1:]
rates = ratesTot[1:]

### set up nodes

nodes = get_nodes(amountNodes = amountNodes, interpolation="linear", startNode = startNode, timeseries = timeseries)

### calculate initial values

if not wlNatIn:
    wlNatIn = get_initial_wlNat(waterlevel = waterlevelTot)
yIn = rates
zIn = get_initial_responses(nodes = nodes, waterlevel = waterlevel, rates = yIn, wlNat = wlNatIn, timeseries = timeseries)
xIn = np.insert(yIn, 0, wlNatIn)
xIn = xIn.reshape((len(xIn), 1))

### get weighting functions

weights = dict()

## own weights

if weightingFunctions["tidals"]:
    weights["tidals"] = get_tidal_weighting(timeseries)
if weightingFunctions["holidays"]:
    weights["holidays"] = get_holiday_weighting(timeseries)
if weightingFunctions["beginning"]:
    weights["beginning"] = get_beginningWeight(beginningWeightPar)
if weightingFunctions["pumping"]:
    weights["pumping"] = get_pumping_weigthing(yIn)

## general weights

weights["rew"] = get_rateError_weight(wlNat = wlNatIn, waterlevel = waterlevel, rates = rates)
weights["dw"] = get_derivate_weight(wlNat = wlNatIn, waterlevel = waterlevel)

### solve the non-linear LTS-Problem

y, z, wlNat = variable_projection(nodes, waterlevel, xIn, rates, zIn, stoppingCriterion, weights, stepsize, timeseries)

### show result



### optional scripts:
check_p1_contribtution(nodes, z, q, pNat, waterlevel)
check_rank()
