### main script for bootstrapping

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


## bootstrapping parameters
starttimeTotal = datetime.strptime("01/12/2005","%d/%m/%Y")
endtimeTotal = datetime.strptime("30/06/2017","%d/%m/%Y")
samplesize = 14 #in days
iterationNumber = 1000



# amount of time nodes of the response function
# increasing the amount increases its resolution but makes the whole TLS-problem "less overdetermined"
# if no amount is given by the user the amount will be equal to the amount of rates
amountNodes = 60

## refining parameters

# directory where to find the rate- and water level data
dataDir = "/home/leonard/Documents/Praktikum/hydraulic_data"
# the critical ratio of decrease in error after one VPA-cycle and the original error telling VPA when to stop
stoppingCriterion = 10**-8
# time of the first node in percentage of one day (so 1 > startNode > 0)
startNode = 0.8#10**(-3) / 24
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
maxGaps = 4
# provides the assumed natural water level. If "None" the maximum water level value of that period will be taken
wlNatIn = 2.8


totalDays = endtimeTotal - starttimeTotal
availableDays = totalDays.days - samplesize

zValuesMat = np.empty((iterationNumber,amountNodes,))
zValuesMat[:] = np.nan
wlNatValuesVec = np.empty((iterationNumber,))
wlNatValuesVec[:] = np.nan
errorVec = np.empty((iterationNumber,))
errorVec[:] = np.nan

starttimesDict = dict()
endtimesDict = dict()


for itnum in xrange(iterationNumber):
    print("#"*50+"\nReached iteration: "+str(itnum)+"\n"+"#"*50)
    startDay = np.random.randint(0, availableDays + 1)

    # start time of time series to look at
    starttime = starttimeTotal + timedelta(days=startDay)
    # end time of time series to look at
    endtime = starttime + timedelta(days=samplesize)
    

    ### prepare time series

    ## get barometrically corrected data
    waterlevelRaw = get_waterlevel(starttime = starttime, endtime = endtime, dataDir = dataDir)
    ratesRaw = get_rates(starttime = starttime, endtime = endtime, dataDir = dataDir)

    ## clean and merge data
    try:
        timeseries, waterlevelTot, ratesTot, isInterWl, isInterRates, timeRange = process_data(waterlevelFrame = waterlevelRaw, ratesFrame = ratesRaw, maxGaps = maxGaps)
        if not timeseries[-1] == samplesize:
            print("Had to skip...(timeseries)")
            continue
    except: 
        print("Had to skip...(gaps)")
        continue
        
    waterlevel = waterlevelTot[:-1,]
    rates = ratesTot[:-1,]


    ### set up nodes
    nodes = get_nodes(amountNodes = amountNodes, interpolation="linear - node domain", startNode = startNode, timeseries = timeseries)

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
    weights["dw"] = get_derivate_weight(wlNat = wlNatIn, waterlevel = waterlevel) * 10 ** (-7)

    ### get initial response values
    zIn = get_initial_responses(nodes, waterlevel, yIn, wlNatIn, timeseries, weights["totalWeightMatrix"])


    ### set up monitoring objects
    fMatOut = dict()
    entriesConvMat = dict()


    ### solve the non-linear LTS-Problem
    try:
        y, z, wlNat, finalError = variable_projection(nodes, waterlevel, xIn, rates, zIn, stoppingCriterion, weights, timeseries)
    except ValueError:
        print("Had to skip...(vpa)")
        continue
    zValuesMat[itnum,:] = np.reshape(z, (1,len(z)))
    wlNatValuesVec[itnum] = wlNat
    errorVec[itnum] = finalError
    starttimesDict[str(itnum)] = starttime
    endtimesDict[str(itnum)] = endtime

    
zValuesMean = np.nanmean(zValuesMat, axis=0)
zValuesStd = np.nanstd(zValuesMat, axis=0)

wlNatValuesMean = np.nanmean(wlNatValuesVec)
wlNatValuesStd = np.nanstd(wlNatValuesVec)
