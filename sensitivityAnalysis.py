### sensitivity analysis

## import build-in
from scipy.stats import norm
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

starttime = datetime.strptime("01/01/2006","%d/%m/%Y")
endtime = datetime.strptime("01/01/2007","%d/%m/%Y")

amountNodes = 500

dataDir = "/home/leonard/Documents/Praktikum/hydraulic_data"

stoppingCriterion = 10**-8

startNode = 0.08#10**(-2) #0.8

weighingFunctions = {
        "tidals":{"Consider":False, "Envelope":True},
        "holidays":{"Consider":False, "Variance":1},
        "beginning":{"Consider":False, "Sharpness":10000,"Shift":1.5},
        "pumping":{"Consider":False},
        "derivativeContrain":{"Consider":False,"unconstrainedTime":30,"relWeight":2000}
        }

maxGaps = 10

waterlevelRaw = get_waterlevel(starttime = starttime, endtime = endtime, dataDir = dataDir)
ratesRaw = get_rates(starttime = starttime, endtime = endtime, dataDir = dataDir)

timeseries, waterlevelTot, ratesTot, isInterWl, isInterRates, timeRange = process_data(waterlevelFrame = waterlevelRaw, ratesFrame = ratesRaw, maxGaps = maxGaps)

waterlevel = waterlevelTot[:-1,]
lengthDiffWlMat = len(waterlevel)
del waterlevel
rates = ratesTot[:-1,]
lengthRates = len(rates)
mu, std = norm.fit(rates)
del rates

zError = dict()
wlError = dict()
crashes = dict()


scaleStepSize = 0.01
realisations = 5  # per counter

for scale in np.arange(scaleStepSize, 0.3, scaleStepSize):
    print("#"*50+"\nReached scale: "+str(scale)+"\n"+"#"*50)
    
    diffWlMat = np.empty((realisations**2,lengthDiffWlMat,))
    diffWlMat[:] = np.nan
    
    diffzValMat = np.empty((realisations**2,amountNodes,))
    diffzValMat[:] = np.nan
    
    noiseScaleWlNat = noiseScaleWL = noiseScaleRates = scale
    crashNumber = 0
    idxMat = 0
    
    
    limitedInfluence = True
    lastRelevantDay = 4
    
    for counter1 in xrange(1,realisations+1):
        for counter2 in xrange(1,realisations+1):
            #
            ### randomise inital conditions (rates, natural water level)
            #

            np.random.seed(counter1)
            theorRates = np.random.normal(mu, std,(lengthRates, 1)) # maybe gumbel distibution better
            theorWLNat = float(np.round(np.random.uniform(0.5,5,(1,1)), 2))
            theorZ = np.full((amountNodes, 1), -7.758)
            #a = 1#0.5
            #lastRelevantDay = 5
            #constantValueFrom = 8
            #x = np.arange(-(amountNodes-lastRelevantDay), lastRelevantDay)
            #x = x[::-1]
            #theorZ = -7.758 - np.exp(-x*a)
            #theorZ[constantValueFrom:] = theorZ[constantValueFrom]

            #
            ### forward calculation
            #

            np.random.seed(counter2)

            ## generate water level

            theorNodes = get_nodes(amountNodes = amountNodes, interpolation="linear - node domain", startNode = startNode, timeseries = timeseries)
            convMat = generate_convMatrix(theorZ, lengthRates, theorNodes, timeseries)
            if limitedInfluence:
                convMat = numpy.triu(convMat, -(lastRelevantDay - 1))
            theorDD = convMat.dot(theorRates)
            theorWL = theorWLNat - theorDD
            if np.isnan(np.sum(convMat)):
                raise Exception("Theoretical z-Values are too low!")

            ### add noise to actual sets
            
            wlNatInStd = theorWLNat * noiseScaleWlNat
            wlNatIn = theorWLNat + float(np.round(np.random.normal(0,wlNatInStd), 2))
            
            waterlevelStd = noiseScaleWL * np.linalg.norm(theorWL) / np.sqrt(len(theorWL))
            waterlevel = theorWL.copy() + np.random.normal(0, waterlevelStd, (lengthRates, 1))
            
            nodes = theorNodes.copy()
            
            ratesStd = noiseScaleRates * np.linalg.norm(theorRates) / np.sqrt(lengthRates)
            rates = theorRates + np.random.normal(0, ratesStd, (lengthRates, 1))

            #
            ### solving backwards
            #

            yIn = rates
            xIn = np.insert(yIn, 0,wlNatIn, axis=0)

            weights = dict()

            weights["totalWeightMatrix"] = get_total_weight(weighingFunctions, nodes, rates, len(timeseries), timeRange)

            weights["rew"] = 0#get_rateError_weight(wlNat = wlNatIn, waterlevel = waterlevel, rates = rates)
            weights["dw"] = get_derivate_weight(wlNat = wlNatIn, waterlevel = waterlevel) * 10 ** (-7)

            zIn = get_initial_responses(nodes, waterlevel, yIn, wlNatIn, timeseries, weights["totalWeightMatrix"])
            
            try:
                y, z, wlNat, finalError = variable_projection(nodes, waterlevel, xIn, rates, zIn, stoppingCriterion, weights, timeseries)
            except ValueError:
                crashNumber = crashNumber + 1
                continue
            #
            ### compare results
            #

            convMat = generate_convMatrix(z, len(rates), nodes, timeseries)
            checkDD = convMat.dot(rates)
            checkWL = wlNat - checkDD
            
            diffWL = waterlevel - checkWL
            diffWlMat[idxMat,:] = np.reshape(diffWL, (1,len(diffWL)))
            
            diffzVal = theorZ - z
            diffzValMat[idxMat,:] = np.reshape(diffzVal, (1,len(diffzVal)))

            idxMat = idxMat + 1
            
    totalDiffWL = np.nanmean(diffWlMat, axis=0)
    totalDiffzVal = np.nanmean(diffzValMat, axis=0)
    
    wlError[str(scale)] = totalDiffWL
    zError[str(scale)] = totalDiffzVal
    crashes[str(scale)] = crashNumber
