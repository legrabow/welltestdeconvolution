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

startNode = 0.8#10**(-3) / 24

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
mu, std = norm.fit(rates)

zError = dict()
wlError = dict()


scaleStepSize = 0.05
realisations = 10  # per counter

for scale in np.arange(scaleStepSize, 1.5, scaleStepSize):
    print("#"*50+"\nReached scale: "+str(scale)+"\n"+"#"*50)
    
    diffWlMat = np.empty((realisations**2,lengthDiffWlMat,))
    diffWlMat[:] = np.nan
    idxMat = 0
    
    diffzValMat = np.empty((realisations**2,amountNodes,))
    diffzValMat[:] = np.nan
    
    noiseScaleWlNat = noiseScaleWL = noiseScaleRates = scale
    
    for counter1 in xrange(1,realisations+1):
        for counter2 in xrange(1,realisations+1):
            #
            ### randomise inital conditions (rates, natural water level)
            #

            np.random.seed(counter1)
            rates = np.random.normal(mu, std,(len(rates), 1)) # maybe gumbel distibution better
            theorWLNat = float(np.round(np.random.uniform(0.5,5,(1,1)), 2))
            theorZ = np.full((amountNodes, 1), -7.758)

            #
            ### forward calculation
            #

            np.random.seed(counter2)

            ## generate water level

            theorNodes = get_nodes(amountNodes = amountNodes, interpolation="linear - node domain", startNode = startNode, timeseries = timeseries)
            convMat = generate_convMatrix(theorZ, len(rates), theorNodes, timeseries)
            theorDD = convMat.dot(rates)
            theorWL = theorWLNat - theorDD

            ### add noise to actual sets

            wlNatIn = theorWLNat + float(np.round(np.random.normal(0,theorWLNat * noiseScaleWlNat), 2))
            waterlevel = theorWL.copy() + np.random.normal(0,0.2 * noiseScaleWL,(len(rates), 1))
            nodes = theorNodes.copy()
            rates = rates + np.random.normal(0,200 * noiseScaleRates,(len(rates), 1))

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

            y, z, wlNat, finalError = variable_projection(nodes, waterlevel, xIn, rates, zIn, stoppingCriterion, weights, timeseries)

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
            
    totalDiffWL = diffWlMat.mean(axis=0)
    totalDiffzVal = diffzValMat.mean(axis=0)
    
    wlError[str(scale)] = totalDiffWL
    zError[str(scale)] = totalDiffzVal
