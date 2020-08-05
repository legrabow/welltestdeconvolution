### main script for deconvolution

from datetime import datetime
import warnings

## important parameters
starttime = datetime.Date("01/01/2001","%d/%m/%Y")
endtime = datetime.Date("01/01/2008","%d/%m/%Y")
amountNodes = 

## refining parameters
dataDir = "/home/leonard/Documents/Praktikum/hydraulic_data"
stoppingCriterion = 
resolution = 
startNode = 
weightingFunctions = {
	"tidals":True,
	"holidays":True,
	"beginning":True,
	"pumping":True
}
beginningWeightPar = 

### prepare time series

## get barometrically corrected data

waterlevelRaw = get_waterlevel(starttime = starttime, endtime = endtime)
ratesRaw = get_rates(starttime = starttime, endtime = endtime)

## clean and merge data

timeseries, waterlevel, rates, isInterWl, isInterRates = merge_sets(waterlevelFrame = waterlevelRaw, ratesFrame = ratesRaw)

### set up nodes

timeLength = timeseries[-1] - timeseries[0]
nodes = get_nodes(amountNodes = amountNodes, interpolation="linear", startNode = startNode, timeLength = timeLength)

### calculate initial values

wlNatIn = get_initial_wlNat(waterlevel = waterlevel)
yIn = rates
zIn = get_initial_responses(nodes = nodes, waterlevel = waterlevel, rates = yIn, pNat = pNat)


### get weighting functions

get_tidal_weighting(timeseries)
get_holiday_weighting(timeseries)
get_beginningWeight(beginningWeightPar)
get_pumping_weigthing(yIn)

get_rateError_weight()
get_derivate_weight()

### solve the non-linear LTS-Problem

y, z, pNat = variable_projection(nodes = nodes, waterlevel = waterlevel, rates = qIn, pNat = pNatIn, z = zIn, sc = stoppingCriterion, weights = weights)

### show result



### optional scripts:
check_p1_contribtution(nodes, z, q, pNat, waterlevel)
check_rank()
