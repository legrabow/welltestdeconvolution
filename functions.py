import warnings
import numpy as np

def get_nodes(amountNodes, interpolation, startNode, timeseries):
    ### calculate and return the nodes corresponding to the logarithmic time steps of the response function
    if not amountNodes:
        amountNodes = len(timeseries)
    n0 = np.log(startNode * timeseries[1])
    nEnd = np.log(timeseries[-1])
    nodes = list()
    if interpolation == "linear":
        ## use linear interpolation scheme, between the start - and end node
        for k in xrange(amountNodes):
            n = n0 + k * (nEnd - n0) / (amountNodes - 1)
            nodes.append(n)
    return nodes

def get_initial_wlNat(waterlevel):
    ### return the maximum waterlevel value as the natural waterlevel
    message = "\nNo inital guess for the natural water level given." \
    "The maximum water level in the time range will be taken instead." \
    "\nNote that this makes only sense if recovery periods exist in the time range!"
    warnings.warn(message)
    wlNat = max(waterlevel)
    return wlNat

def get_initial_responses(nodes, waterlevel, rates, wlNat, timeseries, wMat):
    ### calculate and return an inital response function as a start value for the
    ### iterative VP-algorithm
    ### first guess: well bore storage before and radial flow after the first node
    ## calculate the best fit for a radial flow
    wlLength = len(waterlevel)
    wlError = wMat[:wlLength,:wlLength]
    radialFlow = np.zeros(shape = (wlLength,len(rates)))
    firstIntegral = np.log(timeseries[1]) - nodes[0]
    fullIntegrals = [np.log(tEnd) - np.log(tStart) for tStart, tEnd in zip(timeseries[1:], timeseries[2:])]
    v1 = np.insert(fullIntegrals, 0, firstIntegral)
    for i in xrange(wlLength):
        radialFlow[i:, i] = v1[:(wlLength - i)]
    convolution = radialFlow.dot(rates)
    convolution = wlError.dot(convolution)
    drawdown = wlNat - waterlevel
    drawdown = wlError.dot(drawdown)
    rfCoef = np.vdot(convolution, drawdown) / np.linalg.norm(convolution) ** 2
    ## calculate well bore storage coefficient
    wbsCoef = rfCoef * np.exp(-nodes[0])
    ## create response array
    responsesRadial = np.full(shape = ((len(nodes) - 1), 1), fill_value = np.log(rfCoef))
    responsesTotal = np.insert(responsesRadial, 0, np.log(wbsCoef), axis = 0)
    return responsesTotal
