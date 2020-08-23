import warnings
import numpy as np

def get_nodes(amountNodes, interpolation, startNode, timeseries):
    ### calculate and return the nodes corresponding to the logarithmic time steps of the response function
    if not amountNodes:
        amountNodes = len(timeseries)
    if interpolation == "linear - node domain":
        ## use linear interpolation scheme, between the start - and end node
        n0 = np.log(startNode * timeseries[1])
        nEnd = np.log(timeseries[-1])
        nodes = np.linspace(n0, nEnd, amountNodes)
    if interpolation == "linear - time domain":
        ## use linear interpolation scheme, between the start - and end time
        t0 = startNode * timeseries[1]
        tEnd = timeseries[-1]
        times = np.linspace(t0, tEnd, amountNodes)
        nodes = np.log(times)
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
    v2 = np.insert(fullIntegrals, 0, firstIntegral)
    for i in xrange(wlLength):
        radialFlow[i:, i] = v2[:(wlLength - i)]
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
