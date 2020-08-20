import numpy as np

def get_tidal_weight(takeEnv, starttime, endtime):
    return

def get_holiday_weight(var, starttime, endtime):
    return

def get_beginning_weight(shs, shift, length):
    if shs <= 0:
        raise Exception("The sharpness must be positive! No inverting allowed..")
    if shift < 0:
        raise Exception("The shift must be positive!")
    
    weightValues = [np.arctan((d - shift) * shs) / np.pi + 0.5 for d in xrange(length)]
    weightMat = np.zeros((length, length))
    row,col = np.diag_indices(weightMat.shape[0])
    weightMat[row,col] = weightValues
    return weightMat

def get_pumping_weight(rates):
    return

def get_derivative_contrains(uncTime, relWeight):
    return

def get_total_weight(weighingFunctions, nodes, rates, timeLength, timeRange):
    ## define shapes
    nrowWl = timeLength - 1
    nrowRate = len(rates)
    nrowSmt = len(nodes) - 1
    
    ## get weights for water level errors
    wlError = np.identity(nrowWl)
    fillMat = np.zeros((nrowRate + nrowSmt, nrowWl))
    if weighingFunctions["tidals"]["Consider"]:
        tidMat = get_tidal_weight(weighingFunctions["tidals"]["Envelope"], timeRange[0], timeRange[1])
        wlError = tidMat.dot(wlError)
    if weighingFunctions["holidays"]["Consider"]:
        holMat = get_holiday_weight(weighingFunctions["holidays"]["Variance"], timeRange[0], timeRange[1])
        wlError = holMat.dot(wlError)
    if weighingFunctions["beginning"]["Consider"]:
        begMat = get_beginning_weight(weighingFunctions["beginning"]["Sharpness"], weighingFunctions["beginning"]["Shift"], nrowWl)
        wlError = begMat.dot(wlError)
    wlTotal = np.vstack([wlError, fillMat])
        
    ## get weights for rate errors
    rateError = np.identity(nrowRate)
    fillMat1 = np.zeros((nrowWl, nrowRate))
    fillMat2 = np.zeros((nrowSmt, nrowRate))
    if weighingFunctions["pumping"]["Consider"]:
        rateMat = get_pumping_weight(rates)
        rateError = rateMat.dot(rateError)
    rateTotal = np.vstack([fillMat1, rateError, fillMat2])
        
    ## get weights for smoothness measure
    smtError = np.identity(nrowSmt)
    fillMat = np.zeros((nrowWl + nrowRate, nrowSmt))
    if weighingFunctions["derivativeContrain"]["Consider"]:
        uncTime = weighingFunctions["derivativeContrain"]["unconstrainedTime"]
        relWeight = weighingFunctions["derivativeContrain"]["relWeight"]
        derMat = get_derivative_constrains(uncTime, relWeight, nodes)
        smtError = derMat.dot(smtError)
    smtTotal = np.vstack([fillMat, smtError])
    
    ## combine all
    wMat = np.hstack([wlTotal, rateTotal, smtTotal])
    return wMat

def get_rateError_weight(wlNat, waterlevel, rates):
    ### calculate weight based on average square drawdown and rate
    ### no GCV altertnative in use as errors are allocated for sure
    drawdown = wlNat - waterlevel
    rew = len(rates) * np.linalg.norm(drawdown) ** 2 / (len(waterlevel) * np.linalg.norm(rates) ** 2)
    return rew

def get_derivate_weight(wlNat, waterlevel):
    ### calculate weight based on average square drawdown
    drawdown = wlNat - waterlevel
    dw = np.linalg.norm(drawdown) ** 2 / len(waterlevel)
    return dw

