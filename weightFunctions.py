import numpy as np

def get_tidal_weight(takeEnv, starttime, endtime):
    return

def get_holiday_weight(var, starttime, endtime):
    return

def get_beginning_weight(shs, shift, length):
    [for d in xrange(length) np.arctan(d)]
    return

def get_pumping_weight(yIn):
    return

def get_total_weight(weighingFunctions, nodeLength, rates, timeLength, timeRange):
    ## get weights for water level errors
    nrow = timeLength - 1
    wlError = np.identity(nrow)
    if weighingFunctions["tidals"]["Consider"]:
        tidMat = get_tidal_weight(weighingFunctions["tidals"]["Envelope"], timeRange[0], timeRange[1])
        wlError = tidMat.dot(wlError)
    if weighingFunctions["holidays"]["Consider"]:
        holMat = get_holiday_weight(weighingFunctions["holidays"]["Variance"], timeRange[0], timeRange[1])
        wlError = holMat.dot(wlError)
    if weighingFunctions["beginning"]["Consider"]:
        begMat = get_beginning_weight(weighingFunctions["beginning"]["Sharpness"], weighingFunctions["beginning"]["Shift"], nrow)
        wlError = begMat.dot(wlError)
        
    ## get weights for rate errors
    nrow = len(rates)
    rateError = np.identity(nrow)
    if weighingFunctions["pumping"]["Consider"]:
        rateMat = get_pumping_weight(rates)
        rateError = rateMat.dot(rateError)
        
    ## get weights for smoothness measure
    nrow = len(nodes) - 1
    smtError = np.identity(nrow)
    
    ## combine all
    wMat = np.vstack([wlError, rateError, smtError])
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
