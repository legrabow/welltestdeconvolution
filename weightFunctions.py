def get_tidal_weighting(timeseries):
    return

def get_holiday_weighting(timeseries):
    return

def get_beginningWeight(beginningWeightPar):
    return

def get_pumping_weigthing(yIn):
    return

def get_rateError_weight(wlNat, waterlevel, rates):
    ### calculate weight based on average square drawdown and rate
    ### no GCV altertnative in use as errors are allocated for sure
    drawdown = wlNat - waterlevel
    rew = len(rates) * np.linalg.norm(drawdown) ** 2 / (len(waterlevel) * np.linalg.norm(rates) ** 2)
    return rew

def get_derivate_weight(wlNat, waterlevel, rates):
    ### calculate weight based on average square drawdown
    drawdown = wlNat - waterlevel
    dw = np.linalg.norm(drawdown) ** 2 / len(waterlevel)
    return dw
