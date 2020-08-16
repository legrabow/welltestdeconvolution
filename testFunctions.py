### run tests on the code

import numpy as np
import pandas as pd
from scipy.optimize import check_grad

def check_example():
    ## test 1: Uniform spacing of nodes in the time domain
    a = 20.4
    timeseries = np.arange(0, 20,dtype=float)
    nodes = np.log(timeseries[2:])
    nodes = np.insert(nodes, 0, np.log(0.5))
    rates = np.arange(100,119)
    z = np.full(len(nodes), np.log(a))

    convMat = generate_convMatrix(z, len(rates), nodes, timeseries)
    print(pd.DataFrame(convMat))

    v1Truncated = list()
    for t in xrange(2,len(timeseries)):
        res = (np.log(timeseries[t]) - np.log(timeseries[t - 1])) * a
        v1Truncated.append(res)

    print(v1Truncated)

    jacobian = generate_jacobianConvolution(nodes, z, rates, len(rates), timeseries)
    print(pd.DataFrame(jacobian))
      
    
def check_gradients_convMat():
    ## test 2: Check jacobian by evaluating each row of the error function based on how 
    ## good its corresponding gradient approximates the numerical gradient
    row = 0
    timeseries = np.arange(0, 31)
    z = np.random.normal(loc = -12, size = len(timeseries))
    rates = np.random.normal(loc = 900,scale = 100, size = len(timeseries) - 1)
    waterlevel = np.random.normal(loc = 0.6, size = len(timeseries) - 1)
    nodes = get_nodes(startNode = 0.4, timeseries = timeseries)
    wlNatIn = max(waterlevel)
    weights["rew"] = get_rateError_weight(wlNat = wlNatIn, waterlevel = waterlevel, rates = rates)
    weights["dw"] = get_derivate_weight(wlNat = wlNatIn, waterlevel = waterlevel)
    rates = rates.reshape((len(rates),1))
    waterlevel = waterlevel.reshape((len(waterlevel),1))
    z = z.reshape((len(z), 1))
    x0 = np.vstack([wlNatIn, rates, z])
    resid = check_grad(error_measure, grad_measure, x0, rates, weights, waterlevel, nodes, timeseries, row)
    print(resid)
    
    
    
def grad_measure(total, rates, weights, waterlevel, nodes, timeseries, row):
    x = total[:(len(rates) + 1)]
    z = total[(len(rates) + 1):]
    fMat = generate_fMatrix(weights["rew"], z, len(rates), len(waterlevel), nodes, timeseries)
    zJacobian = generate_jacobian(nodes, z, x[1:], weights["dw"], len(rates), timeseries)
    totalJacobian = np.hstack((fMat, zJacobian))
    return totalJacobian[row,:]
    
def error_measure(total, rates, weights, waterlevel, nodes, timeseries, row):
    x = total[:(len(rates) + 1)]
    z = total[(len(rates) + 1):]
    fMat = generate_fMatrix(weights["rew"], z, len(rates), len(waterlevel), nodes, timeseries)
    fMat = fMat[:len(waterlevel),:]
    vVec = generate_vVector(weights["rew"], weights["dw"], z, waterlevel, rates, nodes)
    vVec = vVec[:len(waterlevel),:]
    error = fMat.dot(x) - vVec
    return error[row]
## test 3: Check determinant of jacobian
