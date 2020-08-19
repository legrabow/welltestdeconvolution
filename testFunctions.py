### run tests on the code

## import build-in
import numpy as np
import pandas as pd
from scipy.optimize import check_grad, approx_fprime

## import other functions
from functions import *
from data_preparation import *
from variableProjection import *
from weightFunctions import *


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
    weights = dict()
    timeseries = np.arange(0, 31)
    z = np.random.normal(loc = -12, size = len(timeseries))
    rates = np.random.normal(loc = 900,scale = 100, size = len(timeseries) - 1)
    waterlevel = np.random.normal(loc = 0.6, size = len(timeseries) - 1)
    nodes = get_nodes(startNode = 0.4, timeseries = timeseries, amountNodes=None, interpolation="linear")
    wlNatIn = max(waterlevel)
    weights["rew"] = get_rateError_weight(wlNat = wlNatIn, waterlevel = waterlevel, rates = rates)
    weights["dw"] = get_derivate_weight(wlNat = wlNatIn, waterlevel = waterlevel)
    yIn = rates + np.random.normal(loc = 0,scale = 10, size = len(timeseries) - 1)
    total = np.concatenate((yIn, z))
    total = np.insert(total, 0 , wlNatIn)
    total = list(total)
    rates = rates.reshape((len(rates),1))
    waterlevel = waterlevel.reshape((len(waterlevel), 1))
    for row in xrange(0, 62):
        resid = check_grad(error_measure, grad_measure, total, rates, weights, waterlevel, nodes, timeseries, row)
        eps = np.sqrt(np.finfo(float).eps)
        jacPython = approx_fprime(total, error_measure, eps, rates, weights, waterlevel, nodes, timeseries, row)
        jacMe = grad_measure(total, rates, weights, waterlevel, nodes, timeseries, row)
        print("################ row: " + str(row))
        print("Residuum: " + str(resid))
        print("Numeric gradient:")
        print(jacPython)
        print("Program's gradient:")
        print(jacMe)
    
    
    
def grad_measure(total, rates, weights, waterlevel, nodes, timeseries, row):
    x = total[:(len(rates) + 1)]
    x = np.array(x)
    x = x.reshape((len(x), 1))
    z = total[(len(rates) + 1):]
    z = np.array(z)
    z = z.reshape((len(z), 1))
    fMat = generate_fMatrix(weights["rew"], z, len(rates), len(waterlevel), nodes, timeseries)
    zJacobian = generate_jacobian(nodes, z, x[1:], weights["dw"], len(rates), timeseries)
    totalJacobian = np.hstack((fMat, zJacobian))
    return list(totalJacobian[row,:])
    
def error_measure(total, rates, weights, waterlevel, nodes, timeseries, row):
    x = total[:(len(rates) + 1)]
    x = np.array(x)
    x = x.reshape((len(x), 1))
    z = total[(len(rates) + 1):]
    z = np.array(z)
    z = z.reshape((len(z), 1))
    fMat = generate_fMatrix(weights["rew"], z, len(rates), len(waterlevel), nodes, timeseries)
    vVec = generate_vVector(weights["rew"], weights["dw"], z, waterlevel, rates, nodes)
    error = fMat.dot(x) - vVec
    return error[row][0]


def check_gradients_convMat_NORATES():
    ## test 3: Same as test2, but check for the case of no rate-errors
    weights = dict()
    timeseries = np.arange(0, 31)
    z = np.random.normal(loc = -12, size = len(timeseries))
    rates = np.random.normal(loc = 900,scale = 100, size = len(timeseries) - 1)
    waterlevel = np.random.normal(loc = 0.6, size = len(timeseries) - 1)
    nodes = get_nodes(startNode = 0.4, timeseries = timeseries, amountNodes=None, interpolation="linear")
    wlNatIn = max(waterlevel)
    weights["rew"] = 0
    weights["dw"] = get_derivate_weight(wlNat = wlNatIn, waterlevel = waterlevel)
    total = np.insert(z, 0 , wlNatIn)
    total = list(total)
    rates = rates.reshape((len(rates),1))
    waterlevel = waterlevel.reshape((len(waterlevel), 1))
    for row in xrange(0, 62):
        resid = check_grad(error_measure_NORATES, grad_measure_NORATES, total, rates, weights, waterlevel, nodes, timeseries, row)
        eps = np.sqrt(np.finfo(float).eps)
        jacPython = approx_fprime(total, error_measure_NORATES, eps, rates, weights, waterlevel, nodes, timeseries, row)
        jacMe = grad_measure_NORATES(total, rates, weights, waterlevel, nodes, timeseries, row)
        print("################ row: " + str(row))
        print("Residuum: " + str(resid))
        print("Numeric gradient:")
        print(jacPython)
        print("Program's gradient:")
        print(jacMe)
    
    
    
def grad_measure_NORATES(total, rates, weights, waterlevel, nodes, timeseries, row):
    z = total[1:]
    z = np.array(z)
    z = z.reshape((len(z), 1))
    fMat = generate_fMatrix(weights["rew"], z, len(rates), len(waterlevel), nodes, timeseries)
    A = fMat[:,0]
    A = A.reshape((len(A), 1))
    zJacobian = generate_jacobian(nodes, z, rates, weights["dw"], len(rates), timeseries)
    totalJacobian = np.hstack((A, zJacobian))
    return list(totalJacobian[row,:])
    
def error_measure_NORATES(total, rates, weights, waterlevel, nodes, timeseries, row):
    x = np.vstack([total[0], rates])
    z = total[1:]
    z = np.array(z)
    z = z.reshape((len(z), 1))
    fMat = generate_fMatrix(weights["rew"], z, len(rates), len(waterlevel), nodes, timeseries)
    vVec = generate_vVector(weights["rew"], weights["dw"], z, waterlevel, rates, nodes)
    error = fMat.dot(x) - vVec
    return error[row][0]

