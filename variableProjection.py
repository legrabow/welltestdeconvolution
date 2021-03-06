import warnings
import numpy as np
from scipy.linalg import svd, pinv, lstsq
np.set_printoptions(threshold=np.inf)

def variable_projection(nodes, waterlevel, x, rates, z, sc, weights, timeseries):
    ### VP-algorithm, each cycle calculates the new z and x (i.e. pNat and y) values
    
    errorRatio = 1
    errorRatioChange = 1
    wMat = weights["totalWeightMatrix"]
    ## calculate initial matrices
    print("Generate the fMatrix")
    fMat = generate_fMatrix(weights["rew"], z, len(rates), len(waterlevel), nodes, timeseries, wMat)
    print("Generate the vVector")
    vVec = generate_vVector(weights["rew"], weights["dw"], z, waterlevel, rates, nodes, wMat)
    
    ## calculate initial error
    errorIn = np.linalg.norm(fMat.dot(x) - vVec) ** 2
    errorBefore = errorIn
    print("Initial error: " + str(errorIn))
    
    while errorRatio > sc:# and errorRatioChange > 10**-10:
        print("Start with a solving cycle...")

        ## solve linear part
        print("Solve the linear sub problem")
        constant = vVec - fMat.dot(x)
        if weights["rew"] == 0:
            A = fMat[:,0]
            A = A.reshape((len(A), 1))
            slicePoint = 0 
            dp = np.vdot(A, constant) / (np.linalg.norm(A) ** 2)
            # set new vector
            x[0] += dp
        else:
            A = fMat
            slicePoint = len(rates)
            dx = lstsq(A, constant)[0]
            # set new vector
            x += dx
        ## solve non-linear part
        print("Solve the non-linear sub problem")
        currentPosition = -(fMat.dot(x) - vVec)
        zJacobian = generate_jacobian(nodes, z, x[1:], weights["dw"], len(rates), timeseries, wMat)

        totalJacobian = np.hstack((A, zJacobian))
        dTotal = lstsq(totalJacobian, currentPosition)[0]
        print("Total gradient vector: " + str(np.linalg.norm(dTotal)))
        
        # find stepsize
        stepsize = find_stepsize(fMat, vVec, dTotal, rates, nodes, waterlevel, timeseries, z, x, weights, totalJacobian, slicePoint, wMat)

        # set new vector
        dTotal = dTotal * stepsize
        x[:(slicePoint + 1)] += dTotal[:(slicePoint + 1)]
        z = z + dTotal[(slicePoint + 1):]
        
        ## calculate new matrices
        print("Generate the fMatrix")
        fMat = generate_fMatrix(weights["rew"], z, len(rates), len(waterlevel), nodes, timeseries, wMat)
        print("Generate the vVector")
        vVec = generate_vVector(weights["rew"], weights["dw"], z, waterlevel, rates, nodes, wMat)
        
        ## calculate error
        error = np.linalg.norm(fMat.dot(x) - vVec) ** 2
        #errorRatioChange = abs(errorRatio - error / errorIn) / errorRatio
        errorRatio = (errorBefore - error) / errorIn
        errorBefore = error
        print("Calculated error ratio: " + str(errorRatio))

    return x[1:], z, x[0], error

def find_stepsize(fMat, vVec, dTotal, rates, nodes, waterlevel, timeseries, z, x, weights, totalJacobian, slicePoint, wMat):
    ### find the optimal stepsize for a given gradient vector to reduce the error
    ### the power of 0.5 will be continuously increased until the error is smaller than the prior error
    print("Find optimal stepsize")
    minCondition = True
    reducingPower = 0
    a = fMat.dot(x) - vVec
    c = totalJacobian.dot(dTotal)
    #global fMatBefore
    #fMatBefore = fMat.copy()
    while minCondition:
        stepsize = 0.5 ** reducingPower
        if weights["rew"] == 0:
            xTest = x.copy()
            xTest[0] += dTotal[0] * stepsize
        else:
            xTest = x + dTotal[:(slicePoint + 1)] * stepsize
        zTest = z + dTotal[(slicePoint + 1):] * stepsize
        fMatTest = generate_fMatrix(weights["rew"], zTest, len(rates), len(waterlevel), nodes, timeseries, wMat)
        if not np.isnan(np.sum(fMatTest)):
            vVecTest = generate_vVector(weights["rew"], weights["dw"], zTest, waterlevel, rates, nodes, wMat)
            b = fMatTest.dot(xTest) - vVecTest
            #global fMatAfter
            #fMatAfter = fMatTest.copy()
            result = np.linalg.norm(a)**2 - np.linalg.norm(b)**2 - 0.5 * stepsize * np.linalg.norm(c)**2
            minCondition = result < 0
        reducingPower += 1
    print("Stepsize found: " + str(stepsize))
    return stepsize

def generate_vVector(rew, dw, z, waterlevel, rates, nodes, wMat):
    ### calculate the vector that contains the waterlevel, the wheighted rates
    ### and the smoothness measure for the total error function
    smoothness = generate_smoothnessMeasure(dw, nodes, z)
    vVec = np.vstack([waterlevel, rates * np.sqrt(rew), smoothness])
    vVec = wMat.dot(vVec)
    return vVec

def generate_smoothnessMeasure(dw, nodes, z):
    ### calulculate the smoothness measure to minimize
    sdMat = generate_secondDerivativeMatrix(nodes)
    expected = np.zeros((len(nodes) - 1, 1))
    expected[0] = 1
    smoothness = (sdMat.dot(z) - expected) * np.sqrt(dw)
    return smoothness

def generate_jacobian(nodes, z, y, dw, rateLength, timeseries, wMat):
    ### calculate the jacobian matrix of the error measure with respect to the 
    ### response values for Gauss-Newton
    print("Generate the jacobian matrix with respect to the response values:")
    jacobianConvolution = - generate_jacobianConvolution(nodes, z, y, rateLength, timeseries)
    sdMat = - generate_secondDerivativeMatrix(nodes) * np.sqrt(dw)
    jacobianRates = np.zeros(shape=(len(y), len(z)))
    zJacobian = np.concatenate((jacobianConvolution, jacobianRates, sdMat))
    zJacobian = wMat.dot(zJacobian)
    return zJacobian

def generate_jacobianConvolution(nodes, z, y, rateLength, timeseries):
    ### calculate the jacobian matrix of the convolution matrix with respect to the response values
    ### Workflow: For the deriavtive for z[k], find the time enclosing the corresponding nodes (no
    ###           de[k-1] until nodes[k]) and split the time by the time points that lie in that ti
    ###           me range. Evaluate the integral for each subinterval. Only do the calculation for 
    ###           one pumping period, as the rest is symmetric
    print("Generate the convolution part of jacobian")
    timeLength = len(timeseries) - 1
    jacobianConvolution = np.zeros(shape=(timeLength, len(z)))
    for k in xrange(len(z)):
        v1 = np.zeros((timeLength, 1))
        if k == 0:
            endTriangle = np.exp(nodes[k + 1])
            timeIdxStart = 0
            timeIdxEnd = np.searchsorted(timeseries,[endTriangle])[0]
            derivativeEntry = evaluate_integral(nodes[0], nodes[-1], 0, endTriangle, z[-1], z[0])
            v1[0] = derivativeEntry
            nodeCurrent = nodes[k + 1]
            nodeBefore = nodes[k]
            zCurrent = z[k + 1]
            zBefore = z[k]
            for timeIdx in xrange(timeIdxStart, timeIdxEnd):
                subStart = timeseries[timeIdx]
                subEnd = timeseries[timeIdx + 1]
                derivativeEntry = evaluate_triangle_integral(nodeCurrent, nodeBefore, subStart, subEnd, zBefore, zCurrent, True)
                v1[timeIdx] += derivativeEntry
        else:
            nodeCurrent = nodes[k]
            nodeBefore = nodes[k - 1]
            startTriangle = np.exp(nodeBefore)
            endTriangle = round(np.exp(nodeCurrent), 6)
            timeIdxStart = np.searchsorted(timeseries,[startTriangle])[0] - 1
            timeIdxEnd = np.searchsorted(timeseries,[endTriangle])[0]
            zCurrent = z[k]
            zBefore = z[k - 1]
            for timeIdx in xrange(timeIdxStart, timeIdxEnd):
                subStart = timeseries[timeIdx]
                subEnd = timeseries[timeIdx + 1]
                derivativeEntry = evaluate_triangle_integral(nodeCurrent, nodeBefore, subStart, subEnd, zBefore, zCurrent, False)
                v1[timeIdx] = derivativeEntry
            if k != (len(z) - 1):
                startTriangle = np.exp(nodes[k])
                endTriangle = round(np.exp(nodes[k + 1]), 6)
                timeIdxStart = np.searchsorted(timeseries,[startTriangle])[0] - 1
                timeIdxEnd = np.searchsorted(timeseries,[endTriangle])[0]
                nodeCurrent = nodes[k + 1]
                nodeBefore = nodes[k]
                zCurrent = z[k + 1]
                zBefore = z[k]
                for timeIdx in xrange(timeIdxStart, timeIdxEnd):
                    subStart = timeseries[timeIdx]
                    subEnd = timeseries[timeIdx + 1]
                    derivativeEntry = evaluate_triangle_integral(nodeCurrent, nodeBefore, subStart, subEnd, zBefore, zCurrent, True)
                    v1[timeIdx] += derivativeEntry
        kthDerivativeMat = np.zeros(shape=(len(v1),rateLength))
        for i in xrange(rateLength):
                kthDerivativeMat[i:, i] = v1[:(timeLength- i),0]
        jacobianConvolution[:,k] = np.transpose(kthDerivativeMat.dot(y))
    return jacobianConvolution

def evaluate_triangle_integral(nodeCurrent, nodeBefore, subStart, subEnd, zBefore, zCurrent, goingDown):
    ### integrate the derivative with respect to the response values over time
    dCurrent = evaluate_integral_derivative(nodeCurrent, nodeBefore, subStart, subEnd, zBefore, zCurrent)
    cCurrent = evaluate_integral(nodeCurrent, nodeBefore, subStart, subEnd, zBefore, zCurrent)
    if goingDown:
        factor = - 1
        xIntersect = nodeCurrent
    else:
        factor = 1
        xIntersect = nodeBefore
    result = factor * (dCurrent - xIntersect * cCurrent) / (nodeCurrent - nodeBefore)
    return result

    
def generate_secondDerivativeMatrix(nodes):
    ### calculate the matrix measuring the sinus of the angle between each interpolating
    ### function of the response estimate
    dimCol = len(nodes)
    sdMat = np.zeros(shape = ((dimCol - 1), dimCol))
    for idx in xrange((dimCol - 1)):
        angleSideAfter = nodes[idx + 1] - nodes[idx]
        if idx == 0:
            sdMat[0, 0] = - 1 / angleSideAfter
            sdMat[0, 1] = 1 / angleSideAfter
        else:
            # sign inverted? Typo in schroeter et al 2004?
            angleSideBefore =  nodes[idx] - nodes[idx - 1]
            sdMat[idx, idx + 1] = 1 / angleSideAfter
            sdMat[idx, idx] = - (angleSideAfter + angleSideBefore) / (angleSideAfter * angleSideBefore)
            sdMat[idx, idx - 1] = 1 / angleSideBefore
    return sdMat

def generate_fMatrix(rew, z, rateLength, wlLength, nodes, timeseries, wMat):
    ### calculate the matrix that will be multiplied with the presumed rates for the total error function
    ## create skeleton for convolution error measure
    convMat = - generate_convMatrix(z, rateLength, nodes, timeseries)
    #print("convMat:")
    #print(np.diagonal(convMat))
    wlNatUnit = np.ones((wlLength, 1))
    convError = np.hstack([wlNatUnit, convMat])
    ## create skeleton for rate error measure
    rateIdentity = np.identity(rateLength) * np.sqrt(rew)
    ratesError = np.hstack([np.zeros((rateLength, 1)), rateIdentity])
    ## create skeleton for smoothness measure
    smthsError = np.zeros((len(nodes) - 1,rateLength + 1))
    ## put them all together
    fMat = np.vstack([convError, ratesError, smthsError])
    fMat = wMat.dot(fMat)
    return fMat
    

def generate_convMatrix(z, rateLength, nodes, timeseries):
    ### calculate the actual convolution matrix (assuming constant rate intervals!)
    v1 = [calculate_entries(tStart, tEnd, z, nodes) for tStart, tEnd in zip(timeseries, timeseries[1:])]
    convMat = np.zeros(shape=(len(v1),rateLength))
    for i in xrange(rateLength):
            convMat[i:, i] = v1[:(len(v1)- i)]
    return convMat

def calculate_entries(start, end, z, nodes):
    ### calculate each entry of the convolution matrix which corresponds to each pumping period in time (only works for constant rate intervals!)
    #print("Calculate conv-entry for the time range:")
    #print("Start: " + str(start) + " - End: " + str(end))
    if start == 0:
        idxs = np.where(nodes <= np.log(end))[0]
    else:
        idxs = np.where((np.log(start) < nodes) & (nodes < np.log(end)))[0]
    if(len(idxs) == 0):
        warnings.warn("One pumping period does not enclose at least one node interval. The resolution of nodes might be too low!") 
        idxs = np.array([min(np.where(nodes >= np.log(end))[0])])
    elif idxs[-1] + 1 != len(nodes):
        idxs = np.append(idxs, idxs[-1] + 1)
    entry = 0
    for idx in idxs:
        idxBefore = idx - 1
        subSum = evaluate_integral(nodes[idx], nodes[idxBefore], start, end, z[idxBefore], z[idx])
        #monitor_integral(subSum, nodes[idxBefore], nodes[idx], start, end, z[idxBefore], z[idx])
        entry += subSum
    return entry
    
def evaluate_integral_derivative(nodeCurrent, nodeBefore, start, end, zBefore, zCurrent):
    ### integrate linearly amplified node-dependent part of the derivative with respect to the response values over time
    upperLim = np.minimum(np.log(end), nodeCurrent)
    if start == 0:
        lowerLim = nodeBefore
    else:
        lowerLim = np.maximum(np.log(start), nodeBefore)
    if upperLim < lowerLim:
        ## should never happen
        raise Exception("\nSomething went wrong in the integration domain!")
        #slope = 1
        #intersect = 0
        #result = np.exp(upperLim * slope) * (upperLim / slope - 1 / slope ** 2)
    else:
        slope = (zCurrent - zBefore) / (nodeCurrent - nodeBefore)
        intersect = zCurrent - slope * nodeCurrent
        if slope != 0: #round(slope, 6) != 0:
            result = ((upperLim - 1 / slope) * np.exp(slope * upperLim) - (lowerLim - 1 / slope) * np.exp(slope * lowerLim)) / slope
        else:
            result = (upperLim ** 2 - lowerLim ** 2) / 2
    result = result * np.exp(intersect)
    return result

def evaluate_integral(nodeCurrent, nodeBefore, start, end, zBefore, zCurrent):
    ### calculate the actual convolution of each entry by integrating the response function over time
    upperLim = np.minimum(np.log(end), nodeCurrent)
    if start == 0:
        lowerLim = nodeBefore
    else:
        lowerLim = np.maximum(np.log(start), nodeBefore)
    if upperLim < lowerLim:
        ## wellbore storage case
        slope = 1
        intersect = zCurrent - nodeCurrent
        result = np.exp(upperLim * slope) / slope
    else:
        slope = (zCurrent - zBefore) / (nodeCurrent - nodeBefore)
        intersect = zCurrent - slope * nodeCurrent
        if slope != 0: #round(slope, 6) != 0:
            result = (np.exp(upperLim * slope) - np.exp(lowerLim * slope)) / slope
        else:
            result = upperLim - lowerLim
    result = result * np.exp(intersect)
    return result

def monitor_integral(subSum, nodeBefore, nodeCurrent, start, end, zBefore, zCurrent):
    ### monitor function for trouble-shooting
    global subSums
    entryKey = str(start) + "..." + str(end)
    subSumResult = {
        "subSum":subSum,
        "nodeBefore":nodeBefore,
        "nodeCurrent":nodeCurrent,
        "zBefore":zBefore,
        "zCurrent":zCurrent
    }
    if entryKey in subSums:
        entryList = subSums[entryKey]
        entryList.append(subSumResult)
    else:
        entryList = [subSumResult]
    subSums[entryKey] = entryList




