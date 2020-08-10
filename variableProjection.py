def variable_projection(nodes, waterlevel, x, rates, z, sc, weights, stepsize, timeseries, errorIn = None):
    ### VP-algorithm, each cycle calculates the new z and x (i.e. pNat and y) values
    ### recursive function

    fMat = generate_fMatrix(weights["rew"], z, len(rates), len(waterlevel), nodes, timeseries)
    vVec = generate_vVector(weights["rew"], weights["dw"], z, waterlevel, rates, nodes)

    ## calculate error
    error = np.linalg.norm(fMat.dot(x) - vVec) ** 2
    if errorIn:
        errorRatio = error / errorIn
        if(errorRatio <= sc):
            return x[:(len(rates) + 1)], z, x[0]
    else:
        errorIn = error

    ## solve linear part
    constant = vVec - fMat.dot(x)
    dx = pinv(fMat).dot(constant)
    x = x + dx

    ## solve non-linear part
    currentPosition = fMat.dot(x) - vVec
    zJacobian = generate_jacobian(nodes, z, y, weights["dw"], len(rates))
    totalJacobian = np.concatenate((fMat, zJacobian), axis=1)
    dTotal = pinv(totalJacobian).dot(-currentPosition) * stepsize
    x = x + dTotal[:(len(rates) + 1)]
    z = z + dTotal[(len(rates) + 1):]

    ## call again
    variable_projection(nodes, waterlevel, x, rates, z, sc, weights, stepsize, timeseries, errorIn)

def generate_vVector(rew, dw, z, waterlevel, rates, nodes):
    ### calculate the vector that contains the waterlevel, the wheighted rates
    ### and the smoothness measure for the total error function
    smoothness = generate_smoothnessMeasure(dw, nodes, z)
    vVec = [waterlevel, rates * np.sqrt(rew), smoothness]
    return vVec

def generate_smoothnessMeasure(dw, nodes, z):
    ### calulculate the smoothness measure to minimize
    sdMat = generate_secondDerivativeMatrix(nodes)
    expected = np.zeros(len(nodes) - 1)
    expected[0] = 1
    smoothness = (sdMat.dot(z) - expected) * np.sqrt(dw)
    return smoothness

def generate_jacobian(nodes, z, y, dw, rateLength):
    ### calculate the jacobian matrix of the error measure with respect to the 
    ### response values for Gauss-Newton
    jacobianConvolution = generate_jacobianConvolution(nodes, z, y, rateLength)
    sdMat = generate_secondDerivativeMatrix(nodes) * np.sqrt(dw)
    jacobianRates = np.zeros(shape=(len(y), len(z)))
    zJacobian = np.concatenate((jacobianConvolution, sdMat, jacobianRates))
    return zJacobian

def generate_jacobianConvolution(nodes, z, y, rateLength):
    ### calculate the jacobian matrix of the convolution matrix with respect to the response values
    jacobianConvolution = np.zeros(shape=(int(np.exp(nodes[-1])), len(z)))
    for k in xrange(len(nodes)):
        v1 = np.zeros(int(np.exp(nodes[-1])))
        if k == 0:
            startTriangle = 0
            endTriangle = np.ceil(np.exp(nodes[k + 1]))
            derivativeEntry = evaluate_integral(nodes[0], nodes[-1], startTriangle, endTriangle, z[-1], z[0])
            v1[0] = derivativeEntry
            nodeCurrent = nodes[k + 1]
            nodeBefore = nodes[k]
            zCurrent = z[k + 1]
            zBefore = z[k]
            for timeIdx in xrange(startTriangle, endTriangle):
                subStart = timeIdx
                subEnd = timeIdx + 1
                derivativeEntry = evaluate_triangle_integral(nodeCurrent, nodeBefore, subStart, subEnd, zBefore, zCurrent, True)
                v1[timeIdx] += derivateEntry
        else:
            startTriangle = np.floor(np.exp(nodes[k - 1]))
            endTriangle = np.ceil(np.exp(nodes[k]))
            nodeCurrent = nodes[k]
            nodeBefore = nodes[k - 1]
            zCurrent = z[k]
            zBefore = z[k - 1]
            for timeIdx in xrange(startTriangle, endTriangle):
                subStart = timeIdx
                subEnd = timeIdx + 1
                derivativeEntry = evaluate_triangle_integral(nodeCurrent, nodeBefore, subStart, subEnd, zBefore, zCurrent, False)
                v1[timeIdx] = derivateEntry
            if k != (len(k) - 1):
                startTriangle = np.floor(np.exp(nodes[k]))
                endTriangle = np.ceil(np.exp(nodes[k + 1]))
                nodeCurrent = nodes[k + 1]
                nodeBefore = nodes[k]
                zCurrent = z[k + 1]
                zBefore = z[k]
                for timeIdx in xrange(startTriangle, endTriangle):
                    subStart = timeIdx
                    subEnd = timeIdx + 1
                    derivativeEntry = evaluate_triangle_integral(nodeCurrent, nodeBefore, subStart, subEnd, zBefore, zCurrent, True)
                    v1[timeIdx] += derivateEntry
        
        kthDerivativeMat = np.zeros(shape=(len(v1),rateLength))
        for i in xrange(rateLength):
                kthDerivativeMat[i:, i] = v1[:- i]
        jacobianConvolution[:,k] = - kthDerivativeMat.dot(y)
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
            sdMat[idx, idx] = (angleSideAfter + angleSideBefore) / (angleSideAfter * angleSideBefore)
            sdMat[idx, idx - 1] = 1 / angleSideBefore
    return sdMat

def generate_fMatrix(rew, z, rateLength, wlLength, nodes, timeseries):
    ### calculate the matrix that will be multiplied with the presumed rates for the total error function
    convMat = generate_convMatrix(z, rateLength, nodes, timeseries)
    wlNatIdentity = np.full(wlLength,1)
    rateIdentity = np.identity(rateLength) * np.sqrt(rew)
    fMat = [[wlNatIdentity, -convMat], [0, rateIdentity], [0, 0]]
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
    if start == 0:
        idxs = np.where(nodes <= np.log(end))[0]
    else:
        idxs = np.where((np.log(start) <= nodes) & (nodes <= np.log(end)))[0]
    if(len(idxs) == 0):
        warnings.warn("One pumping period does not enclose at least one node interval. The resolution of nodes might be too low!") 
        idxs = np.array([min(np.where(nodes >= np.log(end))[0])])
    else:
        idxs = np.append(idxs, idxs[-1] + 1)
    entry = 0
    for idx in idxs:
        idxBefore = idx - 1
        subSum = evaluate_integral(nodes[idx], nodes[idxBefore], start, end, z[idxBefore], z[idx])
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
        slope = 1
        intersect = 0
        result = np.exp(upperLim * slope) * (upperLim / slope - 1 / slope ** 2)
    else:
        slope = (zCurrent - zBefore) / (nodeCurrent - nodeBefore)
        intersect = zCurrent - slope * nodeCurrent
        if round(slope, 6) != 0:
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
        intersect = 0
        result = np.exp(upperLim * slope) / slope
    else:
        slope = (zCurrent - zBefore) / (nodeCurrent - nodeBefore)
        intersect = zCurrent - slope * nodeCurrent
        if round(slope, 6) != 0:
            result = (np.exp(upperLim * slope) - np.exp(lowerLim * slope)) / slope
        else:
            result = upperLim - lowerLim
    result = result * np.exp(intersect)
    return result




