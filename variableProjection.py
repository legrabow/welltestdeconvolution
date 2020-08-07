import numpy as np

def variable_projection(nodes, waterlevel, rates, pNat, z, sc, weights):
	# solve linear part
	# solve non-linear part
	if(lower stoping criterion):
		variable_projection()
	else:
		return()

def generate_vVector(rew, dw, z, waterlevel, rates):
	### calculate the vector that contains the waterlevel, the wheighted rates
	### and the smoothness measure for the total error function
	smoothness = generate_smoothnessMeasure()
	vVec = [waterlevel, rates * np.sqrt(rew), smoothness]
	return(vVec)

def generate_smoothnessMeasure(dw):
	### calulculate the smoothness measure to minimize
	sdMat = generate_secondDerivativeMatrix(nodes)
	expected = np.zeros(len(nodes))
	expected[0] = 1
	smoothness = (sdMat.dot(z) - expected) * np.sqrt(dw)
	return(smoothness)

def generate_jacobian():
	### calculate the jacobian matrix of the error measure with respect to the 
	### response values for Gauss-Newton
	for k in xrange(nodes):
		for i in 
	
def generate_secondDerivativeMatrix(nodes):
	### calculate the matrix measuring the sinus of the angle between each interpolating
	### function of the response estimate
	dimRow = len(nodes)
	sdMat = np.zeros(shape = (dimRow, dimRow)
	for idx in xrange(dimRow):
		angleSideAfter = nodes[idx + 1] - nodes[idx]
		if idx == 0:
			sdMat[0, 0] = - 1 / angleSideAfter
			sdMat[0, 1] = 1 / angleSideAfter
		else:
			angleSideBefore =  nodes[idx] - nodes[idx - 1]
			sdMat[idx, idx - 1] = 1 / angleSideAfter
			sdMat[idx, idx] = (angleSideAfter + angleSideBefore) / (angleSideAfter * angleSideBefore)
			sdMat[idx, idx + 1] = 1 / angleSideBefore
	return(sdMat)

def generate_fMatrix(rew, z, rateLength, wlLength):
	### calculate the matrix that will be multiplied with the presumed rates for the total error function
	convMat = generate_convMatrix(z, rateLength)
	pNatIdentity = np.identity(wlLength)
	rateIdentity = np.identity(rateLength) * np.sqrt(rew)
	fMat = [[pNatIdentity, -convMat], [0, rateIdentity], [0, 0]]
	return(fMat)
	

def generate_convMatrix(z, rateLength):
	### calculate the actual convolution matrix (assuming constant rate intervals!)
	v1 = [calculate_entries(start = t - 1, end = t, z) for t in 1:int(np.exp(nodes[-1]))]
	convMat = np.zeros(shape=(len(v1),rateLength))
	for i in xrange(rateLength):
    		convMat[i:, i] = v1[:- i]
	return(convMat)

def calculate_entries(start, end, z, nodes):
	### calculate each entry of the convolution matrix which corresponds to each pumping period in time (only works for constant rate intervals!)
	if start == 0:
		idxs = np.where(nodes <= np.ln(end))[0]
	else:
		idxs = np.where((np.ln(start) <= nodes) & (nodes <= np.ln(end)))[0]
	if(len(idxs) == 0):
		raise Exception("One pumping period does not enclose at least one node interval. The resolution of nodes is too low!")
	idxs = np.append(idxs, idxs[-1] + 1)
	entry = 0
	for idx in idxs:
		idxBefore = idx - 1
		subSum = evaluate_integral(nodeCurrent = nodes[idx], nodeBefore = nodes[idxBefore], start, end, zBefore = [idxBefore], zCurrent = z[idx])
		entry += subSum
	return(entry)
	
def evaluate_integral_derivative(nodeCurrent, nodeBefore, start, end, zBefore, zCurrent):
	### integrate the node-dependent part of the derivative with respect to the response values over time
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
		if round(slope, 4) != 0:
			result = ((upperLim - 1 / slope) * np.exp(slope * upperLim) - (lowerLim - 1 / slope) * np.exp(slope * lowerLim)) / slope
		else:
			result = (upperLim ** 2 - lowerLim ** 2) / 2
	return(result)

def evaluate_integral(nodeCurrent, nodeBefore, start, end, zBefore, zCurrent):
	### calculate the actual convolution of each entry by integrating the response function over time
	upperLim = np.minimum(np.log(end), nodeCurrent)
	if start == 0:
		lowerLim = nodeBefore
	else:
		lowerLim = np.maximum(np.log(start), nodeBefore)
	if upperLim < lowerLim:
		slope = 1
		intersect = 0
		result = np.exp(upperLim * slope) / slope
	else:
		slope = (zCurrent - zBefore) / (nodeCurrent - nodeBefore)
		intersect = zCurrent - slope * nodeCurrent
		if round(slope, 4) != 0:
			result = (np.exp(upperLim * slope) - np.exp(lowerLim * slope)) / slope
		else:
			result = upperLim - lowerLim
	result = result * np.exp(intersect)
	return(result)




