import numpy as np

def variable_projection(nodes, waterlevel, rates, pNat, z, sc, weights):
	# solve linear part
	# solve non-linear part
	if(lower stoping criterion):
		variable_projection()
	else:
		return()

def generate_Fmatrix(rew, z, rateLength):
	convMat = generate_convMatrix(z, rateLength)
	

def generate_convMatrix(z, rateLength):
	v1 = [calculate_entries(start = t - 1, end = t, z) for t in 1:int(np.exp(nodes[-1]))]
	convMat = np.zeros(shape=(len(v1),rateLength))
	for i in xrange(rateLength):
    		convMat[i:, i] = a[:(len(v1) - i)]
	return(convMat)

def calculate_entries(start, end, z, nodes):
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
	

def evaluate_integral(nodeCurrent, nodeBefore, start, end, zBefore, zCurrent):
	upperLim = np.minimum(np.log(end), nodeCurrent)
	if start == 0:
		lowerLim = nodeBefore
	else:
		lowerLim = np.maximum(np.log(start), nodeBefore)
	if upperLim < lowerLim:
		slope = 1
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




