import numpy as np

def get_nodes(amountNodes, interpolation, startNode, timeLength):
	n0 = np.ln(startNode)
	nEnd = np.ln(timeLength)
	nodes = [n0]
	if interpolation == "linear":
		for k in xrange(amountNodes):
			n = n0 + k * (nEnd - n0) / (amountNodes - 1)
			nodes.append(n)
	nodes.append(nEnd)
	return(nodes)

def get_initial_wlNat(waterlevel):
	wlNat = max(waterlevel)
	return(wlNat)

def get_initial_responses(nodes, waterlevel, rates, pNat):
	### first guess: well bore storage before and radial flow after the first node
	## calculate the best fit for a radial flow
	radialFlow = np.full(shape = (len(waterlevel),len(rates)), fill_value = 1)
	radialFlow = np.tril(radialFlow)
	convolution = radialFlow.dot(rates)
	drawdown = pNat - waterlevel
	rfCoef = drawdown.dot(convolution) / np.linalg.norm(convolution) ** 2
	## calculate well bore storage coefficient
	wbsCoef = rfCoef * np.exp(-nodes[0])
	## create response array
	responsesRadial = np.full(np.ln(radialFlow), len(nodes) - 1)
	responsesTotal = np.insert(responsesRadial, 0, ln(wbsCoef))

	return(responsesTotal)
