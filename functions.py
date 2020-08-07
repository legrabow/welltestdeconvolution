def get_nodes(amountNodes, interpolation, startNode, timeLength):
	### calculate and return the nodes corresponding to the logarithmic time steps of the response function
	n0 = np.ln(startNode)
	nEnd = np.ln(timeLength)
	nodes = [n0]
	if interpolation == "linear":
		## use linear interpolation scheme, between the start - and end node
		for k in xrange(amountNodes):
			n = n0 + k * (nEnd - n0) / (amountNodes - 1)
			nodes.append(n)
	nodes.append(nEnd)
	return nodes

def get_initial_wlNat(waterlevel):
	### return the maximum waterlevel value as the natural waterlevel
	message = "\nNo inital guess for the natural water level given." \
	"The maximum water level in the time range will be taken instead." \
	"\nNote that this makes only sense if recovery periods exist in the time range!"
	warnings.warn(message)
	wlNat = max(waterlevel)
	return wlNat

def get_initial_responses(nodes, waterlevel, rates, pNat):
	### calculate and return an inital response function as a start value for the
	### iterative VP-algorithm
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
	return responsesTotal
