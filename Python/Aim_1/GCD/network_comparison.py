#!/usr/bin/env python3

"""		
	Purpose of the script:
	----------------------
	The code for computing alignment-free network distances between a set of network given in 
	a folder. 
	
	The script automatically searches all relevant network files for the computation of the 
	network distances in a given folder. 
	
	With the script, the following network distances can be computed:	
		1) GCD distance
		2) RGF Distance
		3) GDD-Agreement with Arithmetic Mean (GDD-A)
		4) GDD-Agreement with Geometric Mean (GDD-G)
		5) Degree Distribution Distance
		6) Average Degree Distance
		7) Average Clustering Coefficient Distance
		8) Diameter Distances of the network
	
	The script outputs the computed distances between all pairs of networks in CSV format
		
	Notes on using the script:
	--------------------------
	
	The script excepts a folder. Using python's os.walk() function, it automatically 
	identifies all networks under the given folder.
	
	Each network files should contain a set of networks (in LEDA format -- ending 
	with extension .gw) and their graphlet degree vector (.ndump2) files.
	
	The code processes them through the following steps:
	
	1) Reads the networks from '.gw' files and the graphlet signatures from '.ndump' files
	2) Depending on the requested distance measure, computes the relevant network properties 
	(e.g., degree distribution, diameter, clustering coefficient) if needed
	3) For computing graphlet degree Compute the log-scaled signatures of all nodes in all files in the folder (log(gd + 1) for all graphlet degrees)
	
	Run as:
		python network_comparison.py <network_folder> <distance_type> <process_count>
		
		<network_folder> : The folder in which all networks are going to be compared.
				Contains the '.ndump2' and '.gw' files of the networks in the cluster
				The content of the networks folder is as follows:
					1) '.ndump2' files - contains the graphlet signatures of all the nodes
					2) '.gw' files - contains the network in LEDA format (required only in <test_mode> == 3)
				The names of the '.gw' files and '.ndump2' files should exactly match.
		
		<distance_type>: 
			'rgf'			- RGF distance
			'gdda'			- GDD Agreement (Both Arithmetic \& Geometric)
			'degree' 		- Degree distribution & Average degree distances
			'clustering'	- Clustering Coefficient
			'diameter' 		- Diameter
			'gcd11'  		- Graphlet Correlation distance with non-redundant 2-to-4 node graphlet orbits
			'gcd15'  		- Graphlet Correlation distance with all 2-to-4 node graphlet orbits
			'gcd58'  		- Graphlet Correlation distance with non-redundant 2-to-5 node graphlet orbits
			'gcd73'  		- Graphlet Correlation distance with all 2-to-5 node graphlet orbits
			'spectral'		- Spectral distance using the eigenvalues of the Laplacian representation of the network
		
		<process_count>:
			Any number higher than or equal to 1. Determines the number of processes to use for computing the 
			distances.
			
	Implemented by:
		Omer Nebil Yaveroglu
		
		First implementation 	= 01.06.2012 - 11:53
		Revision 1 				= 10.08.2012 - 16:48
		Revision 2 				= 23.08.2012 - 11:34
		Parallelized			= 09.10.2012 - 17:00
		Clean Version 			= 06.05.2014 - 15:21
		Python 3 Update         = 13.09.2024 - 12:00
"""

import sys
import os
import math
import numpy
import time
import networkx as nx
import multiprocessing
import queue

from scipy import stats

"""
	Functions
"""

# Read the signatures from ndump2 files
def readSignatures(file):
	signDict = []
	
	with open(file, 'r') as fRead:
		for line in fRead:
			signDict.append([float(value) for value in line.strip().split(' ')[-73:]])
	
	return signDict

# Remove the redundant orbits and return the log scaled graphlet degrees
def formatSignatures(signList, testMode):
	formattedSignatures = []
	
	for sign in signList:
		# Eliminate the orbits that we are not interested
		if testMode == 16:		# GCD-15
			log = sign[:15]
		elif testMode == 14:	# GCD-73
			log = sign
		elif testMode == 7:		# GCD-58
			eliminateList = [3, 5, 7, 14, 16, 17, 20, 21, 23, 26, 28, 38, 44, 47, 69, 71, 72]
			log = [sign[i] for i in range(73) if i not in eliminateList]
		elif testMode == 10:	# GCD-11
			eliminateList = [3, 12, 13, 14]
			log = [sign[i] for i in range(15) if i not in eliminateList]
		
		formattedSignatures.append(log)
	
	return formattedSignatures

# Compute the correlation matrix without isnan values by adding a dummy signature
def computeCorrelMat(formattedSigns):
	
	length = len(formattedSigns[0])
	
	# Add the dummy signature for some noise
	formattedSigns.append([1] * length)
	
	# Compute the ranking for the Spearman's correlation coefficient computation
	rankList = []
	for i in range(length):
		rankList.append(stats.mstats.rankdata([val[i] for val in formattedSigns]))
		
	correlMat = numpy.corrcoef(rankList, rowvar=1)
		
	return correlMat


# The parallel reading class to compute the orbit correlation matrices depending on the test mode
class MatrixReader(multiprocessing.Process):
	def __init__(self, work_queue, result_queue, testMode):
		multiprocessing.Process.__init__(self)
		
		self.work_queue 	= work_queue
		self.result_queue 	= result_queue
		self.testMode 		= testMode
		self.kill_received 	= False
		
	
	def run(self):
		while not self.kill_received:
			# Get a task
			try:
				ndumpName = self.work_queue.get_nowait()
				
				signatures = readSignatures(f'{ndumpName}.ndump2')
				formatted = formatSignatures(signatures, self.testMode)				
				correlMat = computeCorrelMat(formatted)
								
				self.result_queue.put((ndumpName, correlMat))
				
			except queue.Empty:
				pass
			

# Computes the orbit correlation matrices for all the correlation matrices provided in allIndexes 
def getCorrelationMatrices(allIndexes, testMode):
	# Prepare the list of files to be processed
	file_queue = multiprocessing.Queue()
	result_queue = multiprocessing.Queue()
	
	processList = []
	for i in range(num_processes):
		reader = MatrixReader(file_queue, result_queue, testMode)
		reader.start()
		processList.append(reader)
	
	# Put the jobs to be consumed
	jobCount = len(allIndexes)
	submitCount = 0
	
	for index in allIndexes:
		file_queue.put(index)
		submitCount += 1
		
		#if submitCount % 100 == 0:
		#	print(f'Submitted correlation computation for: {float(submitCount) / jobCount * 100}%')

	# Process the results of computation
	correlMats = {}
	
	finishedCount = 0
	while finishedCount < len(allIndexes):
		try:
			matrix = result_queue.get_nowait()
			correlMats[matrix[0]] = matrix[1]
			finishedCount += 1
		except queue.Empty:
			time.sleep(1)
		
		#if finishedCount % 100 == 0:
		#	print(f'Finished reading: {float(finishedCount) / jobCount * 100}%')
	
	for proc in processList:
		proc.terminate()
	
	return correlMats

# Computes the euclidean distance between two correlation matrices
def computeMatrixDist(matrix1, matrix2):
	differenceSum = 0
	
	for i in range(len(matrix1) - 1):
		for j in range(i + 1, len(matrix1)):
			differenceSum += pow(matrix1[i][j] - matrix2[i][j], 2)
	
	eucDist = math.sqrt(differenceSum)
	
	return eucDist

# The parallel reading class to compute the orbit correlation distances
class correlDistanceComputer(multiprocessing.Process):
	def __init__(self, work_queue, result_queue):
		multiprocessing.Process.__init__(self)
		
		self.work_queue = work_queue
		self.result_queue = result_queue
		self.kill_received = False
	
	def run(self):
		while not self.kill_received:
			# Get a task
			try:
				# matrixPair : 0,1 holds names; 2, 3 holds matrices
				matrixPair = self.work_queue.get_nowait()
				distance = computeMatrixDist(matrixPair[2], matrixPair[3])
				self.result_queue.put((matrixPair[0], matrixPair[1], distance))
			except queue.Empty:
				pass


# Given a matrix, writes the matrix with the network names into the output file
def saveDistanceMatrix(matrix, networkNames, outputFile):
	with open(outputFile, 'w') as fWrite:
		# Write the names of the networks
		fWrite.write('\t' + '\t'.join(networkNames) + '\n')
		
		# Write the distances among networks
		for i, name in enumerate(networkNames):
			fWrite.write(f"{name}\t" + '\t'.join(map(str, matrix[i])) + '\n')

# The function to compute all the distances between the provided correlation matrices in parallel
def computeCorrelDist(corrMats, outputName):
	# Start the processes
	pair_queue = multiprocessing.Queue()
	result_queue = multiprocessing.Queue()
	processList = []
	
	for i in range(num_processes):
		computer = correlDistanceComputer(pair_queue, result_queue)
		computer.start()
		processList.append(computer)
	
	# Put the jobs to be consumed
	totalJobCount = len(corrMats) * (len(corrMats) - 1) // 2
	matList = list(corrMats.keys())
	matValList = [corrMats[mat] for mat in matList]
	
	pairCount = 0
	for i in range(len(matValList) - 1):
		corrMat1 = matValList[i]
		
		for j in range(i+1, len(matValList)):
			corrMat2 = matValList[j]
			
			pair_queue.put((i, j, corrMat1, corrMat2))
			pairCount += 1
	
	# Consume the results of completed computation
	distances = [[0] * len(corrMats) for i in range(len(corrMats))]
	
	computedCount = 0
	while computedCount < pairCount:
		try:
			results = result_queue.get_nowait()
			distances[results[0]][results[1]] = results[2] 
			distances[results[1]][results[0]] = results[2]
			computedCount += 1
		except queue.Empty:
			time.sleep(1)
		
	for proc in processList:
		proc.terminate()
	
	# Save the results in the output file
	saveDistanceMatrix(distances, matList, outputName)

# Function to compute the graphlet counts from ndump2 files
def getGraphletFreq(signList):
	orbits = [2, 3, 5, 7, 8, 9, 12, 14, 17, 18, 23, 25, 27, 33, 34, 35, 39, 44, 45, 50, 52, 55, 56, 61, 62, 65, 69, 70, 72]
	weights = [1, 3, 2, 1, 4, 1, 2, 4, 1, 1, 1, 1, 1, 1, 5, 1, 1, 1, 1, 2, 1, 2, 1, 1, 1, 1, 1, 2, 5]
	
	# Derive the graphlet counts from the orbit degrees
	graphletCounts = []
	
	for i in range(len(orbits)):
		orbit = orbits[i]
		sumCount = sum(val[orbit] for val in signList)
		graphletCounts.append(sumCount / weights[i])
	
	return graphletCounts

# Normalize and scale the graphlet distributions for the computation of GDD Agreement
def scaleGraphletDists(signatures):
	distributions = []
	
	for i in range(73):
		# Get the distribution
		values = {}
		for val in signatures:
			values[val[i]] = values.get(val[i], 0) + 1
	
		if 0 in values:
			del values[0]
		
		# Scale the distribution values for GDD agreement
		total = sum(float(values[val]) / val for val in values)
	
		# Normalize the distributions
		for val in values:
			values[val] = (float(values[val]) / val) / total
			
		distributions.append(values)
	
	return distributions

# Write the distributions for the network
def writeDistributions(outputName, distribution):
	with open(outputName, 'w') as fWrite:
		for dictin in distribution:
			toPrint = ','.join(f"{val}_{dictin[val]}" for val in dictin)
			fWrite.write(toPrint + '\n')

# The parallel running class for reading the graphlet counts from ndump files
class GraphletCountGetter(multiprocessing.Process):
	def __init__(self, work_queue, result_queue, mode):
		multiprocessing.Process.__init__(self)
		
		self.work_queue = work_queue
		self.result_queue = result_queue
		self.mode = mode
		self.kill_received = False
	
	def run(self):
		while not self.kill_received:
			# Get a task
			try:
				ndumpName = self
