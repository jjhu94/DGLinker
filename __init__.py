# -*- coding: utf-8 -*-
import matplotlib																										 
matplotlib.use('Agg')

from collections import OrderedDict
from collections import defaultdict
from collections import Counter
from flask import flash, Flask, make_response, redirect, render_template, request, send_file, url_for
from pyvis.network import Network
from threading import Thread
from werkzeug.utils import secure_filename

import csv
import json
import getpass
import igraph
import json
import numpy as np
import math
import matplotlib.pyplot as plt
import os
import pandas as pd 
import random
import requests
import shutil
import smtplib, ssl
import subprocess
import sys
import time
import uuid
import yagmail
import zipfile


import sys, csv, json, igraph, itertools
import numpy as np
import scipy.stats as stats #for stats.fisher_exact
import yagmail
from Objective import Objective

#for multiple testing correction
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector


ALLOWED_EXTENSIONS = set(['csv'])


app = Flask(__name__)

app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

def allowed_file(filename):
	return '.' in filename and filename.rsplit('.', 1)[1] in ALLOWED_EXTENSIONS

def send_error(email, tempid, error1, error2):
	if email != '':
		yag = yagmail.SMTP()
		subject = "DGLinker: There is something wrong with your project."
		contents = "Dear user,\n\n  There is something wrong with your project： \n https://dglinker.rosalind.kcl.ac.uk/result/%s \n" % tempid
		contents = contents + str(error1) + '\n'+ str(error2)
		contents = contents + '\n Your project will not continue. Please submit a new one. \n\n Best wishes \n DGLinker Service'
		yag.send(email, subject, contents)

def error_page(tempid, error1, error2):
	f1 = open('/var/www/kgp/FlaskApp/cache/error_file%s.txt' % tempid,'a')
	f1.write(error1)
	f1.write('\n')
	f1.write(error2)

class EdgePrediction:
	"""Edge Prediction class
	Implements the aglorithm described in Bean et al. 2017
	All parameters have sensible defaults except to_predict, which must be specified by the user at some point
	before the algorithm can run.
	Parameters
	----------
	min_weight : float
		minimum allowed feature weight, default 0.0
	max_weight : float
		maximum allowed feature weight is max_weight - step. Default 1.1 gives a  max feature weight is 1.0 with the 
		default step size 0.1
	
	step : float
		step size used in parameter grid search for feature weights. Default 0.1
	
	to_predict : str
		Type of edge to predict. Must exactly match an edge type in the graph and must be set by user
		before the prediction algorithm will run. default None.
	
	pval_significance_threshold : float
		threshold applied to select features from the enrichment test. Features with p-value for enrichment (after
		multiple testing correction if used) < threshold are considered enriched. Default 0.05
	
	require_all_predictors : bool
		Optional. If true, a model will only be trained for a given target if at least one predictor for each type of edge 
		in the graph is found with the enrichment test (after applying multiple testing correction if used). Default True
	
	objective_function : str
		Optional. This is parameterised to allow convenient extension. Only the J objective is used in Bean et al 2017.
		J means Youden's J statistic, J = sensitivity + specificity - 1. Default is "J", options are {'J', 'F1', 'F2', 'F05', 'ACC'}.
	
	ties : str
		Optional. Method used to select which set of weights to use where two sets give identical performance of the objective 
		function. 'first' uses the set found first, 'minL2norm' uses L2 normalisation to prefer balanced weights. If
		two sets of weights give idential performance and have the same L2 normalisation, the weights found first are
		kept. This means the network_order is important for both methods. Default 'minL2norm', options are {'first', 'minL2norm'}
	network_order : list
		Optional. The order in which the features are iterated over. This can be important as ultimately any ties are broken by keeping
		the parameters found first in the search. The best order may not ultimately matter, and will be context-specific. In
		general it is recommended to specify the order so it will remain consistent. The default order is determined by the 
		keys of the internal dict object, which is not guaranteed. 
	
	randomise_folds : bool
		Optional. Whether to randomise the order of items in the full dataset before splitting all items into train-test folds.
		If False, the folds will always be identical between runs. Default True
	correct_pval : str
		Optional. Correct p-values from enrichment test for multiple testing. Currently setting anything other than "BH" will
		result in no correction being applied. 'BH' or None, default 'BH'.
	Attributes
	----------
	graphs : dict 
		Internal representation of the input graph. There is one element per edge type. {'graph': igraph.Graph, 'sourcenodes': list, 
		'sourcetype': str, 'targetnodes': list, 'targettype': str}
	
	can_analyse : bool
		Flag used to keep track of any conditions that mean the network cannot be analysed
	
	optimisation_method : str
		Currently only "graph" is available, which implements the method of Bean et al. 2017. 
	"""
	def __init__(self, email, tempid, to_predict = None):
		
		self.min_weight = 0.0 
		self.max_weight = 1.1
		self.step = 0.1
		self.to_predict = None 
		self.email = email
		self.tempid = tempid
		
		self.pval_significance_threshold = 0.05
		self.require_all_predictors = False #build a model even if predictors are not found in some input graphs 
		self.objective_function = "J"
		
		self.ties = "minL2norm" #method to break tied weights
		self.randomise_folds = True #put known source nodes in random order before generating folds
		self.network_order = None

		#internal parameters
		self.graphs = {}
		self.can_analyse = True
		self.correct_pval = "BH"
		self.optimisation_method = "graph"

	def CSV_to_graph(self,fname, srcNameCol = "Source node name", 
				  srcTypeCol = "Source node type", tgtNameCol = "Target node name", 
				  tgtTypeCol = "Target node type", edgeTypeCol = "Relationship type"):
		"""Parse csv file to internal graph representation
		
		The parsed graph is stored internally in self.graphs and is not returned. 
		Parameters
		----------
		fname : str
			Input file name or path. Must be a csv file with a header.
		srcNameCol : str
			Column in input file containing source node names.
		srcTypeCol : str
			Column in input file containing source node type. 
		tgtNameCol : str
			Column in input file containing target node names. 
		tgtTypeCol : str
			Column in input file containing target node types. 
		edgeTypeCol : str
			Column in input file containing edge types.
		Returns
		-------
		bool : bool 
			True for success, False otherwise.
		"""

		edge_types = {}
		with open(fname, 'rU') as f:
			reader = csv.reader(f)
			header = next(reader, None)
			for line_raw in reader:
				line = dict(zip(header, line_raw))
				etype = line[edgeTypeCol]
				if not etype in edge_types:
					edge_types[etype] = set()
				edge_types[etype].add((line[srcNameCol], line[srcTypeCol], line[tgtNameCol], line[tgtTypeCol]))
		
		for etype in edge_types:
			sourcenodes = set()
			targetnodes = set()
			sourcetypes = set()
			targettypes = set()
			edges = set()
			for edge in edge_types[etype]:
				if edge[0] != "" and edge[2] != "":
					edges.add((edge[0], edge[2]))
					sourcenodes.add(edge[0])
					targetnodes.add(edge[2])
					sourcetypes.add(edge[1])
					targettypes.add(edge[3])
			if len(sourcetypes) > 1 or len(targettypes) > 1:
				print "ERROR: too many source or target node types"
				email = self.email
				tempid = self.tempid
				error1 = "\"ERROR: too many source or target node types.\""
				error2 = "One possible reason can be that databases has more than one source node type. Please make sure all of the source node types are 'gene'."
				send_error(email, tempid, error1, error2)
				error_page(tempid, error1, error2)
				return render_template('error.html', tempid=tempid, error1=error1, error2=error2)
			else:
				sourcetype = list(sourcetypes)[0]
				targettype = list(targettypes)[0]
			g = igraph.Graph(directed = True)
			nodes = list(sourcenodes) + list(targetnodes)
			g.add_vertices(nodes)
			nodeProp = [sourcetype]*len(sourcenodes) + [targettype]*len(targetnodes)
			g.vs['type'] = nodeProp
			g.add_edges(list(edges))
			data = {'graph': g, 'sourcenodes': list(sourcenodes), 'sourcetype': sourcetype, 'targetnodes': list(targetnodes), 'targettype': targettype}
			self.graphs[etype] = data
		return True

	def preprocess(self):
		"""Automates the first two steps to prepare data for training loop
		Does not need to be called manually.
		Parameters
		----------
		self : object
		Returns
		-------
		None : 
			Internal network representation is updated
		"""
		self.filterNetworksByCommonSource()
		self.sparseAdjacency()


	## filter all graphs, keeping only source nodes common to all graphs (with all of their edges)
	## filtering is not necessary for the algorithm, but the current implementation expects filtering and won't work otherwise
	def filterNetworksByCommonSource(self):
		"""Delete source nodes that don't have at least one edge of every type.
		Parameters
		----------
		self : object
		Returns
		-------
		None :
			Internal network representation is updated. 
		"""
		graph_names = self.graphs.keys()

		#get the common source nodes
		common_source_nodes = set(self.graphs[graph_names[0]]['sourcenodes'])
		if len(graph_names) > 1:
			for nw in graph_names[1:]:
				common_source_nodes = common_source_nodes.intersection(set(self.graphs[nw]['sourcenodes']))
		if len(common_source_nodes) == 0:
			self.can_analyse = False
			print "ERROR no source nodes are common to all input graphs. Cannot continue"
			error1 = "\"ERROR: no source nodes are common to all input graphs.\""
			error2 = "One possible reason can be that the input genes are not present in any of the selected databases. Please try to select more or different databases and/or more or other genes upon submission."
			email = self.email
			tempid = self.tempid
			send_error(email, tempid, error1, error2)
			error_page(tempid, error1, error2)
			return render_template('error.html', tempid=tempid, error1=error1, error2=error2)
		
		print "%s source nodes are common to all %s input graphs" % (len(common_source_nodes), len(graph_names))
		self.n_source_nodes = len(common_source_nodes)

		for network_name in self.graphs:
			networkA = self.graphs[network_name]
			A_source_nodes = set(networkA['sourcenodes'])
			#get a list of the drug nodes in A that are not in the common set
			delete_from_A = A_source_nodes.difference(common_source_nodes)

			#get the indices of these nodes in this graph
			A_graph_del_idx = [v.index for v in networkA['graph'].vs() if v['name'] in delete_from_A]

			#delete these nodes
			print "deleting %s nodes that don't overlap between networks" % len(A_graph_del_idx)
			networkA['graph'].delete_vertices(A_graph_del_idx)

			#check for nodes that now have degree zero
			networkA['graph'].vs['degree'] = networkA['graph'].degree()
			A_graph_del_idx = [v.index for v in networkA['graph'].vs() if v['degree'] == 0]
			#delete these nodes
			print "deleting another %s nodes that now have degree zero" % len(A_graph_del_idx)
			networkA['graph'].delete_vertices(A_graph_del_idx)

			#reset sourcenodes and targetnodes
			networkA['sourcenodes'] = networkA['graph'].vs(type_eq=networkA['sourcetype'])['name']
			networkA['targetnodes'] = networkA['graph'].vs(type_eq=networkA['targettype'])['name']

		#the network objects are mutable, so they'll be modified in the calling scope
		#return None so that this is clear
		return None

	def sparseAdjacency(self):
		"""Efficient representation of sparse adjacency matrix.
		Updates self.graphs generated from csv input with a sparse adjacency matrix. Edges are stored in both directions:
		source to target (ST) and target to source (TS). The representation is a dict where keys are node names and values are
		sets of other nodes connected with an edge of each type. There is one dict per edge type in the input data.
		Parameters
		----------
		self : object
		Returns
		-------
		None : 
			Internal network representation is updated.
		"""
		for network_name in self.graphs:

			#map graph ids to node names
			graph_names = self.graphs[network_name]['graph'].vs()['name']
			graph_idToName = dict(zip(range(0,len(graph_names)), graph_names))
			
			adj_TS = {}
			for t in self.graphs[network_name]['targetnodes']:
				adj = self.graphs[network_name]['graph'].neighbors(t, "IN")
				adj_TS[t] = set([graph_idToName[x] for x in adj])


			adj_ST = {}
			for s in self.graphs[network_name]['sourcenodes']:
					adj = self.graphs[network_name]['graph'].neighbors(s, "OUT")
					adj_ST[s] = set([graph_idToName[x] for x in adj])

			self.graphs[network_name]['ST'] = adj_ST
			self.graphs[network_name]['TS'] = adj_TS

	def groupSparseAdjacency(self, target):
		"""Adjacency for all nodes with known edges to the target vs all others.
		Parameters
		----------
		target : str
			The name of the target node that we're predicting edges to.
		Returns
		-------
		grouped : dict
			The grouped adjacency matrix. Each element of the dict is one type of edge in the network.
			The output is the full (sparse) matrix.
		"""
		known = self.graphs[self.to_predict]['TS'][target]
		grouped = {}
		for network_name in self.graphs:
			sparseAdj = self.graphs[network_name]['TS']
			nrows = len(sparseAdj.keys())
			output_matrix = np.zeros(shape = (nrows, 2), dtype=int)
			rownames = sparseAdj.keys()
			known_adj = [len(known.intersection(sparseAdj[x])) for x in rownames]
			other_adj = [len(sparseAdj[x]) - len(known.intersection(sparseAdj[x])) for x in rownames]
			output_matrix[:,0] = known_adj
			output_matrix[:,1] = other_adj
			grouped[network_name] = {'matrix':output_matrix, 'rownames':rownames}

		return grouped

	def filterSparseAdjacency(self, pvals, ignore = None):
		"""Filter a sparse adjacency matrix, keeping only the target nodes that are significantly enriched
		Parameters
		----------
		pvals : dict
			Output from from self.enrichment, see return value for self.enrichment
		ignore : bool
			name of the target node that that edges are predicted for, so it should removed from the enrichment calculation
		Returns
		-------
		all_filtered : dict
			keys are edge types, values are {'overlap':list,'colnames':list, 'predictors': list}
			'overlap' : adjacency of each source node with all predictors
			'colnames' : source nodes in the graph
			'predictors' : all enriched predictor names
		"""
		all_filtered = {}
		for network_name in self.graphs:
			predictors = pvals[network_name][:,0] < self.pval_significance_threshold
			rownames = np.array(self.graphs[network_name]['TS'].keys())
			if ignore != None and network_name == self.to_predict:
				if ignore in rownames:
					predictors[rownames == ignore] = False
					print "ignoring %s as a predictor in network %s" % (ignore, network_name)
				else:
					print "ERROR %s not found in row names for network %s, continuing" % (ignore, network_name)
			predictors = set(rownames[predictors])
			colnames = self.graphs[network_name]['ST'].keys()
			overlaps = [len(predictors.intersection(self.graphs[network_name]['ST'][x])) for x in colnames]

			all_filtered[network_name] = {'overlap':overlaps,'colnames':colnames, 'predictors': predictors}
		return all_filtered

	def enrichment(self, grouped, n_known, n_other):
		"""Fisher's exact test for enrichment to identify features (predictors)
		Parameters
		----------
		grouped : dict 
			output from from self.groupSparseAdjacency
		n_known : int
			number of source nodes with an edge to the target node
		n_other : int
			number of source nodes without an edge to the target node
		Returns
		-------
		all_pvals : dict
			Keys are edge types, values are numpy arrays. Array columns are [p, known_present, other_present, known_absent, other_absent ]
		"""
		all_pvals = {}
		for network_name in grouped:
			nrows = grouped[network_name]['matrix'].shape[0]
			pvals = np.zeros(shape=(nrows,5))
			for i in range(nrows):
				known_present = grouped[network_name]['matrix'][i,0]
				other_present = grouped[network_name]['matrix'][i,1]
				known_absent = n_known - known_present
				other_absent = n_other - other_present
				#				  known		  other
				# edge present [known_present, other_present],
				# edge absent  [known_absent,  other_absent ]
				#
				odds, p = stats.fisher_exact([ [known_present, other_present], [known_absent, other_absent] ], alternative = "greater")
				pvals[i,:] = [p, known_present, other_present, known_absent, other_absent ]
			if self.correct_pval == "BH":
				r_stats = importr('stats')
				p_list = pvals[:,0].tolist()
				p_adjust = r_stats.p_adjust(FloatVector(p_list), method = self.correct_pval)
				pvals[:,0] = list(p_adjust)
			else:
				print "WARNING - NOT correcting p-values for multiple comparisons"
			all_pvals[network_name] = pvals
		return all_pvals

	def createWeightsGenerator(self, min_weight = None, max_weight = None, step = None):
		"""Generate weights for parameter grid search.
		All parameters are required. They must be set on self for some omtimisation methods (min_weight, max_weight, step). 
		Parameters
		----------
		min_weight : float
			lower bound of search space
		max_weight : float
			upper bound of search space, should be intended bound + step
		step : float
			granularity of search space
		Returns
		-------
		weights_generator : generator
			Instance of a generator that returns all combinations of parameters in the specified range
		"""
		weights_generator = itertools.product(np.arange(min_weight, max_weight, step), repeat = len(self.graphs))
		return weights_generator
		
	
	def getKnown(self, target):
		"""Convenience function to list all nodes with an edge to the target.
		self.to_predict must be set to a valid edge type.
		Parameters
		----------
		target: str
			Node name.
		Returns
		-------
		list : list
			all nodes with an edge to the target
			returns empty list if the target is found but has no edges or is not found
		"""
		if target in self.graphs[self.to_predict]['TS']:
			return list(self.graphs[self.to_predict]['TS'][target])
		return []

	def normalisePredictorOverlap(self, filtered):
		"""Perform feature normalisation to range 0-1.
		The raw adjacencies for each feature are divided by the max value for that feature.
		Parameters
		----------
		filtered: dict 
			output of self.filterSparseAdjacency
		Returns
		-------
		all_normalised : dict
			Keys are edge types, values are dicts. The nexted dict is keyed by source node name and values are normalised
			adjacencies.
		all_overlap_max : dict
			Keys are edge types, values are the max adjacency in for that edge type.
		"""
		all_normalised = {}
		all_overlap_max = {}
		for network_name in filtered:
			norm = {}
			overlap_max = float(max(filtered[network_name]['overlap']))
			all_overlap_max[network_name] = overlap_max
			if overlap_max == 0:
				overlap_max = 1 #all values must = 0, divide by 1 so unchanged but still keyed by node name
			for i in range(len(filtered[network_name]['colnames'])):
				norm[filtered[network_name]['colnames'][i]] = filtered[network_name]['overlap'][i]/overlap_max
			all_normalised[network_name] = norm
		return all_normalised, all_overlap_max

	def weightPredictorOverlap(self, overlaps, weights):
		"""Multiply each feature by a weight.
		Parameters
		----------
		overlaps : dict 
			output of normalisePredictorOverlap
		weights : dict
			Keys are edge types, values are weights
		Returns
		-------
		weighted : dict
			dict with same structure as input overlaps, but with all values multiplied by their respective weights
		"""
		weighted = {}
		for network_name in overlaps:
			weighted[network_name] = {}
			for node_name in overlaps[network_name]:
				weighted[network_name][node_name] = overlaps[network_name][node_name] * weights[network_name]
		return weighted

	def score(self, overlaps):
		"""Calculate the final score for each source node from weighted features.
		Parameters
		----------
		overlaps : dict 
			output from self.weightPredictorOverlap, normalised and weighted features for each source node
		Returns
		-------
		scores : dict
			Keys are edge types, values are dicts keyed by source node name and values are scores.
		"""
		networks = overlaps.keys()
		nodes = overlaps[networks[0]].keys()
		scores = {}
		for n in nodes:
			scores[n] = 0
			for nw in networks:
				scores[n] += overlaps[nw][n]
		return scores

	def findOptimumThreshold(self, score, known, calculate_auc = False):
		"""Set the prediction threshold according to the objective function
		The objective function is set by self.objective_function
		Parameters
		----------
		score : dict 
			output from from self.score, keys are edge types, values are dicts keyed by source node name and values are scores.
		known : list
			source nodes with an edge to the target of type self.to_predict
		calculate_auc : bool
			Whether or not to calculare and return the AUC. Default True.
		Returns
		-------
		best : dict
			Contains many standard metrics for the model, e.g. F1 score, AUC, precision, recall, which have predictable names.
			Important proporties of the output are:
			'threshold' : cutoff value that maximises the objective function
			'unique_threshold' : bool, true if the same performance can be achieved with at least one different threshold
			'hits_known' : hits from the model that are already known in the input graph
			'hits_new' : hits from the model that are not already known in the input graph
			'is_hit' : bool list, hit status for every source node.
		"""
		node_names = score.keys()
		score = np.array([score[x] for x in score])
		known = np.in1d(node_names, known, assume_unique=True)
		thresholds = np.unique(score)
		placeholder = {}
		method = self.objective_function
		placeholder[method] = -1 #all current objective functions are in range 0-1 so the first result always replaces the placeholder
		best_performance = [placeholder]
		obj = Objective(score, known)
		if calculate_auc:
			x = [] #FPR = FP/(population N)
			y = [] #TPR TP/(population P)
			pop_pos = float(np.sum(known))
			pop_neg = float(len(known) - pop_pos)
		for t in thresholds:
			result = obj.evaluate(t)
			if result[method] > best_performance[0][method]:
				result['unique_threshold'] = True
				result['threshold'] = t
				best_performance = [result] #if there's a new best, reset to an array of one element
			elif result[method] == best_performance[0][method]: 
				#keep track of different thresholds that give equivalent results
				result['threshold'] = t
				result['unique_threshold'] = False
				best_performance.append(result)
				best_performance[0]['unique_threshold'] = False 
			#auc
			if calculate_auc:
				x.append(result['contingency']['fp']/ pop_neg)
				y.append(result['contingency']['tp'] / pop_pos)
		best = best_performance[0] #only return one result even if there are ties
		best['all_hits'] = set(itertools.compress(node_names, best['is_hit']))
		del best['is_hit']
		best['auc'] = "NA"
		if calculate_auc:
			x = np.array(x)
			y = np.array(y)
			best['auc'] = self.auc(x, y, True)
		return best

	def L2norm(self, weights):
		"""Regluarisation of weights
		Parameters
		----------
		weights : list
			Model parameters, weights of each feature.
		Returns
		-------
		Float : Float
			L2 regularisation of the weights
		"""
		return sum([x**2 for x in weights])

	def predict(self, target, calculate_auc = False):
		"""Train a predictive model for a given target.
		Optimum parameters are found using a grid search.
		Parameters
		----------
		target : str
			target node name to predict edges of type self.to_predict for
		calculate_auc: bool
			If True, the AUC is calculated and included in the output. Default False.
		
		Returns
		-------
		optimisation_result : dict
			Predictions from the trained model and various standard metrics such as precision, recall, F1, etc. 
			Output contains the model target and objective function so the results are self-describing. The most
			important proporties are:
			'all_hits' : all hit source nodes from the model
			'new_hits' : all hits from the model that are not known in the input graph
			'known_hits' : all hits from the model that are known in the input graph
			'weights' : dict of parameters in the trained model, keys are edge types
			'threshold' : threshold of trained model
		"""
		if self.to_predict == None or self.can_analyse == False:
			print "ERROR can't run prediction. self.to_predict = %s, self.can_analyse = %s" % (self.to_predict, self.can_analyse)
			email = self.email
			tempid = self.tempid
			error1 = "\"ERROR: Can't run prediction. self.to_predict = %s, self.can_analyse = %s.\"" % (self.to_predict, self.can_analyse) 
			error2 = "Please submit a new job."
			send_error(email, tempid, error1, error2)
			error_page(tempid, error1, error2)
			return render_template('error.html', tempid=tempid, error1=error1, error2=error2)
		known = self.getKnown(target)
		if len(known) == 0:
			print "ERROR no edges for target or target not in network: %s" % target
			email = self.email
			tempid = self.tempid
			error1 = "\"ERROR: No edges for target or target not in network: %s.\""  % target
			error2 = "This is probably because the number of input genes is too small to build a model. Please try to select more genes and/or phenotypes."
			send_error(email, tempid, error1, error2)
			error_page(tempid, error1, error2)
			return render_template('error.html', tempid=tempid, error1=error1, error2=error2)
		known_set = set(known)
		n_known = len(known)
		n_other = self.n_source_nodes - n_known
		grouped = self.groupSparseAdjacency(target)
		enrichment_pvals = self.enrichment(grouped, n_known, n_other)
		filtered = self.filterSparseAdjacency(enrichment_pvals, target)
		normalised, overlap_max = self.normalisePredictorOverlap(filtered)

		if self.require_all_predictors:
			no_predictor_overlap = [x for x in overlap_max if overlap_max[x] == 0]
			if len(no_predictor_overlap) > 0:
				print "self.require_all_predictors is %s and %s/%s networks have 0 predictor overlap" % (self.require_all_predictors, len(no_predictor_overlap), len(normalised))
				print "not optimising a model for %s." % (target)
				email = self.email
				tempid = self.tempid
				error1 = "\"Only build a model if predictors are found in all input graphs. %s/%s networks have 0 predictor overlap. Not optimising a model for %s.\"" % (len(no_predictor_overlap), len(normalised), target)
				error2 = "Please try to select more or different databases and/or to input more genes or diseases."
				send_error(email, tempid, error1, error2)
				error_page(tempid, error1, error2)
				optimisation_result = {}
				optimisation_result['model_target'] = target
				optimisation_result['model_built'] = False
				return render_template('error.html', tempid=tempid, error1=error1, error2=error2)

		if self.optimisation_method == "graph":
			weights_generator = self.createWeightsGenerator(self.min_weight, self.max_weight, self.step)
			if self.network_order == None:
				network_names = normalised.keys()
			else:
				network_names = self.network_order
			optimisation_result = {}
			optimisation_result[self.objective_function] = -1 #all current objectives are in range 0-1 so the first result always replaces this
			for weights in weights_generator:
				weights = dict(zip(network_names, weights))
				weighted = self.weightPredictorOverlap(normalised, weights)
				scores = self.score(weighted)
				best_threshold_for_weights = self.findOptimumThreshold(scores, known, calculate_auc)
				if best_threshold_for_weights[self.objective_function] > optimisation_result[self.objective_function]:
					optimisation_result = best_threshold_for_weights
					optimisation_result['weights'] = weights
					optimisation_result['count_equivalent_weights'] = 1
				elif best_threshold_for_weights[self.objective_function] == optimisation_result[self.objective_function]:
					optimisation_result['count_equivalent_weights'] += 1 #equivalent in terms of objective function score, not necessarily predictions made
					if self.ties == "minL2norm":
						lnorm_best = self.L2norm(optimisation_result['weights'].values())
						lnorm_now = self.L2norm(weights.values())
						if lnorm_now < lnorm_best:
							optimisation_result = best_threshold_for_weights
							optimisation_result['weights'] = weights
							optimisation_result['count_equivalent_weights'] = 1
			optimisation_result['known_hits'] = known_set.intersection(optimisation_result['all_hits'])
			optimisation_result['new_hits'] = optimisation_result['all_hits'].difference(known_set)
			optimisation_result['all_hits'] = list(optimisation_result['all_hits'])
			optimisation_result['new_hits'] = list(optimisation_result['new_hits'])
			optimisation_result['known_hits'] = list(optimisation_result['known_hits'])
		
		elif self.optimisation_method == "graph_sparse":
			#first search at half the density
			#print "coarse search"
			weights_generator = self.createWeightsGenerator(self.min_weight, self.max_weight, self.step * 2)
			optimisation_result = self.evaluate_weights(weights_generator, normalised, known, known_set, calculate_auc)
			#starting from the current best, test each weight +/- fine grain step
			start_weights = optimisation_result['weights']
			fine_weights = []
			
			if self.network_order == None:
				network_names = normalised.keys()
			else:
				network_names = self.network_order
			
			for x in network_names:
				mid = start_weights[x]
				out = [mid]
				low = mid - self.step
				high = mid + self.step
				if low >= self.min_weight:
					out.append(low)
				if high <= self.max_weight:
					out.append(high)
				out.sort()
				fine_weights.append(out)
			#print fine_weights
			weights_generator = itertools.product(*fine_weights)
			#print "fine search"
			optimisation_result = self.evaluate_weights(weights_generator, normalised, known, known_set, calculate_auc)
			
			
		else:
			print "No method definied to handle the optimisation method %s" % self.optimisation_method
			raise NameError(self.optimisation_method)
		
		optimisation_result['optimisation_method'] = self.optimisation_method
		optimisation_result['objective'] = self.objective_function
		optimisation_result['model_target'] = target
		optimisation_result['model_built'] = True
		optimisation_result['model_edge_type'] = self.to_predict
		optimisation_result['predictors'] = {}
		for network_name in filtered:
			optimisation_result['predictors'][network_name] = list(filtered[network_name]['predictors'])

		return optimisation_result
	
	def evaluate_weights(self, weights_generator, normalised, known, known_set, calculate_auc):
		if self.network_order == None:
				network_names = normalised.keys()
		else:
			network_names = self.network_order
		optimisation_result = {}
		optimisation_result[self.objective_function] = -1 #all current objectives are in range 0-1 so the first result always replaces this
		for weights in weights_generator:
			weights = dict(zip(network_names, weights))
			#print weights
			weighted = self.weightPredictorOverlap(normalised, weights)
			scores = self.score(weighted)
			best_threshold_for_weights = self.findOptimumThreshold(scores, known, calculate_auc)
			if best_threshold_for_weights[self.objective_function] > optimisation_result[self.objective_function]:
				optimisation_result = best_threshold_for_weights
				optimisation_result['weights'] = weights
				optimisation_result['count_equivalent_weights'] = 1
			elif best_threshold_for_weights[self.objective_function] == optimisation_result[self.objective_function]:
				optimisation_result['count_equivalent_weights'] += 1 #equivalent in terms of objective function score, not necessarily predictions made
				if self.ties == "minL2norm":
					lnorm_best = self.L2norm(optimisation_result['weights'].values())
					lnorm_now = self.L2norm(weights.values())
					if lnorm_now < lnorm_best:
						optimisation_result = best_threshold_for_weights
						optimisation_result['weights'] = weights
						optimisation_result['count_equivalent_weights'] = 1
		optimisation_result['known_hits'] = known_set.intersection(optimisation_result['all_hits'])
		optimisation_result['new_hits'] = optimisation_result['all_hits'].difference(known_set)
		optimisation_result['all_hits'] = list(optimisation_result['all_hits'])
		optimisation_result['new_hits'] = list(optimisation_result['new_hits'])
		optimisation_result['known_hits'] = list(optimisation_result['known_hits'])
		return optimisation_result

	def predictAll(self, calculate_auc=False):
		"""Train predictive models for all target nodes.
		Train predictive model for all target nodes of edges with the type self.to_predict. Not all targets
		will necessarily results in models depending on whether any enriched features are identified, and
		on self.require_all_predictors. The results is the same as manually calling self.predict on each 
		target, this function is for convenience.
		Parameters
		----------
		calculate_auc : bool
			If true, the AUC is calculated and returned for each model. Default False.
		Returns
		-------
		all_results : dict
			Keys are model target node names, values are the output of self.predict()
		"""
		if self.to_predict == None or self.can_analyse == False:
			print "can't run prediction. self.to_predict = %s, self.can_analyse = %s" % (self.to_predict, self.can_analyse)
			email = self.email
			tempid = self.tempid
			error1 = "\"ERROR: Can't run prediction. self.to_predict = %s, self.can_analyse = %s.\"" % (self.to_predict, self.can_analyse) 
			error2 = "Please submit a new job."
			send_error(email, tempid, error1, error2)
			error_page(tempid, error1, error2)
			return render_template('error.html', tempid=tempid, error1=error1, error2=error2)
		all_results = {}
		all_targets = self.graphs[self.to_predict]['TS'].keys()
		n_targets = len(all_targets)
		n = 1
		for target in all_targets:
			print "%s (%s/%s)" % (target, n, n_targets)
			n += 1
			all_results[target] = self.predict(target, calculate_auc)
		return all_results

	def loo(self, target, calculate_auc = False):
		"""Leave-one-out cross validation
		In each iteration, a single edge from a source node to the target node is deleted. A predictive model
		is trained on this modified data to determine whether the model predicts the missing (deleted) edge.
		Parameters
		----------
		target : str
			Target node name to predict edges of type self.to_predict for
		calculate_auc : bool
			If true, the AUC is calculated and returned for each model. Default False.
		Returns
		-------
		loo_results : dict
			Keys are names of known source nodes in the graph. Values are the objective function performance and 
			whether the deleted edge was predicted.
		"""
		if self.to_predict == None or self.can_analyse == False:
			print "ERROR can't run prediction. self.to_predict = %s, self.can_analyse = %s" % (self.to_predict, self.can_analyse)
			email = self.email
			tempid = self.tempid
			error1 = "\"ERROR: Can't run prediction. self.to_predict = %s, self.can_analyse = %s.\"" % (self.to_predict, self.can_analyse) 
			error2 = "Please submit a new job."
			send_error(email, tempid, error1, error2)
			error_page(tempid, error1, error2)
			return render_template('error.html', tempid=tempid, error1=error1, error2=error2)
		known = self.getKnown(target)
		target_node_id = self.graphs[self.to_predict]['graph'].vs.select(name_eq=target)
		loo_results = {}
		for k in known:
			#find the edge from this known source node to the target and delete it
			source_node_id = self.graphs[self.to_predict]['graph'].vs.select(name_eq=k)
			edge_to_delete = self.graphs[self.to_predict]['graph'].es.select(_between=(source_node_id, target_node_id))
			self.graphs[self.to_predict]['graph'].delete_edges(edge_to_delete)
			#update the master adjacency matrix
			self.sparseAdjacency()
			#run the prediction
			res = self.predict(target, calculate_auc)
			loo_results[k] = {'target':target, 'left_out_name':k, 'model_built': res['model_built']}
			if res['model_built']:
				ignored_is_hit = k in res['all_hits']
				loo_results[k]['was_predicted'] = ignored_is_hit
				loo_results[k]['objective_performance'] = res[self.objective_function]				

			#put the edge back
			self.graphs[self.to_predict]['graph'].add_edges([(source_node_id[0], target_node_id[0])])

		#update the master adjacency matrix so make sure it contains all edges again
		self.sparseAdjacency()
		return loo_results

	def k_fold(self, target, k, calculate_auc = False):
		"""Modified k-fold cross validation.
		This is a modidication of a standard k-fold cross validation. In this implementation, edges are deleted from 
		the graph and a predictive model is then trained on this modified data. Therefore the test set is not entirely 
		held out during training, instead it is included as true negative examples. The ability of the trained model 
		to predict the deleted edges is determined in every fold. 
		Parameters
		----------
		target : str
			Target node name to predict edges of type self.to_predict for.
		k : int
			The number of folds.
		calculate_auc : bool
			If true, the AUC is calculated and returned for each model. Default False.
		Returns
		-------
		all_folds : list
			Each item in the list is a dict. The result is the output of self.predict with additional properties.
			'left_out_predicted' : which of the deleted edges was predicted
			'proportion_predicted' : proportion of all deleted edges that was predicted
		"""
		#generate folds
		known = self.getKnown(target)
		if self.randomise_folds:
			#get known source nodes into random order
			np.random.shuffle(known)

		#number of edges to delete per fold
		edges_per_fold = len(known)/k
		remainder = len(known) % k

		if edges_per_fold == 0:
			print "specified fold size %s is too large for %s with %s known sources" % (k, target, len(known))
			email = self.email
			tempid = self.tempid
			error1 = "Specified fold size %s is too large for %s with %s known sources" % (k, target, len(known))
			error2 = "Please submit a new job."
			send_error(email, tempid, error1, error2)
			error_page(tempid, error1, error2)
			return render_template('error.html', tempid=tempid, error1=error1, error2=error2)

		start = 0
		stop = edges_per_fold
		all_folds = []
		target_node_id = self.graphs[self.to_predict]['graph'].vs.select(name_eq=target)
		for fold in range(k):
			if fold < remainder:
				stop += 1

			#find and delete the edges between these source nodes and the target
			delete_this_fold = known[start:stop]
			deleted_source_ids = []
			for source_name in delete_this_fold:
				source_node_id = self.graphs[self.to_predict]['graph'].vs.select(name_eq=source_name)
				edge_to_delete = self.graphs[self.to_predict]['graph'].es.select(_between=(source_node_id, target_node_id))
				self.graphs[self.to_predict]['graph'].delete_edges(edge_to_delete)
				deleted_source_ids.append(source_node_id)
			#update the master adjacency matrix
			self.sparseAdjacency()
			#run the prediction
			res = self.predict(target, calculate_auc)
			fold_result = {'target':target, 'left_out': delete_this_fold, 'model_built': res['model_built']}
			fold_result['n_known_train'] = len(known) - len(delete_this_fold)
			fold_result['n_known_test'] = len(delete_this_fold)
			if res['model_built']:
				ignored_is_hit = []
				for source_name in delete_this_fold:
					ignored_is_hit.append(source_name in res['new_hits'])
				fold_result['left_out_predicted'] = ignored_is_hit
				fold_result['proportion_predicted'] = float(len([x for x in ignored_is_hit if x]))/len(ignored_is_hit)
				fold_result['objective_performance'] = res[self.objective_function]	
				fold_result['contingency'] = res['contingency']
			all_folds.append(fold_result)			

			#put the edges back
			for source_node_id in deleted_source_ids:
				self.graphs[self.to_predict]['graph'].add_edges([(source_node_id[0], target_node_id[0])])

			start = stop
			stop += edges_per_fold

		#update the master adjacency matrix so make sure it contains all edges again
		self.sparseAdjacency()
		return all_folds

	def auc(self, x, y, reorder=False):
		"""Calculate AUC
		Credit to scipy.metrics 
		
		Parameters
		----------
		x : list
		y : list
		reorder : bool
			reorder the data points according to the x axis and using y to break ties.
			Default False.
		Returns
		-------
		area : float
			The area under the curve
		"""
		direction = 1
		if reorder:
			
			order = np.lexsort((y, x))
			x, y = x[order], y[order]
		else:
			dx = np.diff(x)
			if np.any(dx < 0):
				if np.all(dx <= 0):
					direction = -1
				else:
					raise ValueError("Reordering is not turned on, and the x array is not increasing: %s" % x)

		area = direction * np.trapz(y, x)
		if isinstance(area, np.memmap):
			# Reductions such as .sum used internally in np.trapz do not return a
			# scalar by default for numpy.memmap instances contrary to
			# regular numpy.ndarray instances.
			area = area.dtype.type(area)
		return area

	def getScores(self,target,weights):
		"""Calculate the score for all source nodes for a given set of weights.

		Not used internally, but a convenient way to calculate the score distribution for an arbitrary set of weights
		to manually explore how the distribution varies with weight, or to visualise the score distributino with the 
		trained model weights.

		Parameters
		----------
		target : str
			Target node name to predict edges of type self.to_predict for.

		weights : dict
			Keys are edge types, values are weights

		Returns
		-------
		scores : dict
			Keys are edge types, values are dicts keyed by source node name and values are scores.
		"""
		if self.to_predict == None or self.can_analyse == False:
			print "ERROR can't run prediction. self.to_predict = %s, self.can_analyse = %s" % (self.to_predict, self.can_analyse)
			email = self.email
			tempid = self.tempid
			error1 = "\"ERROR: Can't run prediction. self.to_predict = %s, self.can_analyse = %s.\"" % (self.to_predict, self.can_analyse) 
			error2 = "Please submit a new job."
			send_error(email, tempid, error1, error2)
			error_page(tempid, error1, error2)
			return render_template('error.html', tempid=tempid, error1=error1, error2=error2)
		known = self.getKnown(target)
		n_known = len(known)
		n_other = self.n_source_nodes - n_known
		grouped = self.groupSparseAdjacency(target)
		enrichment_pvals = self.enrichment(grouped, n_known, n_other)
		filtered = self.filterSparseAdjacency(enrichment_pvals, target)
		normalised, overlap_max = self.normalisePredictorOverlap(filtered)

		if self.require_all_predictors:
			no_predictor_overlap = [x for x in overlap_max if overlap_max[x] == 0]
			if len(no_predictor_overlap) > 0:
				print "self.require_all_predictors is %s and %s/%s networks have 0 predictor overlap" % (self.require_all_predictors, len(no_predictor_overlap), len(normalised))
				print "not optimising a model for %s." % (target)
				email = self.email
				tempid = self.tempid
				error1 = "\"Only build a model if predictors are found in all input graphs. %s/%s networks have 0 predictor overlap. Not optimising a model for %s.\"" % (len(no_predictor_overlap), len(normalised), target)
				error2 = "Please try to select more or different databases and/or to input more genes or diseases."
				send_error(email, tempid, error1, error2)
				error_page(tempid, error1, error2)
				optimisation_result = {}
				optimisation_result['model_target'] = target
				optimisation_result['model_built'] = False
				return render_template('error.html', tempid=tempid, error1=error1, error2=error2)

		weighted = self.weightPredictorOverlap(normalised, weights)
		scores = self.score(weighted)
		return scores



@app.route('/')
def index():
	return render_template('index.html')

@app.route('/index')
def _index():
	return render_template('index.html')

@app.route('/tutorial')
def _manual():
	return render_template('tutorial.html')

@app.route('/help')
def _help():
	return render_template('help.html')

@app.route('/predict', methods=['POST', 'GET'])
def _prepare():
	if request.method == 'POST':
		tempid = str(uuid.uuid4().hex)
		input_dict = request.form.to_dict()
		f1 = open('/var/www/kgp/FlaskApp/data/log.txt','a')
		f1.write("\ntempid: %s, result:" % tempid)
		f1.write(str(input_dict))

		email = str(input_dict['email'])
		file_name = '/var/www/kgp/FlaskApp/cache/output_file%s.csv' % tempid
		file_list = request.form.getlist('filename')
		uploaded_files = request.files.getlist("file[]")

		for file in uploaded_files:
			if file and allowed_file(file.filename):
				filename = secure_filename(file.filename)
				file.save(os.path.join(app.config['UPLOAD_FOLDER'],filename))
				file_list.append(os.path.join(app.config['UPLOAD_FOLDER'],filename))
		file_name_list = []
		for f in file_list:
			f = f.strip().split('/')
			file_name_list.append(f[-1])

		with open(file_name,'w') as wfd:
			judge = 0
			for f in file_list:
				if os.path.exists(f):
					with open(f,'r') as fd:
						judge = 1
						shutil.copyfileobj(fd, wfd)
			if judge == 0:
				with open("/var/www/kgp/FlaskApp/data/DisGeNet v7.0.csv",'r') as fd:
					shutil.copyfileobj(fd, wfd)

		phenotype = input_dict["phenotype"]
		if phenotype == "pheno": 
			input_name = input_dict['tname'].strip().split(',')
			gene_list = None
		elif phenotype == "defined":
			gname = input_dict['g2name'].strip().split(',')
			gene_list = [str(x).strip(' ') for x in gname]
			input_name = input_dict['t2name'].strip().split(',')
		elif phenotype == "from_genes":
			gname = input_dict['g3name'].strip().split(',')
			gene_list = [str(x).strip(' ') for x in gname]
			input_name = None
		if email != '':
			send_start(email, tempid)

		try:
			keywords = input_dict["keywords"]
		except:
			keywords = 'Flase'
		thr = Thread(target=predict, args=(tempid, email, file_name_list, file_name, gene_list, phenotype, input_name,keywords,))
		thr.start()
		return redirect(url_for('result',tempid=tempid))

def send_start(email, tempid):
	yag = yagmail.SMTP()
	subject = "DGLinker: Your project has been submitted."
	contents = "Dear user,\n\n Your project has been submited. You can find your result at the following link when the job will be completed: \n https://dglinker.rosalind.kcl.ac.uk/result/%s \n\n Best wishes, \n The DGLinker Team" % tempid
	yag.send(email, subject, contents)


def predict(tempid, email, file_name_list, file_name, gene_list, phenotype, input_name,keywords):
	if phenotype == "pheno":
		voca = pd.read_csv('/var/www/kgp/FlaskApp/data/vocabulary.csv',header=0,index_col=0)
		vocadict = open("/var/www/kgp/FlaskApp/data/vocabulary.csv", 'r')
		matched_diseidlist = []
		phenotype_used = []
		for i in input_name:
			if i != '':
				if keywords == "True":
					for row in vocadict:
						row = row.strip().split(',')
						if str(i).strip().lower() in row[1].strip().lower():
							matched_diseidlist.append(row[0])
							phenotype_used.append(row[1])
				else:
					matched_diseid = voca.index[voca['diseaseName'] == str(i)].tolist()
					matched_diseidlist.append(matched_diseid[0])
		f2 = open(file_name,'r') #dataset,5 columns
		merge_filename = "/var/www/kgp/FlaskApp/cache/merged_file%s.csv" % tempid
		with open(merge_filename,'w') as new:
			new.write('Source node name,Source node type,Relationship type,Target node type,Target node name\n')
			for line in f2:
				line_list = line.strip().split(',')
				if line_list[4] in matched_diseidlist:
					new.write('%s,%s,%s,%s,%s\n' % (line_list[0],line_list[1],line_list[2],line_list[3],matched_diseidlist[0]))
				else:
					new.write(line)
	elif phenotype == "defined":
		voca = pd.read_csv('/var/www/kgp/FlaskApp/data/vocabulary.csv',header=0,index_col=0)
		vocadict = open("/var/www/kgp/FlaskApp/data/vocabulary.csv", 'r')
		matched_diseidlist = []
		phenotype_used = []
		f1 = open('/var/www/kgp/FlaskApp/data/log.txt','a')
		for i in input_name:
			if i != '':
				if keywords == "True":
					for row in vocadict:
						row = row.strip().split(',')
						if str(i).strip().lower() in row[1].lower():
							matched_diseidlist.append(row[0])
							phenotype_used.append(row[1])
							f1.write(row[0])
							f1.write(row[1])
				else:
					matched_diseid = voca.index[voca['diseaseName'] == str(i)].tolist()
					matched_diseidlist.append(matched_diseid[0])
		f2 = open(file_name,'r') #dataset,5 columns
		merge_filename = "/var/www/kgp/FlaskApp/cache/merged_file%s.csv" % tempid
		f1 = open('/var/www/kgp/FlaskApp/data/log.txt','a')
		with open(merge_filename,'w') as new:
			new.write('Source node name,Source node type,Relationship type,Target node type,Target node name\n')
			for i in gene_list:
				if str(i) != "":
					new.write('%s,%s,%s,%s,%s\n' % (str(i),'gene','risk_gene','disease',str(matched_diseidlist[0])))
			for line in f2:
				line_list = line.strip().split(',')
				if line_list[4] not in matched_diseidlist:
					new.write(line)
	elif phenotype == "from_genes":
		f2 = open(file_name,'r') #dataset,5 columns
		merge_filename = "/var/www/kgp/FlaskApp/cache/merged_file%s.csv" % tempid
		with open(merge_filename,'w') as new:
			new.write('Source node name,Source node type,Relationship type,Target node type,Target node name\n')
			for i in gene_list:
				new.write('%s,%s,%s,%s,%s\n' % (str(i),'gene','risk_gene','disease','USER_DEFINED_DISEASE'))
			for line in f2:
				new.write(line)

	data = pd.read_csv(merge_filename, header=0, error_bad_lines=False, low_memory=False)
	newDf = data.drop_duplicates()
	newDf.to_csv("/var/www/kgp/FlaskApp/cache/new_merged_file%s.csv" % tempid, index=False)
	merge_filename = "/var/www/kgp/FlaskApp/cache/new_merged_file%s.csv" % tempid

	order_dict = pd.read_csv(merge_filename,header=0,index_col=0)
	order_list = []
	for i in order_dict['Relationship type']:
		if i not in order_list:
			order_list.append(str(i))
	if phenotype == "pheno" or phenotype == 'defined':
		target_name = matched_diseidlist[0]  # input_dict['tname']  # "C0027849"
	else:
		target_name = 'USER_DEFINED_DISEASE'
	
	ep = EdgePrediction(email=email, tempid=tempid)
	ep.CSV_to_graph(fname = merge_filename)
	ep.preprocess()
	ep.to_predict = 'risk_gene'
	# ep.network_order = order_list  # input_dict['order']  # ['HAS_SIDE_EFFECT', 'DRUG_TARGETS', 'INDICATED_FOR']
	result = ep.predict(target=target_name, calculate_auc=True)
	new_predictions = result['new_hits']
	known_predictions = result['known_hits']
	weights = result['weights']
	scores = ep.getScores(target_name, weights)
	# for i in range(len(new_predictions)):
		# ranking.append([scores[i][0],scores[i][1]])
		# ranking.append(scores[i][0])

	result_file = open('/var/www/kgp/FlaskApp/cache/result_file%s.txt' % tempid, 'a')
	job_summary = open('/var/www/kgp/FlaskApp/model_results/job_summary_%s.txt' % tempid, 'a')
	download_results = open('/var/www/kgp/FlaskApp/model_results/results_%s.txt' % tempid, 'a')
	result_file.write("%s\t%s\n" % ('Project ID',str(tempid)))
	job_summary.write("%s\t%s\n" % ('Project ID',str(tempid)))
	if email != '':
		result_file.write("%s\t%s\n" % ('Result Sent to',email))
		job_summary.write("%s\t%s\n" % ('Project ID',str(tempid)))
	result_file.write("%s\t%s\n" % ("Dataset Used", str([i.encode('utf8') for i in file_name_list]).replace('.csv','')))
	job_summary.write("%s\t%s\n" % ("Dataset Used", ', '.join([i.encode('utf8') for i in file_name_list]).replace('.csv','')))
	if keywords == 'True':
		result_file.write("%s\t%s\n" % ("Keywords Used", "True"))
		job_summary.write("%s\t%s\n" % ("Keywords Used", "True"))
	else:
		result_file.write("%s\t%s\n" % ("Keywords Used", "False"))
		job_summary.write("%s\t%s\n" % ("Keywords Used", "False"))
	if phenotype == "pheno":
		if keywords == 'True':
			result_file.write("%s\t%s\n" % ('Predict Mode', 'Predict from Phenotype(s)'))
			result_file.write("%s\t%s\n" % ('Phenotype(s) Used', str([i.encode('utf8') for i in phenotype_used])))
			job_summary.write("%s\t%s\n" % ('Predict Mode', 'Predict from Phenotype(s)'))
			job_summary.write("%s\t%s\n" % ('Phenotype(s) Used', ', '.join(([i.encode('utf8') for i in phenotype_used]))))
		else:
			result_file.write("%s\t%s\n" % ('Predict Mode', 'Predict from Phenotype(s)'))
			result_file.write("%s\t%s\n" % ('Phenotype(s) Used', str([i.encode('utf8') for i in input_name])))
			job_summary.write("%s\t%s\n" % ('Predict Mode', 'Predict from Phenotype(s)'))
			job_summary.write("%s\t%s\n" % ('Phenotype(s) Used', ', '.join(([i.encode('utf8') for i in input_name]))))
		result_file.write("%s\t%s\n" % ('Known Genes',result['known_hits']))
	elif phenotype == "defined":
		if keywords == 'True':
			result_file.write("%s\t%s\n" % ('Predict Mode', "Predict from User-defined Disease(s)"))
			result_file.write("%s\t%s\n" % ('Disease(s) Defined', str([i.encode('utf8') for i in phenotype_used])))
			result_file.write("%s\t%s\n" % ('By Genes', str(gene_list)))
			job_summary.write("%s\t%s\n" % ('Predict Mode', "Predict from User-defined Disease(s)"))
			job_summary.write("%s\t%s\n" % ('Disease(s) Defined', ', '.join(([i.encode('utf8') for i in phenotype_used]))))
			job_summary.write("%s\t%s\n" % ('By Genes', str(gene_list)))
		else:
			result_file.write("%s\t%s\n" % ('Predict Mode', "Predict from User-defined Disease(s)"))
			result_file.write("%s\t%s\n" % ('Disease(s) Defined', str(input_name)))
			result_file.write("%s\t%s\n" % ('By Genes', str(gene_list)))
			job_summary.write("%s\t%s\n" % ('Predict Mode', "Predict from User-defined Disease(s)"))
			job_summary.write("%s\t%s\n" % ('Disease(s) Defined', str(input_name)))
			job_summary.write("%s\t%s\n" % ('By Genes', str(gene_list)))
	elif phenotype == "from_genes":
		result_file.write("%s\t%s\n" % ('Predict Mode', "Predict from Gene(s)"))
		result_file.write("%s\t%s\n" % ('Genes Used', str(gene_list)))
		job_summary.write("%s\t%s\n" % ('Predict Mode', "Predict from Gene(s)"))
		job_summary.write("%s\t%s\n" % ('Genes Used', str(gene_list)))

	result_file.write("%s\t%s\n" % ('Total Number of Genes','%d' %  result['hits_total']))
	result_file.write("%s\t%s\n" % ('Number of Predicted Genes','%d' %  result['hits_new']))
	hits_known = result['hits_total'] - result['hits_new']
	result_file.write("%s\t%s\n" % ('Number of Known Genes', '%d' % hits_known))
	job_summary.write("%s\t%s\n" % ('Known Genes',', '.join(result['known_hits'])))
	job_summary.write("%s\t%s\n" % ('Total Number of Genes','%d' %  result['hits_total']))
	job_summary.write("%s\t%s\n" % ('Number of Predicted Genes','%d' %  result['hits_new']))
	job_summary.write("%s\t%s\n" % ('Number of Known Genes', '%d' % hits_known))

	result_file.write("%s\n" % 'More Information About Genes')
	download_results.write("%s\n" % 'Rank\tAssociation type\tGene name\tScore')
	# rank for all
	ranking = sorted(scores.items(), key=lambda item: item[1], reverse=True)
	num = 0
	for tup in ranking:
		if tup[0] in known_predictions:
			result_file.write("%d\t%s\t%s\t%.3f\n" % (num+1,'Known',tup[0],tup[1]))
			download_results.write("%d\t%s\t%s\t%.3f\n" % (num+1,'Known',tup[0],tup[1]))
		elif tup[0] in new_predictions:
			result_file.write("%d\t%s\t%s\t%.3f\n" % (num+1,'Predicted',tup[0],tup[1]))
			download_results.write("%d\t%s\t%s\t%.3f\n" % (num+1,'Predicted',tup[0],tup[1]))
		num += 1
	result_file.flush()
	result_file.close()
	job_summary.flush()
	job_summary.close()
	download_results.flush()
	download_results.close()

	# make directory
	parent_dir = '/var/www/kgp/FlaskApp/static/'
	directory = 'enrichr%s' % tempid
	path = os.path.join(parent_dir, directory)
	if not os.path.exists(path):
		os.mkdir(path)

	raw_file = open('/var/www/kgp/FlaskApp/cache/output_file%s.csv' % tempid, 'r')

	top_100_new = new_predictions[0:100]
	# graph file for top_50_all
	top_50_all = [tup[0] for tup in ranking][0:50]
	line_duplicate = []
	with open('/var/www/kgp/FlaskApp/cache/raw_graph_file%s.csv' % tempid, 'w') as raw_graph_file:
		for line in raw_file:
			linelist = line.strip().split(',')
			if linelist[0] in top_50_all and line not in line_duplicate:
				raw_graph_file.write(line)

	with open('/var/www/kgp/FlaskApp/cache/raw_graph_file%s.csv' % tempid, 'r') as reader:
		d = OrderedDict([('Source node name', []), ('Source node type', []), ('Target node type', []), ('Target node name', [])])
		for row in reader:
			row = row.strip().split(',')
			d['Source node name'].append(row[0])
			d['Source node type'].append(row[1])
			d['Target node type'].append(row[3])
			d['Target node name'].append(row[4])

	graph_file = open('/var/www/kgp/FlaskApp/cache/graph_file%s.txt' % tempid, 'a')
	download_graph_file = open('/var/www/kgp/FlaskApp/model_results/graph_file%s.txt' % tempid, 'a')
	download_graph_file.write('Gene\tOthers\tGene_class\tOthers_Class\tTotal number of interactions gene\tTotal number of interactions others\n')

	with open("/var/www/kgp/FlaskApp/data/vocabulary.txt", 'r') as f:
		vocadict = {}
		dup = []
		for row in f:
			if row not in dup:
				dup.append(row)
				row = row.strip().split('\t')
				vocadict[row[0]] = row[1]

	genecount = Counter(d['Source node name'])
	othercount = Counter(d['Target node name'])
	for i in range(len(d['Source node name'])):
		if d['Source node name'][i] in new_predictions:
			genetype = 'Predicted'
		else:
			genetype = 'Known'
		if d['Target node type'][i] == 'disease' and othercount[d['Target node name'][i]] > 1:
			graph_file.write('%s\t%s\t%s\t%s\t%d\t%d\n' % (d['Source node name'][i],vocadict[d['Target node name'][i]],genetype,d['Target node type'][i],genecount[d['Source node name'][i]],othercount[d['Target node name'][i]]))
			download_graph_file.write('%s\t%s\t%s\t%s\t%d\t%d\n' % (d['Source node name'][i],vocadict[d['Target node name'][i]],genetype,d['Target node type'][i],genecount[d['Source node name'][i]],othercount[d['Target node name'][i]]))
		elif d['Target node type'][i] != 'disease' and othercount[d['Target node name'][i]] > 1:
			graph_file.write('%s\t%s\t%s\t%s\t%d\t%d\n' % (d['Source node name'][i],str(d['Target node name'][i]),genetype,d['Target node type'][i],genecount[d['Source node name'][i]],othercount[d['Target node name'][i]]))
			download_graph_file.write('%s\t%s\t%s\t%s\t%d\t%d\n' % (d['Source node name'][i],str(d['Target node name'][i]),genetype,d['Target node type'][i],genecount[d['Source node name'][i]],othercount[d['Target node name'][i]]))
	graph_file.close()
	download_graph_file.close()

	if not os.path.exists('/var/www/kgp/FlaskApp/static/graph_dir%s' % tempid):
		os.mkdir('/var/www/kgp/FlaskApp/static/graph_dir%s' % tempid)
	
	for gene in top_50_all:
		with open('/var/www/kgp/FlaskApp/static/graph_dir%s/%s.txt' % (tempid,gene), 'a') as each:   
			each.write('Gene\tOthers\tGene_class\tOthers_Class\tTotal number of interactions gene\tTotal number of interactions others\n')

	graph_dict = {}

	with open('/var/www/kgp/FlaskApp/cache/graph_file%s.txt' % tempid, 'r') as reader:
		allrow = []
		for row in reader:
			rowline = row.strip().split('\t')
			if rowline not in allrow:
				allrow.append(rowline)
				graph_dict.setdefault(str(rowline[0]),[]).append(str(rowline[1]))
				with open('/var/www/kgp/FlaskApp/static/graph_dir%s/%s.txt' % (tempid,rowline[0]), 'a') as each:
					each.write(row)
		
	for gene in graph_dict:
		alldata = []
		for others in graph_dict[gene]:
			for i in allrow:
				if others == i[1] and gene !=i[0] and i not in alldata:
					with open('/var/www/kgp/FlaskApp/static/graph_dir%s/%s.txt' % (tempid,gene), 'a') as each:
						each.write('\t'.join(i))
						each.write('\n')
						alldata.append(i)
	
	folder = '/var/www/kgp/FlaskApp/static/graph_dir%s' % tempid
	
	threa = Thread(target=net_visualization, args=(folder,))
	threa.start()

	# analyze_gene_list: top 100 new
	userID = []
	ENRICHR_URL = 'http://maayanlab.cloud/Enrichr/addList'
	genes_str = '\n'.join(top_100_new)
	description = 'top 100 new Genes'
	payload = {
		'list': (None, genes_str),
		'description': (None, description)
	}
	response = requests.post(ENRICHR_URL, files=payload, timeout=80000)
	if not response.ok:
		raise Exception('Error analyzing gene list')
	iddata = json.loads(response.text) # dict()
	userID.append(iddata['userListId'])

	# analyze_gene_list: known
	ENRICHR_URL = 'http://maayanlab.cloud/Enrichr/addList'
	genes_str = '\n'.join(known_predictions)
	description = 'known Genes'
	payload = {
		'list': (None, genes_str),
		'description': (None, description)
	}
	response = requests.post(ENRICHR_URL, files=payload, timeout=80000)
	if not response.ok:
		raise Exception('Error analyzing gene list')
	iddata = json.loads(response.text) # dict()
	userID.append(iddata['userListId'])

	# analyze_gene_list: top 100 new
	ENRICHR_URL = 'http://maayanlab.cloud/Enrichr/addList'
	gene_list = top_100_new + known_predictions
	genes_str = '\n'.join(gene_list)
	description = 'top 100 new Genes + known Genes'
	payload = {
		'list': (None, genes_str),
		'description': (None, description)
	}
	response = requests.post(ENRICHR_URL, files=payload, timeout=80000)
	if not response.ok:
		raise Exception('Error analyzing gene list')
	iddata = json.loads(response.text) # dict()
	userID.append(iddata['userListId'])

	# get_enrichment_results
	gene_set_list = ['ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X', 'KEGG_2019_Human', 'GO_Biological_Process_2018', 'GO_Cellular_Component_2018', 'GO_Molecular_Function_2018', 'OMIM_Disease', 'GWAS_Catalog_2019','DSigDB', 'Human_Gene_Atlas']
	namelist = ['New', 'Known', 'Union']
	for library in gene_set_list:
		i = -1   
		for ID in userID:
			i += 1
			try:
				ENRICHR_URL = 'http://maayanlab.cloud/Enrichr/export'
				query_string = '?userListId=%s&filename=%s&backgroundType=%s'
				user_list_id = ID
				filename = '/var/www/kgp/FlaskApp/static/enrichr%s/%s_%s_%s.txt' % (tempid,tempid,namelist[i],library)
				gene_set_library = library # str(library)
				url = ENRICHR_URL + query_string % (user_list_id, filename, gene_set_library)
				response = requests.get(url, stream=True)
				with open(filename, 'wb') as f:
					for chunk in response.iter_content(chunk_size=1024): 
						if chunk:
							f.write(chunk)
			except Exception as e:
				pass
			# thread = Thread(target=get_enrichment_results, args=(tempid,ID,library,))
			# thread.start()
			# time.sleep(5)
	# time.sleep(30)	 
	f = open ('/var/www/kgp/FlaskApp/static/enrichr%s/apidone.txt' % tempid, 'w')
	for ID in userID:
		f.write('%s,' % ID)
	f.close()

	if not os.path.exists('/var/www/kgp/FlaskApp/static/enrichr%s/graphdone.txt' % tempid):
		f = open ('/var/www/kgp/FlaskApp/static/enrichr%s/apidone.txt' % tempid, 'r')
		userID = []
		for line in f:
			line = line.strip().split(',')
			for i in line:
				if i != '':
					userID.append(i)
		gene_set_list = ['ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X', 'KEGG_2019_Human', 'GO_Biological_Process_2018', 'GO_Cellular_Component_2018', 'GO_Molecular_Function_2018', 'OMIM_Disease', 'GWAS_Catalog_2019','DSigDB', 'Human_Gene_Atlas']
		parent_dir = '/var/www/kgp/FlaskApp/static/'
		directory = 'enrichr%s' % tempid
		path = os.path.join(parent_dir, directory)

		i = -1
		idlist = []
		for ID in userID:
			namelist = ['New', 'Known', 'Union']
			if ID not in idlist:
				i += 1
				idlist.append(ID)
			for library in gene_set_list:
				graphname = '/var/www/kgp/FlaskApp/static/enrichr%s/%s_%s_%s.png' % (tempid,tempid,namelist[i],library)
				filename = '/var/www/kgp/FlaskApp/static/enrichr%s/%s_%s_%s.txt' % (tempid,tempid,namelist[i],library)
				enrichr_grpah(tempid,ID,library,graphname,filename)
		
		f = open ('/var/www/kgp/FlaskApp/static/enrichr%s/graphdone.txt' % tempid, 'w')
		f.write('This file is used to claim that all of the enrichment and corresponding visulization are done.')
		f.close()

	send_result(email, tempid)


def send_result(email, tempid):
	if email != '':
		yag = yagmail.SMTP()
		subject = "Your result from DGLinker."
		contents = "Dear user,\n\n Please find attached a summary of your results. You can download the complete results at the following link: \n https://dglinker.rosalind.kcl.ac.uk/result/%s \n\n Best wishes, \n The DGLinker Team" % tempid
		attlist = ['/var/www/kgp/FlaskApp/cache/result_file%s.txt' % tempid]
		yag.send(email, subject, contents, attlist)


# @user.route('/result', defaults={'tempid': None})
@app.route('/result/<string:tempid>')
def result(tempid):
	if os.path.exists('/var/www/kgp/FlaskApp/cache/error_file%s.txt' % tempid):
		f1 = open('/var/www/kgp/FlaskApp/cache/error_file%s.txt' % tempid, 'r')
		errors = []
		for line in f1:
			errors.append(line)
		error1 = errors[0]
		error2 = errors[1]
		return render_template('error.html', tempid=tempid, error1=error1, error2=error2)
	if os.path.exists('/var/www/kgp/FlaskApp/static/enrichr%s/graphdone.txt' % tempid):
		final_file = open('/var/www/kgp/FlaskApp/cache/result_file%s.txt' % tempid,'r')
		f1 = open('/var/www/kgp/FlaskApp/data/log.txt','a')
		result = OrderedDict()
		all_hits = []
		new_hits = []
		known_hits = []
		all_hits_list = []
		new_hits_list = []
		known_hits_list = []
		judge = 0
		count_for_new = 0
		count_for_known = 0
		count_for_all = 0
		for line in final_file:
			line = line.strip().split('\t')
			if line[0] != 'More Information About Genes' and judge==0:
				if '[' in line[1]:
					new_line = line[1].replace('"','')
					new_line = new_line.replace('[','')
					new_line = new_line.replace(']','')
					new_line = new_line.replace('\'','')
					result[line[0]] = new_line
				else:
					result[line[0]] = line[1]
			elif line[0] == 'More Information About Genes':
				judge = 1
			elif line[0] != 'More Information About Genes' and judge==1:
				count_for_all += 1
				all_hits.append([count_for_all]+line[1:])
				all_hits_list.append(line[2])
				if line[1] == 'Predicted':
					count_for_new += 1
					new_hits.append([count_for_new]+line[1:])
					new_hits_list.append(line[2])
				else:
					count_for_known += 1
					known_hits.append([count_for_known]+line[1:])
					known_hits_list.append(line[2])
		return render_template('result.html',result=result, all_hits=all_hits, tempid=tempid)
	else:
		return render_template('redirect.html', tempid=tempid)

def net_visualization(folder):
	file_list = os.listdir(folder)
	for gene_file in file_list :
		# set the physics layout of the network
		reference = gene_file.split(".")[0]
		got_net = Network(height="75%", width="75%", font_color="black", heading=reference)
		got_net.barnes_hut()
		
		got_data = pd.read_csv(folder+"/"+gene_file,sep='\t')

		sources = got_data['Gene']
		targets = got_data['Others']
		gene_class = got_data['Gene_class']
		others_Class = got_data['Others_Class']
		weight_gene = got_data['Total number of interactions gene']
		weight_other = got_data['Total number of interactions others']

		edge_data = zip(sources, targets, weight_gene, weight_other, gene_class, others_Class)

		for e in edge_data:
			src = e[0]
			dst = e[1]
			weight_1 = e[2]
			weight_2 = e[3]
			gene_class = e[4]
			others_class = e[5]
			if gene_class == "Known" :
				if src == reference :
					got_net.add_node(src, src, title=src, group=1, shape = "dot", size = int(weight_1)*4+1000 , physics = True , borderWidth = 3 , borderWidthSelected = 7)
				else:
					got_net.add_node(src, src, title=src, group=2, shape = "dot", size = int(weight_1)*4 , physics = True , borderWidth = 3 , borderWidthSelected = 7)
			if gene_class == "Predicted" :
				if src == reference :
					got_net.add_node(src, src, title=src, group=1, shape = "dot", size = int(weight_1)*4+1000, physics = True , borderWidth = 3 , borderWidthSelected = 7)
				else:
					got_net.add_node(src, src, title=src, group=3, shape = "dot", size = int(weight_1)*4, physics = True , borderWidth = 3 , borderWidthSelected = 7)
			if others_class == "ANOTHERGEN" :
				got_net.add_node(dst, dst, title=dst, group=4, shape = "box", size = int(weight_2)+50, physics = True , level = 2)
			if others_class == "GOterm" :
				got_net.add_node(dst, dst, title=dst, group=5, shape = "box", size = int(weight_2)+50, physics = True , level = 2)
			if others_class == "disease" :
				got_net.add_node(dst, dst, title=dst, group=6, shape = "box", size = int(weight_2)+50, physics = True , level = 2)
			if others_class == "pubmedID" :
				got_net.add_node(dst, dst, title=dst, group=7, shape = "box", size = int(weight_2)+50, physics = True , level = 2)
			if others_class == "expression" :
				got_net.add_node(dst, dst, title=dst, group=8, shape = "box", size = int(weight_2)+50, physics = True , level = 2)
			
			got_net.add_edge(src, dst, physics = True)

			neighbor_map = got_net.get_adj_list()

	# add neighbor data to node hover data
		for node in got_net.nodes:
			node["title"] += "<br>Neighbors:<br>" + "<br>".join(neighbor_map[node["id"]])
			node["value"] = len(neighbor_map[node["id"]])

		got_net.show(folder+'/graph_'+gene_file.split(".")[0]+".html")
		with open(folder+'/graph_'+gene_file.split(".")[0]+".html", 'r') as origin:
			modified = []
			for line in origin:
				modified.append(line)
			flag = 0
		# os.remove(folder+'/'+gene_file.split(".")[0]+".html")
		with open(folder+'/'+gene_file.split(".")[0]+".html", 'a') as newer:
			for i in range(len(modified)-10):
				if modified[i] == '<script type="text/javascript">\n':
					flag = 1
				if flag == 1:
					newer.write(modified[i])
			
			if '		network.on("stabilizationProgress", function(params) {\n' not in modified:
				newer.write('''
						network.on("stabilizationProgress", function(params) {
						document.getElementById('loadingBar').removeAttribute("style");
						var maxWidth = 496;
						var minWidth = 20;
						var widthFactor = params.iterations/params.total;
						var width = Math.max(minWidth,maxWidth * widthFactor);

						document.getElementById('bar').style.width = width + 'px';
						document.getElementById('text').innerHTML = Math.round(widthFactor*100) + '%';
					});
					network.once("stabilizationIterationsDone", function() {
						document.getElementById('text').innerHTML = '100%';
						document.getElementById('bar').style.width = '496px';
						document.getElementById('loadingBar').style.opacity = 0;
						// really clean the dom element
						setTimeout(function () {document.getElementById('loadingBar').style.display = 'none';}, 500);
					});
					

						return network;

				}

				drawGraph();

			</script>
			''')
			else:
				newer.write('''
					return network;

				}

				drawGraph();

			</script>
			''')



def get_enrichment_results(tempid,ID,library): 
	try:
		ENRICHR_URL = 'http://maayanlab.cloud/Enrichr/export'
		query_string = '?userListId=%s&filename=%s&backgroundType=%s'
		user_list_id = ID
		filename = '/var/www/kgp/FlaskApp/static/enrichr%s/%s_%s.txt' % (tempid,ID,library)
		gene_set_library = library # str(library)
		if not os.path.exists(filename):
			url = ENRICHR_URL + query_string % (user_list_id, filename, gene_set_library)
			response = requests.get(url, stream=True)
			with open(filename, 'wb') as f:
				for chunk in response.iter_content(chunk_size=1024): 
					if chunk:
						f.write(chunk)
	except Exception as e:
		pass 
	
	
def enrichr_grpah(tempid,ID,library,graphname,filename):
	title = ['ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X','KEGG_2019_Human','GO_Biological_Process_2018','GO_Cellular_Component_2018','GO_Molecular_Function_2018','OMIM_Disease','GWAS_Catalog_2019','DSigDB','Human_Gene_Atlas']
	if os.path.exists(filename):
		bar_data = pd.read_csv(filename,sep="\t",header=0)
		gene_name_list = list(bar_data["Genes"])
		pvalue = list(bar_data["Adjusted P-value"])
		rawp = list(bar_data["P-value"])
		odds_radio = list(bar_data["Odds Ratio"])
		overlap = list(bar_data["Overlap"])
		term = list(bar_data["Term"])
		zipped = zip(gene_name_list,pvalue,rawp,odds_radio,overlap,term)
		sort_zipped = sorted(zipped,key=lambda x:x[1])
		result = zip(*sort_zipped)
		gene_name_list,pvalue,rawp,odds_radio,overlap,term = [list(x)[0:10][::-1] for x in result]
		logp = 0 - np.log10(pvalue)
		fig, ax = plt.subplots()
		color_list = []
		for i in range(len(pvalue)):
			if float(pvalue[i]) < 0.05:
				color_list.append('#8dd3c7')
			elif float(pvalue[i]) > 0.05 and float(rawp[i]) < 0.05:
				color_list.append('#ffffb3')
			elif float(pvalue[i]) > 0.05 and float(rawp[i]) > 0.05:
				color_list.append('#EA9999')

		b = ax.barh(range(len(gene_name_list)), logp, color=color_list)
		i = 0
		for rect in b:
			w = rect.get_width()
			# ax.text(w, rect.get_y()+rect.get_height()/1.5, '%.1f(%s)' % (float(odds_radio[i]),str(overlap[i])), ha='left', va='center', fontsize=7)
			ax.text(0, rect.get_y()+rect.get_height()/2, ' %s' % str(term[i]), ha='left', va='center',fontsize=8)
			i += 1
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.set_yticks(range(len(gene_name_list)))
		ylabel_list = ['%.1f(%s)' % (float(odds_radio[i]),str(overlap[i])) for i in range(len(pvalue))]
		ax.set_yticklabels(ylabel_list)
		plt.xlabel("-log($\mathregular{P_{adj}}$)")
		for i in title:
			if library in i:
				plt.title(str(i))
		plt.tight_layout()
		fig.savefig(graphname, dpi=400)
		plt.clf()
		plt.close()



@app.route('/download/<tempid>')
def downloadFile(tempid):
	path = '/var/www/kgp/FlaskApp/cache/zip%s.zip' % tempid
	if not os.path.exists(path):
		filelist1 = ['/var/www/kgp/FlaskApp/model_results/job_summary_%s.txt' % tempid, '/var/www/kgp/FlaskApp/model_results/results_%s.txt' % tempid, '/var/www/kgp/FlaskApp/model_results/graph_file%s.txt' % tempid]
		dirpath1 = '/var/www/kgp/FlaskApp/static/enrichr%s/' % tempid
		dirpath2 = '/var/www/kgp/FlaskApp/static/graph_dir%s/' % tempid
		filelist2 = search(dirpath1, '_')
		filelist3 = search(dirpath2, '.txt')
		filelist4 = search(dirpath2, '_')

		for f in filelist1 + filelist2 + filelist3 + filelist4:
			addfile('/var/www/kgp/FlaskApp/cache/zip%s.zip' % tempid, f)
	return send_file(path, as_attachment=True)

def addfile(zipfilename, filename):
	dir = "/".join([i for i in filename.strip().split("/")][:-2])
	new_filename = "/".join([i for i in filename.strip().split("/")][-2:])
	os.chdir(dir)
	with zipfile.ZipFile(zipfilename, 'a') as z:
		z.write(new_filename)

def search(path, word):
	filelist = []
	for filename in os.listdir(path):
		fp = os.path.join(path, filename)
		if os.path.isfile(fp) and word in filename:
			filelist.append(fp)
	return filelist


@app.route('/find', methods=['POST', 'GET'])
def _find():
	if request.method == 'POST':
		input_dict = request.form.to_dict()
		tempid = input_dict['jobname']
		return redirect(url_for('result',tempid=tempid))


@app.route('/downloads')
def downloads():
	return render_template('downloads.html')


if __name__ == '__main__':
	app.static_folder = 'static'
	app.run(debug=True)


