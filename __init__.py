
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
import EdgePrediction
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


UPLOAD_FOLDER = '/var/www/kgp/FlaskApp/uploads'
ALLOWED_EXTENSIONS = set(['csv'])

app = Flask(__name__)
app.secret_key='DGLinker'
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

def allowed_file(filename):
	return '.' in filename and filename.rsplit('.', 1)[1] in ALLOWED_EXTENSIONS

def send_error(email, tempid, error):
	if email != '':
		yag = yagmail.SMTP( user="dglinker.service@gmail.com", password="DGLinker!", host='smtp.gmail.com')
		subject = "DGLinker: There is something wrong with your project."
		contents = "Dear user,\n\n  There is something wrong with your project： \n https://dglinker.rosalind.kcl.ac.uk/result/%s \n" % tempid
		contents = contents + str(error)
		contents = contents + '\n Your project will not continue. Please submit a new one. \n\n Best wishes \n DGLinker Service'
		yag.send(email, subject, contents)

def error_page(tempid, err):
	f1 = open('/var/www/kgp/FlaskApp/cache/error_file%s.txt' % tempid,'w')
	f1.write(str(err))
	f1.close()

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
		f1.write(str(request.form.getlist('filename')))

		email = str(input_dict['email'])
		file_name = '/var/www/kgp/FlaskApp/cache/output_file%s.csv' % tempid
		file_list = request.form.getlist('filename')
		uploaded_files = request.files.getlist("file[]")
		k_value = str(input_dict['kvalue'])


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
		thr = Thread(target=predict, args=(tempid, email, file_name_list, file_name, gene_list, phenotype, input_name,keywords,k_value,))
		thr.start()
		return redirect(url_for('result',tempid=tempid))

def send_start(email, tempid):
	yag = yagmail.SMTP( user="dglinker.service@gmail.com", password="DGLinker!", host='smtp.gmail.com')
	subject = "DGLinker: Your project has been submitted."
	contents = "Dear user,\n\n Your project has been submited. You can find your result at the following link when the job will be completed: \n https://dglinker.rosalind.kcl.ac.uk/result/%s \n\n Best wishes, \n The DGLinker Team" % tempid
	yag.send(email, subject, contents)

def predict(tempid, email, file_name_list, file_name, gene_list, phenotype, input_name,keywords,k_value):
	if phenotype == "pheno":
		voca = pd.read_csv('/var/www/kgp/FlaskApp/data/vocabulary.csv',header=0,index_col=0)
		f1 = open('/var/www/kgp/FlaskApp/data/log.txt','a')
		matched_diseidlist = []
		phenotype_used = []
		for i in input_name:
			if i != '':
				if keywords == "True":
					vocadict = open("/var/www/kgp/FlaskApp/data/vocabulary.csv", 'r')
					for row in vocadict:
						row = row.strip().split(',')
						if str(i).strip().lower() in row[1].strip().lower():
							matched_diseidlist.append(row[0])
							phenotype_used.append(row[1])
				else:
					matched_diseid = voca.index[voca['diseaseName'] == str(i)].tolist()
					matched_diseidlist.append(matched_diseid[0])

		f2 = open(file_name,'r') #dataset, 5 columns
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
		matched_diseidlist = []
		phenotype_used = []
		f1 = open('/var/www/kgp/FlaskApp/data/log.txt','a')
		for i in input_name:
			if i != '':
				if keywords == "True":
					vocadict = open("/var/www/kgp/FlaskApp/data/vocabulary.csv", 'r')
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
		with open(merge_filename,'w') as new:
			new.write('Source node name,Source node type,Relationship type,Target node type,Target node name\n')
			for i in gene_list:
				if str(i) != "":
					new.write('%s,%s,%s,%s,%s\n' % (str(i),'gene','gene_disease_association','disease',str(matched_diseidlist[0])))
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
				new.write('%s,%s,%s,%s,%s\n' % (str(i),'gene','gene_disease_association','disease','USER_DEFINED_DISEASE'))
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
	
	ep = EdgePrediction.EdgePrediction()

	try:
		ep.CSV_to_graph(fname = merge_filename)
	except Exception as err:
		send_error(email, tempid, err)
		f1 = open('/var/www/kgp/FlaskApp/cache/error_file%s.txt' % tempid,'w')
		f1.write(str(err))
		f1.close()
		return render_template('error.html', tempid=tempid, error=err)

	try:
		ep.preprocess()
	except Exception as err:
		send_error(email, tempid, err)
		f1 = open('/var/www/kgp/FlaskApp/cache/error_file%s.txt' % tempid,'w')
		f1.write(str(err))
		f1.close()
		return render_template('error.html', tempid=tempid, error=err)

	ep.to_predict = 'gene_disease_association'
	# ep.network_order = order_list  # input_dict['order']  # ['HAS_SIDE_EFFECT', 'DRUG_TARGETS', 'INDICATED_FOR']

	if k_value:
		xv_results = []
		xv = ep.k_fold(target = str(target_name), k = int(k_value), calculate_auc = False) # , int(k_value)
		for fold in xv:
			fold['objective'] = "J"

		built = [x for x in xv if x['model_built'] == True]
		model_performance = open('/var/www/kgp/FlaskApp/model_results/model_performance_%s.txt' % tempid, 'a')
		model_performance.write("%s\t%s\n" % ('The number of successful folds','%d' % len(built)))
		for fold in built:
			fold.update(fold['contingency'])
			fold['n_deleted_predicted'] = sum(fold['left_out_predicted'])
			del fold['contingency']
			del fold['left_out_predicted']
			del fold['left_out']

		df = pd.DataFrame(built)
		df['n_pred_train'] = df['tp'] + df['fp']
		df['train_prec'] = df['tp'] / df['n_pred_train']
		df['train_rec'] = df['tp'] / df['n_known_train']
		n_neg = df['tn'] + df['fp']
		df['train_fpr'] = df['fp'] / n_neg

		# a good model has low FPR and high test recall (proportion_predicted)
		onemins_fpr = 1 - df['train_fpr']
		su = onemins_fpr + df['proportion_predicted']
		pr = onemins_fpr * df['proportion_predicted']
		df['overall'] = 2 * (pr/su)
		df['prediction_hit_rate'] = df['n_deleted_predicted'] / df['fp'] #of FP how many are really linked
		total_neg = df['tn'] + df['fp']
		df['expected_hit_rate'] = df['n_known_test'] / total_neg
		df['fc_hit_rate'] = df['prediction_hit_rate'] / df['expected_hit_rate']
		# output_per_fold = '/var/www/kgp/FlaskApp/cache/xv_per_fold_%s.csv' % tempid
		# df.to_csv(output_per_fold, index=False)

		#average results over all folds
		gr = df.groupby('target')

		keep = ['objective_performance', 'train_prec', 'train_rec', 'train_fpr', 
				'proportion_predicted', 'overall', 'prediction_hit_rate', 'expected_hit_rate',
				'fc_hit_rate']
		xv_avg = gr[keep].mean()
		xv_avg.reset_index(inplace=True)
		model_performance.write("%s\t%s\n" % ('The average performance improvement vs random', '%.3f' % xv_avg['fc_hit_rate']))

	try:
		result = ep.predict(target=target_name, calculate_auc=True, return_scores=True)
	except Exception as err:
		send_error(email, tempid, err)
		f1 = open('/var/www/kgp/FlaskApp/cache/error_file%s.txt' % tempid,'w')
		f1.write(str(err))
		f1.close()
		return render_template('error.html', tempid=tempid, error=err)
	
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
	model_performance = open('/var/www/kgp/FlaskApp/model_results/model_performance_%s.txt' % tempid, 'a')
	result_file.write("%s\t%s\n" % ('Project ID',str(tempid)))
	job_summary.write("%s\t%s\n" % ('Project ID',str(tempid)))
	if email != '':
		result_file.write("%s\t%s\n" % ('Result Sent to',email))
		job_summary.write("%s\t%s\n" % ('Project ID',str(tempid)))
	result_file.write("%s\t%s\n" % ("Dataset Used", str([i for i in file_name_list]).replace('.csv','')))
	job_summary.write("%s\t%s\n" % ("Dataset Used", ', '.join([i for i in file_name_list]).replace('.csv','')))
	if keywords == 'True':
		result_file.write("%s\t%s\n" % ("Keywords Used", "True"))
		job_summary.write("%s\t%s\n" % ("Keywords Used", "True"))
	else:
		result_file.write("%s\t%s\n" % ("Keywords Used", "False"))
		job_summary.write("%s\t%s\n" % ("Keywords Used", "False"))
	if phenotype == "pheno":
		if keywords == 'True':
			result_file.write("%s\t%s\n" % ('Predict Mode', 'Predict from Phenotype(s)'))
			result_file.write("%s\t%s\n" % ('Phenotype(s) Used', str([i for i in phenotype_used])))
			job_summary.write("%s\t%s\n" % ('Predict Mode', 'Predict from Phenotype(s)'))
			job_summary.write("%s\t%s\n" % ('Phenotype(s) Used', ', '.join(([str(i) for i in phenotype_used]))))
		else:
			result_file.write("%s\t%s\n" % ('Predict Mode', 'Predict from Phenotype(s)'))
			result_file.write("%s\t%s\n" % ('Phenotype(s) Used', str([i for i in input_name])))
			job_summary.write("%s\t%s\n" % ('Predict Mode', 'Predict from Phenotype(s)'))
			job_summary.write("%s\t%s\n" % ('Phenotype(s) Used', ', '.join(([str(i) for i in input_name]))))
		result_file.write("%s\t%s\n" % ('Known Genes',result['known_hits']))
	elif phenotype == "defined":
		if keywords == 'True':
			result_file.write("%s\t%s\n" % ('Predict Mode', "Predict from User-defined Disease(s)"))
			result_file.write("%s\t%s\n" % ('Disease(s) Defined', str([i for i in phenotype_used])))
			result_file.write("%s\t%s\n" % ('By Genes', str(gene_list)))
			job_summary.write("%s\t%s\n" % ('Predict Mode', "Predict from User-defined Disease(s)"))
			job_summary.write("%s\t%s\n" % ('Disease(s) Defined', ', '.join(([str(i) for i in phenotype_used]))))
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

	# result_file.write("%s\t%s\n" % ('ACC of the Model','%.3f' % result['ACC']))
	# result_file.write("%s\t%s\n" % ('AUC of the Model','%.3f' % result['auc']))
	# result_file.write("%s\t%s\n" % ('F1 score of the Model','%.3f' % result['F1']))
	model_performance.write("%s\t%s\n" % ('ACC of the Model','%.3f' % result['ACC']))
	model_performance.write("%s\t%s\n" % ('AUC of the Model','%.3f' % result['auc']))
	model_performance.write("%s\t%s\n" % ('F1 score of the Model','%.3f' % result['F1']))


	result_file.write("%s\n" % 'More Information About Genes')
	download_results.write("%s\n" % 'Rank\tAssociation type\tGene name\tScore')
	# rank for all
	ranking = sorted(result['scores']['scores'].items(), key=lambda item: item[1], reverse=True)
	num = 0
	top_new = []
	for tup in ranking:
		if tup[0] in known_predictions:
			result_file.write("%d\t%s\t%s\t%.3f\n" % (num+1,'Known',tup[0],tup[1]))
			download_results.write("%d\t%s\t%s\t%.3f\n" % (num+1,'Known',tup[0],tup[1]))
		elif tup[0] in new_predictions:
			result_file.write("%d\t%s\t%s\t%.3f\n" % (num+1,'Predicted',tup[0],tup[1]))
			download_results.write("%d\t%s\t%s\t%.3f\n" % (num+1,'Predicted',tup[0],tup[1]))
			top_new.append(tup[0])
		num += 1
	result_file.flush()
	result_file.close()
	job_summary.flush()
	job_summary.close()
	download_results.flush()
	download_results.close()

	# top 100 for enrichr
	top_100_new = top_new[0:100]
	# graph file for top_50_all
	top_50_all = [tup[0] for tup in ranking][0:50]

	# pie chart
	parent_dir = '/var/www/kgp/FlaskApp/static/'
	directory = 'pie_chart_%s' % tempid
	path = os.path.join(parent_dir, directory)
	if not os.path.exists(path):
		os.mkdir(path)
	for i in top_50_all:
		info = result['scores']['breakdown'][i]
		pie_chart(tempid, info, i)
	pie_chart(tempid, result['weights'], 'Total')

	# network graph
	line_duplicate = []
	raw_file = open('/var/www/kgp/FlaskApp/cache/output_file%s.csv' % tempid, 'r')
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

	#enrichr
	# make directory for enrichr
	parent_dir = '/var/www/kgp/FlaskApp/static/'
	directory = 'enrichr%s' % tempid
	path = os.path.join(parent_dir, directory)
	if not os.path.exists(path):
		os.mkdir(path)

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
		yag = yagmail.SMTP( user="dglinker.service@gmail.com", password="DGLinker!", host='smtp.gmail.com')
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
		err = errors[0]
		return render_template('error.html', tempid=tempid, error=err)
	if os.path.exists('/var/www/kgp/FlaskApp/static/enrichr%s/graphdone.txt' % tempid):
		final_file = open('/var/www/kgp/FlaskApp/cache/result_file%s.txt' % tempid,'r')
		model_performance = open('/var/www/kgp/FlaskApp/model_results/model_performance_%s.txt' % tempid, 'r')
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

		performance = OrderedDict()
		for line in model_performance:
			line = line.strip().split('\t')
			performance[line[0]] = line[1]

		return render_template('result.html',result=result, all_hits=all_hits, tempid=tempid, performance=performance)
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
	

def pie_chart(tempid, info, i):
	graphname = '/var/www/kgp/FlaskApp/static/pie_chart_%s/%s_%s.png' % (tempid,tempid,i)
	plt.figure()
	labels = [i for i in info.keys()]
	sizes = [i/sum(info.values())*100 for i in info.values()]
	colors = ['lightsalmon', 'peachpuff', 'lightblue', 'bisque', 'plum', 'khaki', 'turquoise'][:len(sizes)]
	patches,text1,text2 = plt.pie(sizes,
						labels=labels,
						colors=colors,
						labeldistance = 1.1,
						autopct = '%3.2f%%',
						shadow = False,
						startangle =90,
						pctdistance = 0.6)
	plt.axis('equal')
	plt.title(str(i))
	# plt.legend()
	plt.tight_layout()
	plt.savefig(graphname, dpi=100)
	plt.clf()
	plt.close()

	
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
		filelist1 = ['/var/www/kgp/FlaskApp/model_results/job_summary_%s.txt' % tempid, '/var/www/kgp/FlaskApp/model_results/results_%s.txt' % tempid, '/var/www/kgp/FlaskApp/model_results/graph_file%s.txt' % tempid, '/var/www/kgp/FlaskApp/model_results/model_performance_%s.txt' % tempid]
		dirpath1 = '/var/www/kgp/FlaskApp/static/enrichr%s/' % tempid
		dirpath2 = '/var/www/kgp/FlaskApp/static/graph_dir%s/' % tempid
		dirpath3 = '/var/www/kgp/FlaskApp/static/pie_chart_%s/' % tempid
		filelist2 = search(dirpath1, '_')
		filelist3 = search(dirpath2, '.txt')
		filelist4 = search(dirpath2, '_')
		filelist5 = search(dirpath3, '.png')

		for f in filelist1 + filelist2 + filelist3 + filelist4 + filelist5:
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



