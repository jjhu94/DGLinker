import os
os.environ['MPLCONFIGDIR'] = "/var/www/kgp/FlaskApp/graph"
import matplotlib
matplotlib.use('Agg')


from collections import OrderedDict
from collections import defaultdict
from collections import Counter
from .Enrichment import Enricher
from flask import flash, Flask, make_response, redirect, render_template, request, send_file, url_for
from functools import reduce
from .net_visualisation import net_visualization
from pyvis.network import Network
from scipy.stats import hypergeom
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
plt.rcParams.update({'figure.max_open_warning': 0})
import os
import operator
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
		contents = "Dear user,\n\n  There is something wrong with your project: \n https://dglinker.rosalind.kcl.ac.uk/result/%s \n" % tempid
		contents = contents + str(error)
		contents = contents + '\n Your project will not continue. Please submit a new one. \n\n Best wishes \n DGLinker Service'
		yag.send(email, subject, contents)

def error_page(tempid, err):
	f1 = open('/var/www/kgp/FlaskApp/cache/error_file%s.txt' % tempid,'w')
	f1.write(str(err))
	f1.close()
	

def hyper_prob_at_least(pop_n, pop_true, draw_n, draw_true):
    #prob of at least h hits is i - cdf(h-1)
    prb = hypergeom.cdf(draw_true-1, pop_n, pop_true, draw_n)
    prb_h_plus = 1 - prb
    return prb_h_plus

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
		# f1.write(str(request.form.getlist('filename')))

		email = str(input_dict['email'])
		file_name = '/var/www/kgp/FlaskApp/cache/output_file%s.csv' % tempid
		file_list = request.form.getlist('filename')
		uploaded_files = request.files.getlist("uploadfile")
		# f1.write(str(uploaded_files))
		try:
			k_value = str(input_dict['kvalue'])
		except:
			k_value = False

		for file in uploaded_files:
			if file and allowed_file(file.filename):
				filename = secure_filename(file.filename)
				file.save(os.path.join(app.config['UPLOAD_FOLDER'],filename))
				file_list.append(os.path.join(app.config['UPLOAD_FOLDER'],filename))
		file_name_list = []
		for f in file_list:
			f = f.strip().split('/')
			file_name_list.append(f[-1])

		f1.write(str(file_list))

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
	contents = "Dear user,\n\n Your project has been submited. You can find your result at the following link when the job will be completed: \n https://dglinker.rosalind.kcl.ac.uk/result/%s \nPlease note that the results will be deleted after 14 days. If you want to keep the results, please download them. \n\n Best wishes, \n The DGLinker Team" % tempid
	yag.send(email, subject, contents)

def queue_system(tempid):
	while True:
		queuing = open('/var/www/kgp/FlaskApp/data/queue.txt','r')
		start = 0
		done = 0
		for line in queuing:
			if 'start' in line:
				start += 1
			elif 'done' in line:
				done += 1
		if start - done <= 5:
			queuing.close()
			break
		else:
			queuing.flush()
			time.sleep(120)
			f1 = open('/var/www/kgp/FlaskApp/data/wait.txt','a')
			f1.write(str(tempid)+' waiting...\n')
			f1.close()

			
def predict(tempid, email, file_name_list, file_name, gene_list, phenotype, input_name,keywords,k_value):
	queue_system(tempid)
	queue = open('/var/www/kgp/FlaskApp/data/queue.txt','a')
	queue.write(str(tempid)+' start\n')
	queue.close()

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
							if row[1] not in phenotype_used:
								phenotype_used.append(row[1])
				else:
					try:
						matched_diseid = voca.index[voca['diseaseName'] == str(i)].tolist()
						matched_diseidlist.append(matched_diseid[0])
					except:
						pass

		f2 = open(file_name,'r') #dataset, 5 columns
		merge_filename = "/var/www/kgp/FlaskApp/cache/merged_file%s.csv" % tempid
		with open(merge_filename,'w') as new:
			new.write('Source node name,Source node type,Relationship type,Target node type,Target node name\n')
			for line in f2:
				line_list = line.strip().split(',')
				try:
					if line_list[4] in matched_diseidlist:
						new.write('%s,%s,%s,%s,%s\n' % (line_list[0],line_list[1].lower(),line_list[2].lower(),line_list[3].lower(),matched_diseidlist[0]))
					else:
						new.write('%s,%s,%s,%s,%s\n' % (line_list[0],line_list[1].lower(),line_list[2].lower(),line_list[3].lower(),line_list[4]))
				except:
					pass
	elif phenotype == "defined":
		voca = pd.read_csv('/var/www/kgp/FlaskApp/data/vocabulary.csv',header=0,index_col=0)
		matched_diseidlist = []
		phenotype_used = []
		for i in input_name:
			if i != '':
				if keywords == "True":
					vocadict = open("/var/www/kgp/FlaskApp/data/vocabulary.csv", 'r')
					for row in vocadict:
						row = row.strip().split(',')
						if str(i).strip().lower() in row[1].lower():
							matched_diseidlist.append(row[0])
							if row[1] not in phenotype_used:
								phenotype_used.append(row[1])
							
				else:
					try:
						matched_diseid = voca.index[voca['diseaseName'] == str(i)].tolist()
						matched_diseidlist.append(matched_diseid[0])
					except:
						pass
		f2 = open(file_name,'r') #dataset,5 columns
		merge_filename = "/var/www/kgp/FlaskApp/cache/merged_file%s.csv" % tempid
		with open(merge_filename,'w') as new:
			new.write('Source node name,Source node type,Relationship type,Target node type,Target node name\n')
			for i in gene_list:
				if str(i) != "":
					new.write('%s,%s,%s,%s,%s\n' % (str(i),'gene','gene_disease_association','disease',str(matched_diseidlist[0])))
			for line in f2:
				line_list = line.strip().split(',')
				try:
					if line_list[4] not in matched_diseidlist: # if in, pass
						new.write('%s,%s,%s,%s,%s\n' % (line_list[0],line_list[1].lower(),line_list[2].lower(),line_list[3].lower(),line_list[4]))
				except:
					pass
	elif phenotype == "from_genes":
		f2 = open(file_name,'r') #dataset,5 columns
		merge_filename = "/var/www/kgp/FlaskApp/cache/merged_file%s.csv" % tempid
		with open(merge_filename,'w') as new:
			new.write('Source node name,Source node type,Relationship type,Target node type,Target node name\n')
			for i in gene_list:
				new.write('%s,%s,%s,%s,%s\n' % (str(i),'gene','gene_disease_association','disease','USER_DEFINED_DISEASE'))
			for line in f2:
				line_list = line.strip().split(',')
				try:
					new.write('%s,%s,%s,%s,%s\n' % (line_list[0],line_list[1].lower(),line_list[2].lower(),line_list[3].lower(),line_list[4]))
				except:
					pass

	order_dict = pd.read_csv(merge_filename,header=0,index_col=0, error_bad_lines=False, low_memory=False)
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
		queue = open('/var/www/kgp/FlaskApp/data/queue.txt','a')
		queue.write(str(tempid)+' done\n')
		queue.close()
		time.sleep(20)
		return render_template('error.html', tempid=tempid, error=err)

	try:
		ep.preprocess()
	except Exception as err:
		send_error(email, tempid, err)
		f1 = open('/var/www/kgp/FlaskApp/cache/error_file%s.txt' % tempid,'w')
		f1.write(str(err))
		f1.close()
		queue = open('/var/www/kgp/FlaskApp/data/queue.txt','a')
		queue.write(str(tempid)+' done\n')
		queue.close()
		time.sleep(20)
		return render_template('error.html', tempid=tempid, error=err)

	ep.to_predict = 'gene_disease_association'
	# ep.network_order = order_list  # input_dict['order']  # ['HAS_SIDE_EFFECT', 'DRUG_TARGETS', 'INDICATED_FOR']

	# k-fold validation
	if k_value:
		xv_results = []
		xv = ep.k_fold(target = str(target_name), k = int(k_value), calculate_auc = False) # , int(k_value)
		if xv:
			for fold in xv:
				fold['objective'] = "J"
			built = [x for x in xv if x['model_built'] == True]
			model_performance = open('/var/www/kgp/FlaskApp/model_results/model_performance_%s.txt' % tempid, 'a')
			model_performance.write("%s\t%s\n" % ('The total number of folds','%d' % int(k_value)))
			model_performance.write("%s\t%s\n" % ('The number of folds could be trained','%d' % len(built)))

			for result in built:
				result['n_deleted_predicted'] = sum(result['left_out_predicted'])
				pop_n = result['contingency']['tn'] + result['contingency']['fp']
				pop_true = result['n_known_test']
				draw_n = result['contingency']['fp']
				draw_true = result['n_deleted_predicted']
				prob = hyper_prob_at_least(pop_n, pop_true, draw_n, draw_true)
				result['prob'] = prob
				result['signif'] = prob < 0.05

			is_significant = [x['signif'] for x in built]
			model_performance = open('/var/www/kgp/FlaskApp/model_results/model_performance_%s.txt' % tempid, 'a')
			model_performance.write("%s\t%s\n" % ('The number of successful and significantly enriched folds','%d' % sum(is_significant)))

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
	else:
		model_performance = open('/var/www/kgp/FlaskApp/model_results/model_performance_%s.txt' % tempid, 'a')
		model_performance.write("The total number of folds\t0\nThe number of folds could be trained\t0\nThe number of successful and significantly enriched folds\t0\nThe average performance improvement vs random\t0\n")

	input_genes = ep.getKnown(target=target_name)

	# to predict
	try:
		result = ep.predict(target=target_name, calculate_auc=True, return_scores=True)
	except Exception as err:
		send_error(email, tempid, err)
		f1 = open('/var/www/kgp/FlaskApp/cache/error_file%s.txt' % tempid,'w')
		f1.write(str(err))
		f1.close()
		queue = open('/var/www/kgp/FlaskApp/data/queue.txt','a')
		queue.write(str(tempid)+' done\n')
		queue.close()
		time.sleep(20)
		return render_template('error.html', tempid=tempid, error=err)
	
	if not result['model_built']:
		err = 'The selected networks have 0 predictor overlap. You can select more databases in the same type and have another try.'
		send_error(email, tempid, err)
		f1 = open('/var/www/kgp/FlaskApp/cache/error_file%s.txt' % tempid,'w')
		f1.write(str(err))
		f1.close()
		queue = open('/var/www/kgp/FlaskApp/data/queue.txt','a')
		queue.write(str(tempid)+' done\n')
		queue.close()
		time.sleep(20)
		return render_template('error.html', tempid=tempid, error=err)

	if float(result['auc']) < 0.5:
		err = "The AUC of the model built is less than 0.5."
		send_error(email, tempid, err)
		f1 = open('/var/www/kgp/FlaskApp/cache/error_file%s.txt' % tempid,'w')
		f1.write(str(err))
		f1.close()
		queue = open('/var/www/kgp/FlaskApp/data/queue.txt','a')
		queue.write(str(tempid)+' done\n')
		queue.close()
		time.sleep(20)
		return render_template('error.html', tempid=tempid, error=err)	
	# get result	
	new_predictions = result['new_hits']
	known_predictions = result['known_hits']
	all_predictions = result['all_hits']
	weights = result['weights']
	scores = ep.getScores(target_name, weights)

	# bed files
	directory = '/var/www/kgp/FlaskApp/model_results/results_genomic_regions_%s' % tempid
	if not os.path.exists(directory):
		os.mkdir(directory)
	hg_19 = pd.read_csv('/var/www/kgp/FlaskApp/data/genes_hg19.bed', sep='\t',header=None,names=['hg19','gene'])
	hg_38 = pd.read_csv('/var/www/kgp/FlaskApp/data/genes_hg38.bed', sep='\t',header=None,names=['hg38','gene'])
	input_genes_hg19 = hg_19[hg_19['gene'].isin(input_genes)]
	input_genes_hg38 = hg_38[hg_38['gene'].isin(input_genes)]
	predicted_genes_hg19 = hg_19[hg_19['gene'].isin(all_predictions)]
	predicted_genes_hg38 = hg_38[hg_38['gene'].isin(all_predictions)]
	input_genes_hg19.to_csv('%s/input_genes_hg19.bed' % directory, sep='\t',index=False)
	input_genes_hg38.to_csv('%s/input_genes_hg38.bed' % directory, sep='\t',index=False)
	predicted_genes_hg19.to_csv('%s/predicted_genes_hg19.bed' % directory, sep='\t',index=False)
	predicted_genes_hg38.to_csv('%s/predicted_genes_hg38.bed' % directory, sep='\t',index=False)

	# job summary and model performance
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
		result_file.write("%s\t%s\n" % ('Input Known Genes',input_genes))
		job_summary.write("%s\t%s\n" % ('Input Known Genes',', '.join(input_genes)))
		result_file.write("%s\t%s\n" % ('Number of Input Known Genes',len(input_genes)))
		job_summary.write("%s\t%s\n" % ('Number of Input Known Genes',len(input_genes)))
		

	elif phenotype == "defined":
		if keywords == 'True':
			result_file.write("%s\t%s\n" % ('Predict Mode', "Predict from User-defined Disease(s)"))
			result_file.write("%s\t%s\n" % ('Disease(s) Defined', str([i for i in phenotype_used])))
			result_file.write("%s\t%s\n" % ('By Genes', str(input_genes)))
			result_file.write("%s\t%s\n" % ('Number of Genes Used',len(input_genes)))
			job_summary.write("%s\t%s\n" % ('Predict Mode', "Predict from User-defined Disease(s)"))
			job_summary.write("%s\t%s\n" % ('Disease(s) Defined', ', '.join(([str(i) for i in phenotype_used]))))
			job_summary.write("%s\t%s\n" % ('By Genes', str(input_genes)))
			job_summary.write("%s\t%s\n" % ('Number of Genes Used',len(input_genes)))
			
		else:
			result_file.write("%s\t%s\n" % ('Predict Mode', "Predict from User-defined Disease(s)"))
			result_file.write("%s\t%s\n" % ('Disease(s) Defined', str(input_name)))
			result_file.write("%s\t%s\n" % ('By Genes', str(input_genes)))
			result_file.write("%s\t%s\n" % ('Number of Genes Used',len(input_genes)))
			job_summary.write("%s\t%s\n" % ('Predict Mode', "Predict from User-defined Disease(s)"))
			job_summary.write("%s\t%s\n" % ('Disease(s) Defined', str(input_name)))
			job_summary.write("%s\t%s\n" % ('By Genes', str(input_genes)))
			job_summary.write("%s\t%s\n" % ('Number of Genes Used',len(input_genes)))

	elif phenotype == "from_genes":
		result_file.write("%s\t%s\n" % ('Predict Mode', "Predict from Gene(s)"))
		result_file.write("%s\t%s\n" % ('Genes Used', str(input_genes)))
		result_file.write("%s\t%s\n" % ('Number of Genes Used', len(input_genes)))
		job_summary.write("%s\t%s\n" % ('Predict Mode', "Predict from Gene(s)"))
		job_summary.write("%s\t%s\n" % ('Genes Used', str(input_genes)))
		job_summary.write("%s\t%s\n" % ('Number of Genes Used', len(input_genes)))

	result_file.write("%s\t%s\n" % ('Total Number of Predicted Genes','%d' %  result['hits_total']))
	result_file.write("%s\t%s\n" % ('Number of New Predicted Genes','%d' %  result['hits_new']))
	hits_known = result['hits_total'] - result['hits_new']
	result_file.write("%s\t%s\n" % ('Number of Known Predicted Genes', '%d' % hits_known))
	job_summary.write("%s\t%s\n" % ('Total Number of Predicted Genes','%d' %  result['hits_total']))
	job_summary.write("%s\t%s\n" % ('Number of New Predicted Genes','%d' %  result['hits_new']))
	job_summary.write("%s\t%s\n" % ('Number of Known Predicted Genes', '%d' % hits_known))

	# result_file.write("%s\t%s\n" % ('ACC of the Model','%.3f' % result['ACC']))
	# result_file.write("%s\t%s\n" % ('AUC of the Model','%.3f' % result['auc']))
	# result_file.write("%s\t%s\n" % ('F1 score of the Model','%.3f' % result['F1']))
	model_performance.write("%s\t%s\n" % ('ACC of the Model','%.3f' % result['ACC']))
	model_performance.write("%s\t%s\n" % ('AUC of the Model','%.3f' % result['auc']))
	model_performance.write("%s\t%s\n" % ('F1 score of the Model','%.3f' % result['F1']))
	model_performance.write("%s\t%s\n" % ('J of the Model','%.3f' % result['J']))


	# ENSP ID
	with open("/var/www/kgp/FlaskApp/data/genes_ENSP.txt", 'r') as f:
		enspdict = {}
		for row in f:
			row = row.strip().split('\t')
			enspdict[row[0]] = row[1]


	result_file.write("%s\n" % 'More Information About Genes')
	download_results.write("%s\n" % 'Rank\tAssociation type\tGene name\tEnsemble protein ID\tScore')
	# rank for all
	ranking = sorted(result['scores']['scores'].items(), key=lambda item: item[1], reverse=True)
	num = 0
	top_new = []
	for tup in ranking:
		if tup[0] in known_predictions:
			association = 'Known'
		elif tup[0] in new_predictions:
			association = 'Predicted'
			top_new.append(tup[0])
		try:
			result_file.write("%d\t%s\t%s\t%s\t%.3f\n" % (num+1,association,tup[0],enspdict[tup[0]],tup[1]))
			download_results.write("%d\t%s\t%s\t%s\t%.3f\n" % (num+1,association,tup[0],enspdict[tup[0]],tup[1]))
		except:
			result_file.write("%d\t%s\t%s\t%s\t%.3f\n" % (num+1,association,tup[0],'-',tup[1]))
			download_results.write("%d\t%s\t%s\t%s\t%.3f\n" % (num+1,association,tup[0],'-',tup[1]))
			
		num += 1
	result_file.flush()
	result_file.close()
	job_summary.flush()
	job_summary.close()
	download_results.flush()
	download_results.close()
	model_performance.flush()
	model_performance.close()
	
	# top 100 for enrichr
	top_100_new = top_new[0:100]
	# graph file for top_50_all
	top_50_all = [tup[0] for tup in ranking][0:50]

	
	# get all genes:
	# for key in result['scores']['scores']:

	
	# pie chart
	parent_dir = '/var/www/kgp/FlaskApp/static/'
	directory = 'pie_chart_%s' % tempid
	path = os.path.join(parent_dir, directory)
	if not os.path.exists(path):
		os.mkdir(path)
	
	try:
		pie_chart(tempid, result['weights'], 'Total')
	except:
		pass

	for i in all_predictions:
		info = result['scores']['breakdown'][i]
		try:
			pie_chart(tempid, info, i)
		except:
			pass
	
	
	# network graph
	line_duplicate = []
	raw_file = open('/var/www/kgp/FlaskApp/cache/output_file%s.csv' % tempid, 'r')
	with open('/var/www/kgp/FlaskApp/cache/raw_graph_file%s.csv' % tempid, 'w') as raw_graph_file:
		for line in raw_file:
			linelist = line.strip().split(',')
			if linelist[0] in all_predictions and line not in line_duplicate:
				raw_graph_file.write(line)

	# get relationships that actually used in networks
	info = result['weights']
	actual_edges = []
	a = list(info.keys())
	b = list(info.values())
	for i in range(len(b)):
		if b[i] > 0:
			actual_edges.append(a[i])

	with open('/var/www/kgp/FlaskApp/cache/raw_graph_file%s.csv' % tempid, 'r') as reader:
		d = OrderedDict([('Source node name', []), ('Source node type', []), ('Target node type', []), ('Target node name', [])])
		for row in reader:
			row = row.strip().split(',')
			if row[2] in actual_edges:
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
				vocadict[row[0]] = row[1] # vocadict['C000965'] = 'AD'

	genecount = Counter(d['Source node name'])
	othercount = Counter(d['Target node name'])
	for i in range(len(d['Source node name'])):
		if d['Source node name'][i] in new_predictions:
			genetype = 'Predicted'
		else:
			genetype = 'Known'
		if d['Target node type'][i] == 'disease' and othercount[d['Target node name'][i]] > 1:
			try:
				graph_file.write('%s\t%s\t%s\t%s\t%d\t%d\n' % (d['Source node name'][i],vocadict[d['Target node name'][i]],genetype,d['Target node type'][i],genecount[d['Source node name'][i]],othercount[d['Target node name'][i]]))
				download_graph_file.write('%s\t%s\t%s\t%s\t%d\t%d\n' % (d['Source node name'][i],vocadict[d['Target node name'][i]],genetype,d['Target node type'][i],genecount[d['Source node name'][i]],othercount[d['Target node name'][i]]))
			except:
				pass
		elif d['Target node type'][i] != 'disease' and othercount[d['Target node name'][i]] > 1:
			try:
				graph_file.write('%s\t%s\t%s\t%s\t%d\t%d\n' % (d['Source node name'][i],str(d['Target node name'][i]),genetype,d['Target node type'][i],genecount[d['Source node name'][i]],othercount[d['Target node name'][i]]))
				download_graph_file.write('%s\t%s\t%s\t%s\t%d\t%d\n' % (d['Source node name'][i],str(d['Target node name'][i]),genetype,d['Target node type'][i],genecount[d['Source node name'][i]],othercount[d['Target node name'][i]]))
			except:
				pass
	graph_file.close()
	download_graph_file.close()

	if not os.path.exists('/var/www/kgp/FlaskApp/static/graph_dir%s' % tempid):
		os.mkdir('/var/www/kgp/FlaskApp/static/graph_dir%s' % tempid)
	
	graph = pd.read_csv('/var/www/kgp/FlaskApp/cache/graph_file%s.txt' % tempid,sep='\t', names=['Gene', 'Others', 'Gene_class', 'Others_Class', 'Total number of interactions gene', 'Total number of interactions others'])
	# graph_50 = graph[graph['Gene'].isin(top_50_all)]
	predictor = reduce(operator.add,list(result['predictors'].values()))
	f1 = open('/var/www/kgp/FlaskApp/data/log.txt','a')
	f1.write(str(predictor))
	f1.close()
	predictor_name = []
	for i in predictor:
		try:
			predictor_name.append(vocadict[i])
		except:
			predictor_name.append(i)

	graph_predictor = graph[graph['Others'].isin(predictor_name) & graph['Gene'].isin(input_genes)] # subgraphs to be get
	
	for gene in top_50_all: # result['all_hits']
		# data = graph[graph['Gene'] == gene]
		data = graph[graph['Gene'].isin([gene]) & graph['Others'].isin(predictor_name)] # get interaction related to the selected gene and predictors
		# subgraph = graph_50[graph_50['Others'].isin(list(data['Others']))] # get interaction related to gene ralated others
		subgraph = graph_predictor[graph_predictor['Others'].isin(list(data['Others']))] # get interaction related to input gene ralated others
		data = data.append(subgraph, ignore_index=True)
		data = data.drop_duplicates()
		# data = data[data.groupby('Others').Others.transform(len) > 1] # delete others which only appear once
		data.to_csv('/var/www/kgp/FlaskApp/static/graph_dir%s/%s.txt' % (tempid,gene),index=False,sep='\t')

	folder = '/var/www/kgp/FlaskApp/static/graph_dir%s' % tempid
	threa = Thread(target=net_visualization, args=(folder,False,))
	threa.start()
	
	#enrichr
	# make directory for enrichr
	parent_dir = '/var/www/kgp/FlaskApp/static/'
	directory = 'enrichr%s' % tempid
	path = os.path.join(parent_dir, directory)
	if not os.path.exists(path):
		os.mkdir(path)

	namelist = ['New', 'Known', 'Union']
	gene_list = [top_100_new,known_predictions,all_predictions]
	files = [
			'GWAS v1.2.csv',
			'OMIM v2021-04.csv',
			'KEGG v2021-03.csv',
			'DSigDB v1.0.csv',
			'ENCODE and ChEA Consensus TFs from ChIP-X v1.0.csv',
			'ArrayExpress Atlas (experiment E-MTAB-513 Illumina body map) v2021-02.csv',
		]
	
	enr = Enricher(backgrounds_files=files)
	files = files + ['GO']
	go_names = ['GO Biological Process v2021-04','GO Cellular Component v2021-04','GO Molecular Function v2021-04']
	abstr = ['BP','CC','MF']
	for i in range(3):
		for db in files:
			db = db.replace('.csv','')
			res = enr.enrichment(gene_list[i], return_filter='all', dataset=db)
			enr.to_tsv('/var/www/kgp/FlaskApp/static/enrichr%s/%s_%s_%s.txt' % (tempid,tempid,namelist[i],db), res)
			if db == 'GO':
				for j in range(len(go_names)):
					data = pd.read_csv('/var/www/kgp/FlaskApp/static/enrichr%s/%s_%s_%s.txt' % (tempid,tempid,namelist[i],'GO'),sep='\t')
					data = data[data[data.columns[1]] == abstr[j]]
					data.to_csv('/var/www/kgp/FlaskApp/static/enrichr%s/%s_%s_%s.txt' % (tempid,tempid,namelist[i],go_names[j]),index=False,header=['# GO','NS','enrichment','name','ratio_in_study','ratio_in_pop','p_uncorrected','depth','study_count','p_fdr_bh','study_items'],sep='\t')
					graphname = '/var/www/kgp/FlaskApp/static/enrichr%s/%s_%s_%s.png' % (tempid,tempid,namelist[i],go_names[j])
					filename = '/var/www/kgp/FlaskApp/static/enrichr%s/%s_%s_%s.txt' % (tempid,tempid,namelist[i],go_names[j])
					try:
						enrichr_grpah(tempid,go_names[j],graphname,filename)
					except:
						pass
			else:
				graphname = '/var/www/kgp/FlaskApp/static/enrichr%s/%s_%s_%s.png' % (tempid,tempid,namelist[i],db)
				filename = '/var/www/kgp/FlaskApp/static/enrichr%s/%s_%s_%s.txt' % (tempid,tempid,namelist[i],db)
				try:
					enrichr_grpah(tempid,db,graphname,filename)
				except:
					pass


	f = open ('/var/www/kgp/FlaskApp/static/enrichr%s/graphdone.txt' % tempid, 'w')
	f.write('This file is used to claim that all of the enrichment and corresponding visulization are done.')
	f.close()

	send_result(email, tempid)
	queue = open('/var/www/kgp/FlaskApp/data/queue.txt','a')
	queue.write(str(tempid)+' done\n')
	queue.close()
	

def send_result(email, tempid):
	if email != '':
		yag = yagmail.SMTP( user="dglinker.service@gmail.com", password="DGLinker!", host='smtp.gmail.com')
		subject = "Your result from DGLinker."
		contents = "Dear user,\n\n Please find attached a summary of your results. You can download the complete results at the following link: \n https://dglinker.rosalind.kcl.ac.uk/result/%s \nPlease note that the results will be deleted after 14 days. If you want to keep the results, please download them. \n\n Best wishes, \n The DGLinker Team" % tempid
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

		performance = []
		for line in model_performance:
			line = line.strip().split('\t')
			performance.append(line[1])
		


		return render_template('result.html',result=result, all_hits=all_hits, tempid=tempid, performance=performance)
	else:
		return render_template('redirect.html', tempid=tempid)


def pie_chart(tempid, info, i):
	graphname = '/var/www/kgp/FlaskApp/static/pie_chart_%s/%s_%s.png' % (tempid,tempid,i)
	plt.figure(num=1, clear=True)
	x = [i.replace('_', ' ') for i in info.keys()]
	sizes = [i/sum(info.values())*100 for i in info.values()]
	labels = ['{0} - {1:1.2f} %'.format(i,j) for i,j in zip(x, sizes)]
	colors = ['lightsalmon', 'peachpuff', 'lightblue', 'bisque', 'plum', 'khaki', 'turquoise','lightgreen'][:len(sizes)]
	patches,text = plt.pie(sizes,
						colors=colors,
						shadow = False,
						startangle =90,
						pctdistance = 0.6)
	patches, labels, dummy =  zip(*sorted(zip(patches, labels, sizes),
											key=lambda x: x[2],
											reverse=True))
	plt.axis('equal')
	plt.title(str(i))
	plt.legend(patches, labels, loc='upper right', bbox_to_anchor=(-0.1, 1.),fontsize=14)
	plt.savefig(graphname, dpi=100, transparent=True, bbox_inches='tight')
	plt.cla()
	plt.clf()
	plt.close()
	plt.close('all')



def enrichr_grpah(tempid,db,graphname,filename):
	if os.path.exists(filename):
		try:
			bar_data = pd.read_csv(filename,sep="\t",header=0)
			gene_name_list = list(bar_data["study_items"])
			pvalue = list(bar_data["p_fdr_bh"])
			rawp = list(bar_data["p_uncorrected"])
			# odds_radio = list(bar_data["Odds Ratio"])
			overlap = list(bar_data["ratio_in_study"])
			term = list(bar_data[bar_data.columns[0]])
			zipped = zip(gene_name_list,pvalue,rawp,overlap,term)
			sort_zipped = sorted(zipped,key=lambda x:x[1])
			result = zip(*sort_zipped)
			gene_name_list,pvalue,rawp,overlap,term = [list(x)[0:10][::-1] for x in result]
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
			ylabel_list = ['%s' % str(overlap[i]) for i in range(len(pvalue))]
			ax.set_yticklabels(ylabel_list)
			plt.xlabel("-log($\mathregular{P_{adj}}$)")
			plt.title(str(db).replace(' Illumina body map',''))
			plt.tight_layout()
			fig.savefig(graphname, dpi=400)
			plt.clf()
			plt.close('all')
		except:
			pass


@app.route('/download/<tempid>')
def downloadFile(tempid):
	path = '/var/www/kgp/FlaskApp/cache/zip%s.zip' % tempid
	if not os.path.exists(path):
		filelist1 = ['/var/www/kgp/FlaskApp/model_results/job_summary_%s.txt' % tempid, '/var/www/kgp/FlaskApp/model_results/results_%s.txt' % tempid, '/var/www/kgp/FlaskApp/model_results/graph_file%s.txt' % tempid, '/var/www/kgp/FlaskApp/model_results/model_performance_%s.txt' % tempid]
		dirpath1 = '/var/www/kgp/FlaskApp/static/enrichr%s/' % tempid
		dirpath2 = '/var/www/kgp/FlaskApp/static/graph_dir%s/' % tempid
		dirpath3 = '/var/www/kgp/FlaskApp/static/pie_chart_%s/' % tempid
		dirpath4 = '/var/www/kgp/FlaskApp/model_results/results_genomic_regions_%s/' % tempid
		filelist2 = search(dirpath1, '_')
		filelist3 = search(dirpath2, '.txt')
		filelist4 = search(dirpath2, '_')
		filelist5 = search(dirpath3, '.png')
		filelist6 = search(dirpath4, '.bed')

		for f in filelist1 + filelist2 + filelist3 + filelist4 + filelist5 + filelist6:
			addfile('/var/www/kgp/FlaskApp/cache/zip%s.zip' % tempid, f,-2)
		addfile('/var/www/kgp/FlaskApp/cache/zip%s.zip' % tempid, '/var/www/kgp/FlaskApp/README.TXT',-1)
	return send_file(path, as_attachment=True)

def addfile(zipfilename, filename,i):
	dir = "/".join([i for i in filename.strip().split("/")][:i])
	new_filename = "/".join([i for i in filename.strip().split("/")][i:])
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
		if os.path.exists('/var/www/kgp/FlaskApp/model_results/results_%s.txt' % tempid) or os.path.exists('/var/www/kgp/FlaskApp/cache/output_file%s.csv' % tempid) or os.path.exists('/var/www/kgp/FlaskApp/cache/error_file%s.txt' % tempid):
			return redirect(url_for('result',tempid=tempid))
		else:
			return render_template('none.html')
	else:
		return render_template('none.html')


@app.route('/downloads')
def downloads():
	return render_template('downloads.html')


if __name__ == '__main__':
	app.static_folder = 'static'
	app.run(debug=True)



