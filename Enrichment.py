import pandas as pd
from goatools.obo_parser import GODag
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS
from goatools.go_enrichment import GOEnrichmentRecord
from goatools.anno.gaf_reader import GafReader
from genes_ncbi_9606_proteincoding import GENEID2NT as GeneID2nt_human
import json, sys, os
from goatools.anno.genetogo_reader import Gene2GoReader
from collections import defaultdict, Counter
import collections as cx
from goatools.rpt.goea_nt_xfrm import MgrNtGOEAs
from goatools.rpt.prtfmt import PrtFmt


class Enricher():
	def __init__(self, obo_path="godata/go-basic.obo", gene2go_path='godata/gene2go', symbol2id_path="godata/ncbi_symbol_to_id.txt", go_propagate_counts = True, backgrounds_path="datasets", backgrounds_files=[], log=sys.stdout):
		self.backgrounds = {}
		self.log = log
		obodag = GODag(obo_path)

		# Read NCBI's gene2go. Store annotations in a list of namedtuples
		# gene2go is from ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz
		objanno = Gene2GoReader(gene2go_path, taxids=[9606])

		# Get associations for each branch of the GO DAG (BP, MF, CC)
		ns2assoc = objanno.get_ns2assc()

		with open(symbol2id_path) as f:
			self.symbol2id = json.load(f)

		###
		# many properties of this object are borrowed for other enrichments
		# be very careful if modifying, may need more methods on Enricher if so
		# especially if changing the correction method. Only the statsmodels part of the API is implemented so these should work:
		# 'bonferroni',	 #  0) Bonferroni one-step correction
		#	 'sidak',		  #  1) Sidak one-step correction
		#	 'holm-sidak',	 #  2) Holm-Sidak step-down method using Sidak adjustments
		#	 'holm',		   #  3) Holm step-down method using Bonferroni adjustments
		#	 'simes-hochberg', #  4) Simes-Hochberg step-up method  (independent)
		#	 'hommel',		 #  5) Hommel closed method based on Simes tests (non-negative)
		#	 'fdr_bh',		 #  6) FDR Benjamini/Hochberg  (non-negative)
		#	 'fdr_by',		 #  7) FDR Benjamini/Yekutieli (negative)
		#	 'fdr_tsbh',	   #  8) FDR 2-stage Benjamini-Hochberg (non-negative)
		#	 'fdr_tsbky',	  #  9) FDR 2-stage Benjamini-Krieger-Yekutieli (non-negative)
		#	 'fdr_gbs',		# 10) FDR adaptive Gavrilov-Benjamini-Sarkar
		###
		self.goeaobj = GOEnrichmentStudyNS(
				GeneID2nt_human.keys(), # List of protein-coding genes
				ns2assoc, # geneid/GO associations
				obodag, # Ontologies
				propagate_counts = go_propagate_counts,
				alpha = 0.05, # default significance cut-off
				methods = ['fdr_bh'], # defult multipletest correction method
				log=log) 
		
		self.pval_obj = self.goeaobj.ns2objgoea['BP'].pval_obj
		self.methods = self.goeaobj.ns2objgoea['BP'].methods
		self.alpha = self.goeaobj.ns2objgoea['BP'].alpha
		self.objprtres = self.goeaobj.ns2objgoea['BP'].objprtres
		self.writer = self.goeaobj.ns2objgoea['BP']

		self._run_multitest = {
			'statsmodels':self._run_multitest_statsmodels}
		
		self.make_backgrounds(backgrounds_path, backgrounds_files, log)
	
	def enrichment(self, genes, return_filter="all", dataset='GO'):
		"""
		genes = list of gene symbols
		return_filter = "all" - return all enrichment results
						"significant" - return only significant results after correction
		dataset = which dataset to test enrichment in, must be "GO" or a key in self.backgrounds

		returns:
			list of GOEnrichmentRecord (n.b. save to tsv file using Enricher.to_tsv)
		"""
		log = self.log
		geneids_study = set([self.symbol2id[x] for x in genes if x in self.symbol2id])
		print(f"Got {len(genes)} and mapped {len(geneids_study)}")
		if dataset == "GO":
			goea_results_all = self.goeaobj.run_study(geneids_study)
			goea_results_sig = [r for r in goea_results_all if r.p_fdr_bh < 0.05]
			if return_filter == 'significant':
				return goea_results_sig
			else:
				return goea_results_all
		elif dataset in self.backgrounds:
			#results = self.get_pval_uncorr(geneids_study, dataset)
			results = self.get_pval_uncorr(set(genes), dataset, log)
			if not results:
				return []
			self._run_multitest_corr(results, self.methods, self.alpha, genes)
			results.sort(key=lambda r: [r.enrichment, r.NS, r.p_uncorrected])
			results_sig = [r for r in results if r.p_fdr_bh < 0.05]
			if return_filter == 'significant':
				return results_sig
			else:
				return results
		else:
			if log is not None:
				log.write(f"Dataset not found: {dataset}")
			return []
	
	def get_pval_uncorr(self, study, dataset, log=sys.stdout):
		"""Calculate the uncorrected pvalues for study items.
		Adapted from GOAtools https://github.com/tanghaibao/goatools"""

		if dataset not in self.backgrounds:
			if log is not None:
				log.write(f"Dataset not found: {dataset}")
			return []
		
		results = []
		study_in_pop = self.backgrounds[dataset]['pop'].intersection(study)
		# " 99%	378 of	382 study items found in population"

		go2studyitems = self.get_terms(study_in_pop, self.backgrounds[dataset]['item2terms'])
		pop_n, study_n = self.backgrounds[dataset]['pop_n'], len(study_in_pop)
		allterms = set(go2studyitems).union(set(self.backgrounds[dataset]['term2items']))
		if log is not None:
			# Some study genes may not have been found in the population. Report from orig
			study_n_orig = len(study)
			perc = 100.0*study_n/study_n_orig if study_n_orig != 0 else 0.0
			log.write("{R:3.0f}% {N:>6,} of {M:>6,} study items found in population({P})\n".format(
				N=study_n, M=study_n_orig, P=pop_n, R=perc))
			if study_n:
				log.write("Calculating {N:,} uncorrected p-values using {PFNC}\n".format(
					N=len(allterms), PFNC=self.pval_obj.name))
		# If no study genes were found in the population, return empty GOEA results
		if not study_n:
			return []
		calc_pvalue = self.pval_obj.calc_pvalue

		for goid in allterms:
			study_items = go2studyitems.get(goid, set())
			study_count = len(study_items)
			pop_items = self.backgrounds[dataset]['term2items'].get(goid, set())
			pop_count = len(pop_items)

			one_record = GOEnrichmentRecord(
				goid,
				p_uncorrected=calc_pvalue(study_count, study_n, pop_count, pop_n),
				study_items=study_items,
				pop_items=pop_items,
				ratio_in_study=(study_count, study_n),
				ratio_in_pop=(pop_count, pop_n))

			results.append(one_record)

		return results
	
	def _run_multitest_corr(self, results, usrmethod_flds, alpha, study):
		"""Do multiple-test corrections on uncorrected pvalues.
		Adapted from GOAtools https://github.com/tanghaibao/goatools"""
		assert 0 < alpha < 1, "Test-wise alpha must fall between (0, 1)"
		pvals = [r.p_uncorrected for r in results]
		ntobj = cx.namedtuple("ntobj", "results pvals alpha nt_method study")
		for nt_method in usrmethod_flds:
			ntmt = ntobj(results, pvals, alpha, nt_method, study)
			self._run_multitest[nt_method.source](ntmt)
	
	def _run_multitest_statsmodels(self, ntmt):
		"""Use multitest mthods that have been implemented in statsmodels.
		Adapted from GOAtools https://github.com/tanghaibao/goatools"""
		# Only load statsmodels if it is used
		# print("running correction")
		multipletests = self.methods.get_statsmodels_multipletests()
		results = multipletests(ntmt.pvals, ntmt.alpha, ntmt.nt_method.method)
		pvals_corrected = results[1] # reject_lst, pvals_corrected, alphacSidak, alphacBonf
		# print(f"have {len(pvals_corrected)} corrected values")
		self._update_pvalcorr(ntmt, pvals_corrected)
	
	@staticmethod
	def _update_pvalcorr(ntmt, corrected_pvals):
		"""Add data members to store multiple test corrections."""
		if corrected_pvals is None:
			return
		for rec, val in zip(ntmt.results, corrected_pvals):
			rec.set_corrected_pval(ntmt.nt_method, val)
	
	

	
	def get_terms(self, geneset, assoc):
		"""Get the terms in the study group
		Adapted from GOAtools https://github.com/tanghaibao/goatools"""

		term2itemids = defaultdict(set)
		genes = [g for g in geneset if g in assoc]
		for gene in genes:
			for goid in assoc[gene]:
				term2itemids[goid].add(gene)
		return term2itemids

	def make_backgrounds(self, path, files, log=sys.stdout ):

		for fname in files:
			if log is not None:
				log.write(f"Loading background for {fname}")
			bgname = fname.replace('.csv', '')
			df = pd.read_csv(os.path.join(path, fname), header=None)
			df.columns = ['source', 'source type', 'edge type', 'target type', 'target']
			#bg = [self.symbol2id[x] if x in self.symbol2id else None for x in df['source'] ]
			bg = df['source']
			df['bg_gene_id'] = bg
			df = df.dropna(subset=['bg_gene_id'])
			gr = df.groupby('target')
			term2items = {}
			for name, group in gr:
				term2items[name] = set(group['source'])
			
			gr = df.groupby('source')
			item2terms = {}
			for name, group in gr:
				item2terms[name] = set(group['target'])
			bgobj = {}
			bgobj['term2items'] = term2items
			bgobj['pop'] = set(bg)
			bgobj['pop_n'] = len(bgobj['pop'])
			bgobj['item2terms'] = item2terms
			self.backgrounds[bgname] = bgobj

	def to_tsv(self, filename, results):
		self.writer.wr_tsv(filename, results)

		

if __name__ == '__main__':
	files = [
			'ClinVar v2020-04.csv',
			'DisGeNet v7.0.csv',	
			'ArrayExpress_Atlas(experiment_E-MTAB-513) v2014-07.csv',
		]
	enr = Enricher(backgrounds_files=files)
	
	genes_study =  ['IL6', 'CAT', 'AGT', 'VEGFA', 'IFNG', 'IL1B', 'MAP3K5', 'MMP9', 'HMOX1', 'PPARA', 'AIFM1', 'CCL2', 'IGF1', 'BAX', 'MMP2', 'CRP', 'PRKCD', 'APOE', 'ACSL1', 'PDK4', 'MTHFR', 'GSTM1', 'MAP1LC3B', 'AHR', 'ADIPOQ', 'FASN', 'ALB', 'ACTB', 'SPP1', 'ABCB1', 'PLAT', 'IFNA2', 'CXCL10', 'BRD4', 'NR1H4', 'CP', 'CSF3', 'NOTCH1', 'ENO1', 'PSEN1', 'HSPA5', 'LCN2', 'NFE2L2', 'GSTA4', 'PPARG', 'MAPT', 'CKB', 'TUBA1A', 'HMGB1', 'DLAT', 'ABCC2', 'RGN', 'CYP2C9', 'DHFR', 'CYP2E1', 'FGA', 'HAVCR1', 'GSN', 'TALDO1', 'EPHX1', 'NME2', 'GPX1', 'ATG5', 'MAT1A', 'RIPK3', 'TTR', 'APOA1', 'CANX', 'PNP', 'CCR2', 'VWF', 'PRKDC', 'TUBB4B', 'HADHA', 'CTSE', 'HACL1', 'CPS1', 'CFH', 'ME1', 'ARG1', 'ERN1', 'GSTM3', 'ALDOB', 'SORD', 'IRAK1', 'ANXA2', 'ANXA6', 'APOH', 'SOD3', 'SESN2', 'PC', 'MTPN', 'PYGL', 'SLC22A8', 'PCYT1A', 'NOX4', 'IL4', 'GSTM4', 'FADS2', 'CCT3']
	res = enr.enrichment(genes_study, return_filter='significant', dataset='GO')
	enr.to_tsv('tests_enrichment/test0.tsv', res)
	print(res[:10])
	res = enr.enrichment(genes_study, return_filter='significant', dataset='DisGeNet v7.0')
	print(res[:10])
	res = enr.enrichment(genes_study, return_filter='significant', dataset='ArrayExpress_Atlas(experiment_E-MTAB-513) v2014-07')
	enr.to_tsv('tests_enrichment/test0.1.tsv', res)
	print(res[:10])

	#test a list that should be enriched in disgenet
	genes_study = ["GRK2","COL8A1","SLC31A1","EYA1","FMR1","FOLR1","NOTCH1","SALL1","SIX1","KLF4","DLL3",]
	res = enr.enrichment(genes_study, return_filter='significant', dataset='DisGeNet v7.0')
	print("should be enriched for C0000768")
	print(res[:10])
	enr.to_tsv('tests_enrichment/test1.tsv', res)

	#get a list that must have one member of every disease
	df = pd.read_csv('datasets/DisGeNet v7.0.csv', header=None)
	df.columns = ['source', 'source type', 'edge type', 'target type', 'target']
	smp = df.drop_duplicates(subset=['target'])
	genes_study = smp['source'].tolist()
	res = enr.enrichment(genes_study, return_filter='all', dataset='DisGeNet v7.0')
	print("should have every condition")
	print(res[:10])
	enr.to_tsv('tests_enrichment/test2.tsv', res)