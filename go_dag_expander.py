from collections import defaultdict
import pandas as pd
def read_obo(fname, print_errors=False):
	all_terms = {}

	term_start = False
	typedef_start = False
	term_id = "NOT SET"
	term = {}
	with open(fname) as f:
		for line in f:
			line = line.rstrip()
			if line.lower() == '[term]':
				term_start = True
				term_id = "NOT SET"
				term = {}
				continue
			elif line.lower() == '[typedef]':
				typedef_start = True
			elif line == '':
				if term != {}:
					all_terms[term_id] = term
				term_start = False
				typedef_start = False
				continue
			elif typedef_start:
				continue
			elif term_start:
				l = line.split(': ')
				if len(l) > 2:
					pos = line.find(': ')
					key = line[:pos]
					value = line[(pos+2):]
					if key != 'def' and key != 'comment' and print_errors:
						print("error in line", line)
				else:
					key = l[0]
					value = l[1]
				if key not in term:
					term[key] = []
				term[key].append(value)
				if key == 'id':
					term_id = value
			else:
				if print_errors:
					print("line not handled", line)
	
	return all_terms

def obodict_to_rels(all_terms):
	missing = set()
	isa_ch2pt = defaultdict(set)
	isa_pt2ch = defaultdict(set)
	for term in all_terms:
		if 'is_a' in all_terms[term]:
			for target in all_terms[term]['is_a']:
				tgt_name = target.split(' ! ')[0]
				if tgt_name not in all_terms:
					missing.add(tgt_name)
				else:
					isa_ch2pt[term].add(tgt_name)
					isa_pt2ch[tgt_name].add(term)
	return isa_ch2pt, isa_pt2ch

def get_all_parents(ch2pt):
	maps = {}
	for term in ch2pt:
		maps[term] = get_all_parents_for_term(term, ch2pt)
	return maps

def get_all_parents_from_obo(go_fname):
	terms = read_obo(go_fname)
	ch2pt, pt2ch = obodict_to_rels(terms)
	ch2all = get_all_parents(ch2pt)
	return ch2all

def get_all_parents_for_term(term, ch2pt):
	ret = set([term])
	if term not in ch2pt or ch2pt[term] == []:
		return ret
	else:
		for parent_term in ch2pt[term]:
			ret.add(parent_term)
			ret.update(get_all_parents_for_term(parent_term, ch2pt))
		return ret

def expand_rels(df, all_parents):
	"""
	rels - pandas dataframe in format ready for edgeprediction library
	"""
	rels = df.drop_duplicates(subset=['source', 'target'])
	new_rows = []
	not_mapped = set()
	for index, row in rels.iterrows():
		if row['target'] in all_parents:
			pts = all_parents[row['target']]
			for pt in pts:
				r = row.to_dict()
				r['target'] = pt
				new_rows.append(r)
		else:
			not_mapped.add(row['target'])
			continue

	df2 = pd.DataFrame(new_rows)
	print(df2.shape[0], 'relations added by expansion')
	a = rels.shape[0]
	b = df2.shape[0]
	combined = pd.concat([rels, df2], ignore_index=True, sort=False)
	combined.drop_duplicates(subset=['source', 'target'], inplace=True)
	c = combined.shape[0]
	print(f"had {a}, added {b}, combine and dedup: {c} factor {c/a}. Terms not mapped: {len(not_mapped)}.")
	return combined


def process_rels(rels_fname, go_fname):
	df = pd.read_csv(rels_fname, header=None)
	df.columns = ['source', 'source type', 'edge type', 'target type', 'target']

	ch2all = get_all_parents_from_obo(go_fname)
	
	expanded = expand_rels(df, ch2all)
	return expanded


if __name__ == '__main__':
	terms = read_obo('godata/go-basic.obo')
	ch2pt, pt2ch = obodict_to_rels(terms)
	pts = get_all_parents_for_term('GO:0060362', ch2pt)
	print(pts)

	print("** small test")
	expanded = process_rels('datasets/go_sample.csv', 'godata/go-basic.obo')
	expanded.to_csv('go_tests/expanded.csv')
	#to save in same format used in DGL server
	expanded.to_csv('go_tests/expanded_dgl.csv', header=False, index=False)

	print("** 1k test")
	expanded = process_rels('datasets/go_sample_1k.csv', 'godata/go-basic.obo')
	print("** whole database test")
	expanded = process_rels('datasets/Gene Ontology v2020-11-17.csv', 'godata/go-basic.obo')
