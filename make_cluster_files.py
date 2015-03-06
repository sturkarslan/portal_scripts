#!/usr/bin/env python

import os, re, sys, string

orgs = ['Phatr']

class Cluster:
	def __init__(self,index=0,cluster_label=[],members=[],motif_pngs=[],motif_Evals=[],plots={},metrics={}):
		self.index = index
		self.members = members
		self.cluster_label = cluster_label 
		self.motif_pngs = motif_pngs
		self.motif_Evals = motif_Evals
		self.plots = plots
		self.metrics = metrics
	def __str__(self):
		delim = '\t'
		subdelim = ';'
		out = []
		out.append(str(self.index))
		out.append(str(self.cluster_label))
		out.append(string.join(self.members,subdelim))
#		print self.motif_pngs
		if self.motif_pngs != []: out.append(string.join(self.motif_pngs,subdelim))
		if self.motif_Evals == []: out.append('NA')
		else: out.append(string.join([str(e) for e in self.motif_Evals],subdelim))
		for t,p in self.plots.items(): out.append(p)
		for t,m in self.metrics.items(): out.append(str(m))
		return string.join(out,delim)
	def header(self):
		delim = '\t'
		out = []
		out.append('index')
		out.append('cluster_label')
		out.append('members')
		if self.motif_pngs != []: out.append('motif_pngs')
		out.append('motif_Evals')
		for t in self.plots: out.append(t)
		for t in self.metrics: out.append(t)
		return string.join(out,delim)

def readMEME(memefile):
	res = {
		'widths' : [],
		'llrs'   : [],
		'sites'  : [],
		'Evals'  : []
	}
	re_stats = re.compile('MOTIF\s+(\w)\s+MEME\s+width =\s+([0-9]+)\s+sites =\s+([0-9]+)\s+llr =\s+([0-9+-.]+)\s+E-value =\s+([0-9-+e.]+)')
	for m in re_stats.finditer(open(memefile).read()):
		name,width,sites,llr,Eval = m.groups()
		res['widths'].append(int(width))
		res['llrs'].append(float(llr))
		res['sites'].append(int(sites))
		res['Evals'].append(float(Eval))
#	print res
	return res

def make_cmonkey_cluster_table(dir,base='data-files',sub='cmonkey'):
	print 'make_cmonkey_cluster_table'
	workdir = '%s/%s' %(dir,sub)
	cls = {}
	header = None

	cluster_members = {}
	for l in open('%s/cluster.members.genes.txt' %workdir):
		fields = l.strip().split()
		if org == 'Phatr':
			geneString = 'PHATRDRAFT_'
		cluster_members[fields[0]] = [geneString + x for x in fields[1:]] 

	for l in open('%s/cluster.summary.tsv' %workdir):
		fields = l.strip().split()
		if not header:
			header = fields
			continue
		index = int(fields[ header.index('k')+1 ])
		cluster_label = '%s_bicluster_%04d' %(org,index)
		resid = float(fields[ header.index('resid')+1 ])
		members = cluster_members[str(index)]
		pngdir = '%s/cmonkey_motifs' %workdir
		motif_Evals = []
		# the R tsv output from cmonkey can leave empty fields, confusing split() indexing...
		if len(fields) == 7 and fields[5:7] == ['Inf','Inf']:
			motif_Evals = [float('inf'),float('inf')]
		else:
			for field in header:
				if re.search('e.value',field):
#					print index, field, fields, header.index(field)
					motif_Evals.append(float(fields[ header.index(field)+1 ]))

			motif_pngs = []
			for i in range(20):
				motpng = 'motif_%04d_%i.png' %(index,i)
				if motpng in os.listdir(pngdir): motif_pngs.append('%s/%s/%s' %(base,pngdir,motpng))

		plots = {
			'plot_exp' : '%s/%s/svgs/cluster%04d.svgz' %(base,workdir,index)
		}

		cls[index] = Cluster(index,cluster_label,members,motif_pngs,motif_Evals,plots=plots,metrics={'residue':resid})

	outf = 'clusters.%s.%s.tsv' %(dir,sub)
	outf = open(outf,'w')
	header = True
	for index,cl in cls.items():
		if header:
			outf.write(cl.header()+'\n')
			header = False
		outf.write(str(cl)+'\n')
	outf.close()

def make_hc_cluster_table(dir,base='data-files',sub='hierarchical_clustering'):
	print 'make_hc_cluster_table'
	workdir = '%s/%s' %(dir,sub)
	cls = {}
	header = None
	for l in open('%s/cluster_stats.tsv' %workdir):
		fields = l.strip().split()
		if not header:
			header = fields
			continue
		phatrString = 'PHATRDRAFT_'
                #cluster_members[fields[0]] = [phatrString + x for x in fields[1:]]
		
		index = int( fields[ header.index('clust') ] )
		if org == 'Phatr':
			geneString = 'PHATRDRAFT_'
		members = [phatrString + x for x in fields[ header.index('ids') ].split(';')]
		#members = fields[ header.index('ids') ].split(';')
		plots = {
			'plot_exp'    : '%s/%s/plots/cluster_%i.png' %(base,workdir,index),
			'plot_conds'  : '%s/%s/condplots/condplot_%i.png' %(base,workdir,index),
			'plot_dendro' : '%s/%s/dendro/pvdend.%i.png' %(base,workdir,index)
		}

		memedir = '%s/motifs/%i' %(workdir,index)
		memefile = '%s/meme.txt' %(memedir)
		meme = readMEME(memefile)
		motif_Evals = meme['Evals']
		motif_pngs = []
		for i in range(20):
			motpng = 'logo%i.png' %i
			if motpng in os.listdir(memedir): motif_pngs.append('%s/%s/%s' %(base,memedir,motpng))
		cluster_label = '%s_hcluster_%04d' %(org,index)

#		print motif_pngs
		cls[index] = Cluster(index,cluster_label,members,motif_pngs,motif_Evals,plots,metrics={})

	outf = 'clusters.%s.%s.tsv' %(dir,sub)
	outf = open(outf,'w')
	header = True
	for index,cl in cls.items():
		if header:
			outf.write(cl.header()+'\n')
			header = False
		outf.write(str(cl)+'\n')
	outf.close()

for org in orgs:
	print org
	if not org in os.listdir('.'):
		print '%s not found, skipping' %org
		continue
	make_hc_cluster_table(org)
	#make_cmonkey_cluster_table(org)
