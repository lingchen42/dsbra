# sam_dsbra_grapher.py <database.p> [-indel]
# based off sam_dsbra_interface.py v0.1.1
# For use with mongo databases

import sys
import os
import ast
import csv
import cPickle
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import lines, gridspec
from scipy.stats import mannwhitneyu

VERBOSE = False

COLORS = ['#e1f5fe','#f0f4c3','#81d4fa','#e6ee9c','#dce775']

font = {'family' : 'sans-serif'}
matplotlib.rc('font', **font)

if __name__ == '__main__':
	VERBOSE = '-v' in sys.argv

	# Load mutational databases if they exists
	samples = [ ]
	summary_names = [ ]
	for arg in sys.argv:
		if '.txt' in arg[-4:]:
			summary_names.append(arg)
			
	break_index,margin,last_margin = sys.argv[sys.argv.index('-b')+1].split(',')
	break_index = int(break_index)
	margin = int(margin)
	last_margin = int(last_margin)
			
	if not summary_names:
		print "WARNING: Mutational spectrum database not found. Please try again."
		sys.exit()
		
	sample_count_list = [ ]
		
	for summary_name in summary_names:	
		#TODO: Implement multiple margins.
		#print "WARNING: Multiple margins not supported."		
		
		if os.path.exists(summary_name):
			with open(summary_name,'r') as sf:
				sample_count_list.append(int(sf.readline()[1:]))
			
				data = list(csv.DictReader(sf,delimiter='\t'))
				#formatting of all fields appropriately
				for x in data:
					x['insertion_seqs'] = ast.literal_eval(x['insertion_seqs'])
					x['deletion_lens'] = ast.literal_eval(x['deletion_lens'])
					x['num_transitions'] = int(x['num_transitions'])
					x['num_deletions'] = int(x['num_deletions'])
					x['num_insertions'] = int(x['num_insertions'])
					x['new_mismatch_bases'] = ast.literal_eval(x['new_mismatch_bases'])
					x['insertion_locs'] = ast.literal_eval(x['insertion_locs'])
					x['deletion_is_micro'] = ast.literal_eval(x['deletion_is_micro'])
					x['micro_seq'] = ast.literal_eval(x['micro_seq'])
					x['count'] = int(x['count'])
				
				samples.append(data)
				
				if not data:
					print "WARNING: Mutational spectrum database not found. Please try again."
					sys.exit()
					
	if '-typedist' in sys.argv:
		if len(samples) > 1:
			print summary_names
			print "WARNING: Graph -typedist does not support multiple mutational databases."
			
		ins_count = 0
		del_count = 0
		mismatch_count = 0
		compound_count = 0
		
		sample = samples[0]
		for x in sample:
			if x['num_deletions'] + x['num_insertions'] > 1:
				compound_count += x['count']
			elif x['num_deletions'] == 1:
				del_count += x['count']
			elif x['num_insertions'] == 1:
				ins_count += x['count']
			else:
				mismatch_count += x['count']
		
		wt_count = sample_count_list[0] - (ins_count + del_count + mismatch_count + compound_count)
				
		fig, ax = plt.subplots()
		sizes = [wt_count, ins_count,del_count,mismatch_count,compound_count]
		labels = ['WT','Insertion','Deletion','Mismatch','Compound']
		colors = COLORS[:5]
		ax.pie(sizes,labels=labels,autopct='%1.1f%%',colors=colors)
		fig.set_size_inches(5,5)
		#plt.tight_layout()
		
		if '-s' in sys.argv:
			plt.savefig(sys.argv[sys.argv.index('-s')+1]+'_typedist.png',dpi=100, transparent=True)
		else:
			plt.show()
					
	if '-indel' in sys.argv:
		print "-indel graph shows only non-compound mutations."
	
		if len(samples) > 1:
			print summary_names
			print "WARNING: Graph -indel does not support multiple mutational databases."
			
		sample = samples[0]
	
		indel_max_graph = 35
		indel_dist = [0]*((indel_max_graph*2)+1)
		micro_indel_dist = [0]*((indel_max_graph*2)+1)
	
		for x in sample:
			if x['num_deletions'] + x['num_insertions'] > 1:
				continue
			elif x['num_deletions'] == 1:
				del_len = x['deletion_lens'][0]
				if del_len <= indel_max_graph:
					indel_dist[indel_max_graph-del_len] += x['count']
					if x['deletion_is_micro'][0]:
						micro_indel_dist[indel_max_graph-del_len] += x['count']
						
			elif x['num_insertions'] == 1:
				ins_len = len(x['insertion_seqs'][0])
				if ins_len <= indel_max_graph:
					indel_dist[indel_max_graph+ins_len] += x['count']
					
		read_count = np.sum([int(x['count']) for x in sample]) + sample_count_list[0]
		
		del_pct = float(np.sum(indel_dist[:indel_max_graph]))/read_count
		micro_pct = float(np.sum(micro_indel_dist))/read_count
		ins_pct = float(np.sum(indel_dist[indel_max_graph:]))/read_count
		
		#Add in WT
		indel_dist[indel_max_graph] = sample_count_list[0]
		
		print '%s\t%s\t%s\t%s'%(summary_name,del_pct,micro_pct,ins_pct)
		
		#norm_indel_count = np.array(indel_list_count).astype(float)/sum(indel_list_count)
				
		#print sum(indel_list_count)
		#print sum(indel_dist)
		#print '%s'%sum(indel_dist)
		#print '%s'%len(read_indels)
		
		#sorted_indels = [x for (y,x) in sorted(zip(indel_list_count,indel_lists))]
		#sorted_count = sorted(indel_list_count)
		
		fig, ax = plt.subplots()
		width = 0.9
		ind = np.arange(len(indel_dist))-indel_max_graph
		ax.bar(ind+0.5,np.array(indel_dist).astype(float)*100/read_count,width,color='black')
		ax.set_xticks(ind[0::5]+width/2)
		ax.set_xticklabels(ind[0::5])
		ax.set_xlim([-indel_max_graph,indel_max_graph+1])
		ax.set_ylim([0,100])
		ax.set_xlabel('Indel Sum on Range %s to %s'%(-margin,margin))
		ax.set_ylabel('Percent of Reads')
		plt.text(10,10,'Deletions: %.2f%% \nMicro Dels: %.2f%% \nInsertions: %.2f%%'%(del_pct*100,micro_pct*100,ins_pct*100))
		plt.tight_layout()

		if '-s' in sys.argv:
			plt.savefig(sys.argv[sys.argv.index('-s')+1]+'_indel.png',dpi=100)
		else:
			plt.show()
		
	if '-microdel' in sys.argv:
		print "-microdel graph shows only non-compound mutations."
	
		if len(samples) > 1:
			print "WARNING: Graph -microdel does not support multiple mutational databases."
			
		sample = samples[0]
		
		max_del = int(sys.argv[sys.argv.index('-microdel')+1])
		del_length_list = np.zeros(max_del+1)
		micro_del_length_list = np.zeros(max_del+1)
		read_count = 0
	
		for x in sample:
			if x['num_deletions'] + x['num_insertions'] > 1:
				continue
			elif x['num_deletions'] == 1:
				del_len = x['deletion_lens'][0]
				if del_len <= indel_max_graph:
					indel_dist[indel_max_graph-del_len] += x['count']
					if x['deletion_is_micro']:
						micro_indel_dist[indel_max_graph-del_len] += x['count']
					
		print('Incomplete Data Analysis')
		sys.exit()
		
		for x in dsb_mut_spect:
			#skip if incomplete
			if len(x[1][0]) < 2*margin:
				continue
		
			read_count += len(x[0])
			total_indels = sum(x[1][0])
			if total_indels < 0:
				if abs(total_indels) >= max_del:
					del_length_list[max_del] += len(x[0])
				else:
					del_length_list[abs(total_indels)] += len(x[0])
					
				if x[1][2].count('') < len(x[1][2]):
					micro_del_length_list[abs(total_indels)] += len(x[0])
					
		del_length_list = np.array(del_length_list)/read_count
		micro_del_length_list = np.array(micro_del_length_list)/read_count
	
		fig, ax = plt.subplots()
		ind = np.arange(len(del_length_list)-1)
		width = 0.45

		rects1 = ax.bar(ind, del_length_list[1:]*100, width, color='black')
		rects2 = ax.bar(ind+width, micro_del_length_list[1:]*100, width, color='white', hatch='//')
		ax.set_ylabel('% of Reads')
		ax.set_xlabel('Deletion Length (bp)')
		ax.set_title('Deletions by bp Length')
		ax.set_xticks(ind+width)
		ax.set_xticklabels(range(1,len(del_length_list)))
		ax.legend((rects1[0],rects2[0]), ('Regular','Microhomology'))

		plt.tight_layout()
		if '-s' in sys.argv:
			plt.savefig(sys.argv[sys.argv.index('-s')+1]+'_microdel.png',dpi=100)
		else:
			plt.show()

		
	if '-beeswarm' in sys.argv:
		print "WARNING: Code does not distinguish compound mutants."
		
		w = 0.75 #width of box plot
		
		del_lists = [ ]
		X = [ ] # x scatter
		#multi = 0
		
		for sample in samples:
			del_list = [ ]
			for x in sample:
				#skip if compound mutant
				if x['num_deletions'] + x['num_insertions'] > 1:
						continue
				elif x['num_deletions'] == 1:
					del_list += [x['deletion_lens'][0]]*x['count']
						
			del_lists.append(del_list)
			X = np.concatenate((X, np.random.normal(len(del_lists),float(w)/15,size=len(del_list))))
		
		boxcolor = '#132029'
		meancolor = '#1e3140'
		plotcolor = '#48779e'
		
		fig = plt.figure(1)
		gs1 = None
		if '-p' in sys.argv:
			gs1 = gridspec.GridSpec(5,1)
			fig.set_size_inches(2+len(del_lists)*1.4,5+0.25*len(del_lists))
		else:
			gs1 = gridspec.GridSpec(1,1)
			fig.set_size_inches(2+len(del_lists)*1.4,4+0.25*len(del_lists))
			
		ax = plt.subplot(gs1[:-1,0])
			
		#fig.set_size_inches(2+len(del_lists)*1.4,4+0.25*len(del_lists))
		
		#TODO: generalize this
		box_labels = [ ]
		for p in summary_names:
			box_labels.append(p[:-11])
		
		ax.boxplot(del_lists,widths=0.6*w,labels=box_labels,showmeans=True,showcaps=False,showfliers=False, meanprops = {'marker':'*','markerfacecolor':meancolor}, boxprops = {'linewidth':2.5,'alpha':0.5,'color':boxcolor}, whiskerprops = {'alpha':0}, medianprops = {'color':boxcolor,'alpha':0.75,'linewidth':3.5}) # 3x stdev on either side
		#X = np.random.normal(1,float(w)/15,size=len(del_list))
		ax.scatter(X,[mut for dl in del_lists for mut in dl],alpha=0.5,color=plotcolor,marker='.')
		ax.set_ylabel('Deletion Size (bp)')
		
		#change this
		max_delsize = 60
		
		ax.set_ylim([0,max_delsize])
		#plot mean and average
		for del_list in del_lists:
			plt.text(del_lists.index(del_list)+0.7,np.median(del_list),'%s'%(int(np.median(del_list))),verticalalignment='bottom',\
				horizontalalignment='right', color=boxcolor, fontweight='bold', fontsize=12)
			plt.text(del_lists.index(del_list)+1.45,np.mean(del_list)+0.05*max_delsize,'%.2f'%(np.mean(del_list)),verticalalignment='bottom',\
				horizontalalignment='right', color=meancolor, fontweight='bold', fontsize=12)

		#check mann-whitney and draw statistically significant lines
		if '-p' in sys.argv:
			gs1 = gridspec.GridSpec(5,1)
			gs1.update(wspace=0.025,hspace=0.05)
			ax2 = plt.subplot(gs1[-1,0],axisbg=(1,1,1,0))
			ax2.set_axis_off()
			ax2.set_xlim(ax.get_xlim())
			sig_values = 0
			for i in range(len(del_lists)):
				for j in range(i,len(del_lists)):
					p_val = mannwhitneyu(del_lists[i],del_lists[j])[1]
					if p_val < 0.05:
						#ax.plot([i+1,j+1],[75-sig_values*2.5,75-sig_values*2.5],lw=2,color='black')
						ax2.plot([i+1,j+1],[sig_values*3,sig_values*3],lw=2,color='black')
						#ax2.text(i+0.95,98-sig_values*3,"*",ha='right',va='center')
						sig_values += 1
				
			ax2.set_ylim([-2,(sig_values+1)*3])
		#plt.tight_layout()		

		if '-s' in sys.argv:
			plt.savefig(sys.argv[sys.argv.index('-s')+1]+'_delscatter.png',dpi=100)
		else:
			plt.show()
			
	if '-barswarm' in sys.argv:
		#if len(dsb_mut_spect_list) > 1:
			#print summary_names
		#	print "WARNING: Graph -indel does not support multiple mutational databases."
			
		dsb_mut_spect = dsb_mut_spect_list[0]
		
		boxcolor = '#132029'
		inscolor = '#1e3140'
		delcolor = '#48779e'
	
		read_indels = { }
		read_micro = { }
	
		for x in dsb_mut_spect:
			for y in x[0]:
				#skip incomplete coverage
				if len(x[1][0]) < 2*margin:
					#print 'skip %s'%read_indels[read]
					continue
				read_indels[y] = x[1][0]
				read_micro[y] = x[1][2]
	
		dels = 0
		ins = 0
		fs = 0
		unique = 0
		micro_dels = 0
		read_count = 0
		
		indel_lists = [ ]
		indel_list_count = [ ]
		
		indel_max_graph = 35
		
		indel_dist = [0]*((indel_max_graph*2)+1)
		
		for read in read_indels:
			#check for microhomology
			if read_micro[read].count('') < len(read_micro[read]):
				micro_dels += 1
		
			indel_sum = sum(read_indels[read])
			
			if not read_indels[read] in indel_lists:
				indel_lists.append(read_indels[read])
				indel_list_count.append(1)
				unique += 1
			else:
				indel_list_count[indel_lists.index(read_indels[read])] += 1
			
			if indel_sum < 0:
				dels += 1
			elif indel_sum > 0:
				ins += 1
				
			if abs(indel_sum) < indel_max_graph:
				indel_dist[indel_sum+indel_max_graph] += 1
				
			if indel_sum%3 != 0:
				fs += 1
			
		del_pct = float(dels)/sum(indel_list_count)
		micro_pct = float(micro_dels)/sum(indel_list_count)
		ins_pct = float(ins)/sum(indel_list_count)
		total_pct = (float(dels)+float(ins))/sum(indel_list_count)
		fs_pct = float(fs)/sum(indel_list_count)
		print '%s\t%s\t%s\t%s\t%s\t%s\t%s'%(summary_name,del_pct,micro_pct,ins_pct,total_pct,fs_pct,float(unique)/sum(indel_list_count))
		
		norm_indel_count = np.array(indel_list_count).astype(float)/sum(indel_list_count)
				
		#print sum(indel_list_count)
		#print sum(indel_dist)
		#print '%s'%sum(indel_dist)
		#print '%s'%len(read_indels)
		
		#sorted_indels = [x for (y,x) in sorted(zip(indel_list_count,indel_lists))]
		#sorted_count = sorted(indel_list_count)
		
		fig, axes = plt.subplots(ncols=2, sharey=True)
		ind = np.arange(1,indel_max_graph+1)
		axes[0].barh(ind,np.array(indel_dist[indel_max_graph+1:]).astype(float)*100/ins,color=inscolor,linewidth=0,height=1.0)
		axes[0].invert_xaxis()
		#ax.set_xticks(ind[0::5]+width/2)
		#ax.set_xticklabels(ind[0::5])
		#ax.set_xlim([-indel_max_graph,indel_max_graph+1])
		#ax.set_ylim([0,100])
		#ax.set_xlabel('Indel Sum on Range %s to %s'%(-margin,margin))
		#ax.set_ylabel('Percent of Reads')
		#plt.text(10,10,'Deletions: %.2f%% \nMicro Dels: %.2f%% \nInsertions: %.2f%% \nFrameshift: %.2f%%'%(del_pct*100,micro_pct*100,ins_pct*100,fs_pct*100))
		axes[1].barh(ind,np.array(indel_dist[:indel_max_graph][::-1]).astype(float)*100/dels,height=1.0,linewidth=0,color=delcolor)
		fig.subplots_adjust(wspace=0)
		axes[0].set_xlim(axes[1].get_xlim()[::-1])
		#plt.tight_layout()

		if '-s' in sys.argv:
			plt.savefig(sys.argv[sys.argv.index('-s')+1]+'_indel.png',dpi=100)
		else:
			plt.show()
		
	if '-micro' in sys.argv:
		# Use '-micro <max_micro>'
		# Make stacked bar chart 0-1,2,3,4+ microhomology
		# TODO: handle multiple microhomologies?
		print "WARNING: graph 'micro' unable to graph compound mutations."

		max_micro = -1
		
		try:
			max_micro = int(sys.argv[sys.argv.index('-micro')+1])
		except:
			print "Max microhomology must be specified after '-micro'"
			sys.exit()
			
		micro_lists = [ ]
			
		sample = samples[0]
		
		for sample in samples:
			micro_list = np.zeros(max_micro+1)
			read_count = 0
		
			#Count number of deletions of >2 bp in length
			for x in sample:
				if x['num_deletions'] + x['num_insertions'] > 1:
					continue
				elif x['num_deletions'] == 1:
					del_len = x['deletion_lens'][0]
					if del_len >= 2:
						read_count += x['count']
						if x['deletion_is_micro'][0]:
							micro_list[max(1,min(max_micro,len(x['micro_seq']),del_len))] += x['count']
							#if len(x['micro_seq']) > del_len:
							#	print "WARNING: Microhomology identified greater than deletion size."
							#	sys.exit()
						#No microhomology		
						else:
							micro_list[0] += x['count']
								
			micro_lists.append(np.array(micro_list).astype(float)*100/read_count)

		stacked_bars = [ ]

		for i in range(max_micro):
			m_list = [ ]
			for j in micro_lists:
				m_list.append(j[i])
			stacked_bars.append(m_list)

		#graphing
		fig,ax = plt.subplots()
		fig.set_size_inches(2+len(micro_lists)*1,6)
		width = 0.9

		ax.set_ylim([0,100])
		cs = [ '#3a828c', '#A1CDA8', '#B5DFCA', '#C5E7E2', '#AD9BAA' ]
		last_height = np.array([float(0) for x in range(len(micro_lists))]) #necessary for stacked bar

		bar_charts = [ ]
		ind = np.arange(len(micro_lists)) + (1-width)/2

		for x in range(max_micro):
			bar_charts.append(ax.bar(ind,stacked_bars[x],width,bottom=last_height, color = cs[x%len(cs)], edgecolor='black', lw=2))
			last_height += stacked_bars[x]
			print last_height

		chart_legends = ['0-1 bp']+['%s bp'%x for x in range(2,max_micro)]+['%s+ bp'%max_micro]
		ax.legend(bar_charts[::-1],chart_legends[::-1],loc=4)
		
		ax.set_xticks(ind+width/2)
		ax.set_xticklabels(summary_names,rotation=45, ha='right')
		ax.set_ylabel('Percent of Reads with Microhomology')
		plt.tight_layout()

		if '-s' in sys.argv:
			plt.savefig(sys.argv[sys.argv.index('-s')+1]+'_microhomology.png',dpi=100)
		else:
			plt.show()

