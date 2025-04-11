#!/usr/bin/python
import os,sys
import json
import argparse

parser = argparse.ArgumentParser(description="Parse fastp output")
parser.add_argument('--fastp-json','-i',help='Fastp json file',required=True)
parser.add_argument('--outdir','-o',help='Output directory',required=True)
parser.add_argument('--samplename','-s',help='sample name',required=True)
argv = parser.parse_args()

json_file = argv.fastp_json
sample = argv.samplename
outdir = argv.outdir

prefix = sample

def qs2error(qual_score):
	return 10**(-qual_score/10)

def mean(list_data):
	if len(list_data) == 0:
		return 0
	return sum(list_data)/float(len(list_data))

data = json.load(open(json_file))


raw_len1 =data['read1_before_filtering']['total_cycles']
raw_len2 =data['read2_before_filtering']['total_cycles']
clean_len1 =data['read1_after_filtering']['total_cycles']
clean_len2 =data['read2_after_filtering']['total_cycles']

## raw left read
for i in range(raw_len1):
	qs = data['read1_before_filtering']['quality_curves']['mean'][i]
	gc_content = data['read1_before_filtering']['content_curves']
## raw right read
for i in range(raw_len2):
	qs = data['read2_before_filtering']['quality_curves']['mean'][i]
	gc_content = data['read2_before_filtering']['content_curves']
## clean left read
for i in range(clean_len1):
	qs = data['read1_after_filtering']['quality_curves']['mean'][i]
	gc_content = data['read1_after_filtering']['content_curves']
## clean right read
for i in range(clean_len2):
	qs = data['read2_after_filtering']['quality_curves']['mean'][i]
	gc_content = data['read2_after_filtering']['content_curves']

q20_raw1 = float(data['read1_before_filtering']['q20_bases'])/(data['read1_before_filtering']['total_bases'])*100
q20_raw2 = float(data['read2_before_filtering']['q20_bases'])/(data['read2_before_filtering']['total_bases'])*100
q20_clean1 = float(data['read1_after_filtering']['q20_bases'])/(data['read1_after_filtering']['total_bases'])*100
q20_clean2 = float(data['read2_after_filtering']['q20_bases'])/(data['read2_after_filtering']['total_bases'])*100
q30_raw1 = float(data['read1_before_filtering']['q30_bases'])/(data['read1_before_filtering']['total_bases'])*100
q30_raw2 = float(data['read2_before_filtering']['q30_bases'])/(data['read2_before_filtering']['total_bases'])*100
q30_clean1 = float(data['read1_after_filtering']['q30_bases'])/(data['read1_after_filtering']['total_bases'])*100
q30_clean2 = float(data['read2_after_filtering']['q30_bases'])/(data['read2_after_filtering']['total_bases'])*100
gc_raw1 = mean(data['read1_before_filtering']['content_curves']['GC'])*100
gc_raw2 = mean(data['read2_before_filtering']['content_curves']['GC'])*100
gc_clean1 = mean(data['read1_after_filtering']['content_curves']['GC'])*100
gc_clean2 = mean(data['read2_after_filtering']['content_curves']['GC'])*100
qs_raw1 = mean(data['read1_before_filtering']['quality_curves']['mean'])
qs_raw2 = mean(data['read2_before_filtering']['quality_curves']['mean'])
qs_clean1 = mean(data['read1_after_filtering']['quality_curves']['mean'])
qs_clean2 = mean(data['read2_after_filtering']['quality_curves']['mean'])

total_reads = data['summary']['before_filtering']['total_reads']/2
clean_reads = data['summary']['after_filtering']['total_reads']/2
low_quality_reads = data['filtering_result']['low_quality_reads']/2
too_many_N_reads = data['filtering_result']['too_many_N_reads']/2
adapter_reads = data['filtering_result']['too_short_reads']/2
datasize = float(data['summary']['before_filtering']['total_bases'])/1e9
effective = float(clean_reads)/total_reads*100
err_rate = qs2error(mean([qs_clean1,qs_clean2]))
q20 = float(data['summary']['after_filtering']['q20_bases'])/data['summary']['after_filtering']['total_bases']*100
q30 = float(data['summary']['after_filtering']['q30_bases'])/data['summary']['after_filtering']['total_bases']*100
gc_content = data['summary']['after_filtering']['gc_content']*100
low_qual_reads = float((low_quality_reads)/total_reads*100)
adapter_reads = float((adapter_reads)/total_reads*100)
summary = open(os.path.join(outdir,"%s.summary"%prefix),'w') 
#summary = open(os.path.join(outdir,"%s"%prefix),'w') # took the ".summary" out here
#summary.write("Sample\tRaw reads\tRaw data(G)\tEffective(%s)\tError(%)\tQ20(%)\tQ30(%)\tLowQualReads\tAdapterReads\tGC(%)\n")
summary.write("%s\t%d\t%0.1f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\n"%(prefix,total_reads,datasize,effective,err_rate,q20,q30,low_qual_reads,adapter_reads,gc_content))
summary.close()
