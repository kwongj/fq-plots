#!/usr/bin/env python
# Script by Jason Kwong
# Estimates and plots insert size histogram and/or read depth coverage of paired-end reads

# Use modern print function from python 3.x
from __future__ import print_function

# Modules
import argparse
from argparse import RawTextHelpFormatter
import sys
import os
import shutil
import re
import StringIO
import subprocess
from subprocess import Popen, PIPE
from subprocess import check_output
import operator
import pandas as pd
from pandas import DataFrame

# Functions
# Log a message to stderr
def msg(*args, **kwargs):
	print(*args, file=sys.stderr, **kwargs)

# Log an error to stderr and quit with non-zero error code
def err(*args, **kwargs):
	msg(*args, **kwargs)
	sys.exit(1);

# Check dependencies
def check_dep(dep):
	devnull = open(os.devnull, 'w')
	checkdep = subprocess.Popen(['which', dep], stdout=devnull, stderr=subprocess.PIPE, close_fds=True)
	output, errcode = checkdep.communicate()
	if checkdep.returncode != 0:
		err('ERROR: Check "{}" is installed correctly and in $PATH.'.format(dep))

# Print banner
def banner():
	print(40*'-')

# Check and format --coords
def check_coords(coords, interval):
	try:
		start = coords.split(':')[0]
		end = coords.split(':')[1]
	except:
		end = ''
	if start == '':
		start = None
	else:
		try:
			start =	int(start)
		except:
			err('ERROR: Please specify coordinates as start:end integers eg. --coords 1000:1200.')
	if end == '':
		end = None
	else:
		try:
			end = int(end)
		except:
			err('ERROR: Please specify coordinates as start:end integers eg. --coords 1000:1200.')
	if start and end:
		if (end - start)/interval < 1:
			err('ERROR: Coordinates span must be greater than interval.')
	return start, end

# Subsample dataframe at specified intervals
# Need to fix index lookup to avoid "UserWarning: Boolean Series key will be reindexed to match DataFrame index."
# integer vs label index lookup
def intervalDF(df, locus, start, end, n):
	if len(df_DEP1[df_DEP1['Locus'] == locus]) < end:
		err('ERROR: End coordinate is greater than length of locus sequence.')
	if start == None and end == None:
		s = n-1
		e = None
	elif start == None:
		s = n-1
		e = df_DEP1[df_DEP1['Locus'] == locus][df_DEP1['Key'] == end].index[0] + 1
	elif end == None:
		s = df_DEP1[df_DEP1['Locus'] == locus][df_DEP1['Key'] == start].index[0]
		e = None
	else:
		s = df_DEP1[df_DEP1['Locus'] == locus][df_DEP1['Key'] == start].index[0]
		e = df_DEP1[df_DEP1['Locus'] == locus][df_DEP1['Key'] == end].index[0] + 1
	return df.iloc[s:e:n, :]

# Run samtools depth
def samtools_depth(bam):
	print('Calculating read depths and plotting histograms ...')
	samtoolsDEPTH = check_output(['samtools', 'depth', '-aa', bam])
	readDEPTHS = StringIO.StringIO()
	readDEPTHS.write(samtoolsDEPTH)
	readDEPTHS.seek(0)
	# Setup Pandas dataframe
	df_DEP1 = pd.read_csv(readDEPTHS, delimiter='\t', header=None, names=['Locus', 'Key', 'Count'])
	return df_DEP1

# Trim outliers and include top centile specified
def percentile(cent, df):
	total_keys = df['Count'].sum()
	cent_keys = (float(cent)/100) * total_keys
	n = 1
	df_sorted = df.sort(['Count'], ascending=False)
	df_cent = df_sorted.head(n)
	while df_cent['Count'].sum() < cent_keys:
		n = n+1
		df_cent = df_sorted.head(n)
	return df_cent.sort(['Key'], ascending=True)

# Run samtools stats
def samtools_stats(bam):
	print('Calculating insert size stats and plotting histogram ...')
	samtoolsSTATS = check_output(['samtools', 'stats', bam])
	readSTATS = StringIO.StringIO()
	readSTATS.write(samtoolsSTATS)
	readSTATS.seek(0)
	dict = {}
	for line in readSTATS:
		line = line.rstrip()
		if re.match('^IS', line, flags=re.MULTILINE):
			row = line.split('\t')
			dict[int(row[1])] = int(row[2])
	# Setup Pandas dataframe
	df_IS1 = pd.DataFrame([[key,value] for key,value in dict.iteritems()],columns=['Key', 'Count'])
	if args.centile:
		df_IS1 = percentile(args.centile[0], df_IS1)
	# Convert frequency dataframe into list
	readLIST = []
	for index,row in df_IS1.iterrows():
		readLIST.append([row['Key']]*row['Count'])
	readLIST = [n for e in readLIST for n in e]
	df_IS2 = pd.DataFrame([[n] for n in readLIST],columns=['Key'])
	return df_IS1, df_IS2

# Calculate stats
def stats(df, col):
	TOTAL = len(df.index)
	try:
		MODE = int(df[col].mode())
	except:
		MODE = 0
	MEAN = float(df[col].mean())
	MEDIAN = int(df[col].median())
	Q25 = int(df[col].quantile(.25))
	Q50 = int(df[col].quantile(.5))
	Q75 = int(df[col].quantile(.75))
	return TOTAL, MODE, MEAN, MEDIAN, Q25, Q50, Q75

# Print insert size stats
def insert_size_stats(results):
	banner()
	print('    Total read pairs:    {}'.format(results[0]))
	print('    Insert size MODE:    {}'.format(results[1]))
	print('                MEAN:    {}'.format(results[2]))
	print('              MEDIAN:    {}'.format(results[3]))
	print('         Q25 Q50 Q75:    {} {} {}\n'.format(results[4], results[5], results[6]))

# Print read depth stats
def read_depth_stats(locus, results):
	banner()
	print('    Ref genome locus:    {}'.format(locus))
	print('     Locus size (bp):    {}'.format(results[0]))
	print('     Read depth MEAN:    {}'.format(results[2]))
	print('              MEDIAN:    {}'.format(results[3]))
	print('         Q25 Q50 Q75:    {} {} {}\n'.format(results[4], results[5], results[6]))

# Draw histogram
def draw_hist(df):
	# Set histogram bar character
	bar = unichr(9608)
	# Retrieve number of columns in terminal
	columns = int(os.popen('stty size', 'r').read().split()[1])
	# Set width of histogram accordingly
	width = columns - 20
	# Retrieve maximum count value
	maxsize = df['Count'].max()
	# Draw histogram scaled to terminal size
	for index, row in df.iterrows():
		key = row['Key']
		count = row['Count']
		scale = float(count)/maxsize
		size = int(scale*width)
		print(key, size*bar, count)
	print('\n')

# Usage
parser = argparse.ArgumentParser(
	formatter_class=RawTextHelpFormatter,
	description='Estimates and plots insert sizes and read depth coverage of paired-end reads\n',
	usage='\n  %(prog)s BAMFILE')
parser.add_argument('input', metavar='BAMFILE', nargs=1, help='BAM file eg. snps.bam from Snippy')
parser.add_argument('--plot', metavar='insert|depth', help='Plot insert sizes ("--plot insert") or read depth ("--plot depth")')
parser.add_argument('--centile', metavar='%', nargs=1, type=int, help='Percentile filter for inserts eg. 95 = 95%% most frequent insert sizes')
parser.add_argument('--locus', metavar='LOCUS', help='Locus to display depth plots')
parser.add_argument('--coords', metavar='START:END', help='Locus coordinates to display depth plots')
parser.add_argument('--interval', metavar='LEN', nargs=1, type=int, help='Interval (in bp) to draw depth plots (default=1000)')
parser.add_argument('--version', action='version', version=
	'=====================================\n'
	'%(prog)s v0.1\n'
	'Updated 30-May-2016\n'
	'Dependencies: Python 2.7.x Samtools\n'
	'=====================================\n')
args = parser.parse_args()

# Calculate stats on dataframe and print histogram
bam = str(args.input[0])
if str.lower(args.plot) == 'depth':
	if args.interval:
		interval = args.interval[0]
	else:
		interval = 1000
	start = None
	end = None
	if args.coords:
		if not args.locus:
			err('ERROR: Please specify locus with --locus.')
		coords = str(args.coords)
		start = check_coords(coords, interval)[0]
		end = check_coords(coords, interval)[1]
	df_DEP1 = samtools_depth(bam)
	results = {}
	loci = sorted(set(df_DEP1['Locus']))
	if args.locus:
		loci = [args.locus]
	df_DEP2 = intervalDF(df_DEP1, loci[0], start, end, interval)
	for locus in loci:
		banner()
		print('LOCUS:' + '\t' + locus)
		banner()
		df_DEP3 = df_DEP2.loc[df_DEP2['Locus'] == locus]
		draw_hist(df_DEP3)
		df_DEP4 = df_DEP1.loc[df_DEP1['Locus'] == locus]
		results[locus] = stats(df_DEP4, 'Count')
	results['All'] = stats(df_DEP1, 'Count')
	for k in sorted(results):
		read_depth_stats(k, results[k])
elif str.lower(args.plot) == 'insert':
	df_IS = samtools_stats(bam)
	df_IS1 = df_IS[0]
	df_IS2 = df_IS[1]
	results_IS = stats(df_IS2, 'Key')
	draw_hist(df_IS1)
	insert_size_stats(results_IS)
else:
	print('Calculating read depths ...')
	df_DEP = samtools_depth(bam)
	df_DEP1 = df_DEP[0]
	df_DEP2 = df_DEP[1]
	results_DEP = stats(df_DEP1, 'Count')
	print('Calculating insert size stats ...')
	df_IS = samtools_stats(bam)
	df_IS1 = df_IS[0]
	df_IS2 = df_IS[1]
	results_IS = stats(df_IS2, 'Key')
	read_depth_stats(results_DEP)
	insert_size_stats(results_IS)
