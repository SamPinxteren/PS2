#!/usr/bin/python3

from helper_functions import *
from itertools import product
from collections import defaultdict
import tempfile
import subprocess
import os
import time
import resource

output_filename = "all_on_all"
p_command = "algorithms/p -P"
WY_command = "algorithms/p -W 0.05"
sigspan_command = "algorithm/p -i"

input_files = [
	{
		'filename': "datasets/ATaleofTwoCities.txt",
		'support':  "0.1%"
	},
	{
		'filename': "datasets/JMLR.txt",
		'support':  "4%"
	},
	{
		'filename': "datasets/TheNewPhysicsandItsEvolution.txt",
		'support':  "0.1%"
	},
	{
		'filename': "datasets/ABookAboutLawyers.txt",
		'support':  "0.1%"
	},
	{
		'filename': "datasets/AdventuresofHuckleberryFinn.txt",
		'support':  "0.1%"
	},
	{
		'filename': "datasets/TheAdventuresofTomSawyer.txt",
		'support':  "0.1%"
	},
	{
		'filename': "datasets/TheIndustriesofAnimals.txt",
		'support':  "0.1%"
	},
	{
		'filename': "datasets/moby_dick_gutenberg_sentences_clean.txt",
		'support':  "0.1%"
	}
]

algorithms = [
	{
		'command': "java -jar algorithms/promise.jar {in_file} {out_file} 100 10048 {support_nopercent} 32 1",
		'numeric':  True,
		'sort':     True,
		'sort_p':   False,
		'WY':       False,
		'calc_p':   True,
		'sigspan':  False,
		'string':  "ProMiSe"
	},
	{
		'command': "java -jar algorithms/spmf.jar run GoKrimp {in_file} {out_file}",
		'numeric':  True,
		'sort':     False,
		'sort_p':   False,
		'WY':       False,
		'calc_p':   True,
		'sigspan':  False,
		'string':  "GoKrimp"
	},
	{
		'command': "java -jar algorithms/spmf.jar run SKOPUS {in_file} {out_file} 100 true 10 true 1",
		'numeric':  False,
		'sort':     False,
		'sort_p':   False,
		'WY':       False,
		'calc_p':   True,
		'sigspan':  False,
		'string':  "SKOPUS"
	},
	{
		'command': "java -jar algorithms/spmf.jar run CloFast {in_file} {out_file} {support}",
		'numeric':  True,
		'sort':     True,
		'sort_p':   True,
		'WY':       True,
		'calc_p':   True,
		'sigspan':  True,
		'string':  "Support"
	}
]
if __name__ == "__main__":
	with open(output_filename, 'w') as output:
		for in_file, algo in product(input_files, algorithms):
			# Convert to numeric input if needed
			if algo['numeric']:
				in_file_conv = file_to_numeric(in_file['filename'])
			else:
				in_file_conv = in_file['filename']

			# Run algorithm
			fd, out_file = tempfile.mkstemp()
			os.close(fd)
			command = algo['command'].format(
				in_file=in_file_conv,
				out_file=out_file,
				support=in_file['support'],
				support_nopercent=str(float(in_file['support'].replace("%","")) / 100)
			)
			print("Running: {}".format(command))
			start_time = resource.getrusage(resource.RUSAGE_CHILDREN)
			start_wall_time = time.time()
			os.system(command)
			end_time = resource.getrusage(resource.RUSAGE_CHILDREN)
			end_wall_time = time.time()
			algo_time = "real: {}, usr: {}, sys: {}, tot: {}".format(
				end_wall_time - start_wall_time,
				end_time.ru_utime - start_time.ru_utime,
				end_time.ru_stime - start_time.ru_stime,
				end_time.ru_utime - start_time.ru_utime + end_time.ru_stime - start_time.ru_stime
			)

			# Convert back from numeric if needed
			if algo['numeric']:
				out_file_conv = file_from_numeric(out_file)
			else:
				out_file_conv = out_file

			# Sort if requested
			if algo['sort']:
				f_tmp = out_file_conv
				out_file_conv = sort_and_remove_singles(f_tmp)
				os.remove(f_tmp)

			# Calculate p-values
			if algo['calc_p']:
				start_time = resource.getrusage(resource.RUSAGE_CHILDREN)
				start_wall_time = time.time()
				scores_file = get_p_scores(out_file_conv, in_file['filename'], p_command)
				end_time = resource.getrusage(resource.RUSAGE_CHILDREN)
				end_wall_time = time.time()
				p_time = "real: {}, usr: {}, sys: {}, tot: {}".format(
					end_wall_time - start_wall_time,
					end_time.ru_utime - start_time.ru_utime,
					end_time.ru_stime - start_time.ru_stime,
					end_time.ru_utime - start_time.ru_utime + end_time.ru_stime - start_time.ru_stime
				)

			# Get westfall-young
			if algo['WY']:
				start_time = resource.getrusage(resource.RUSAGE_CHILDREN)
				start_wall_time = time.time()
				WY_file = get_p_scores(out_file_conv, in_file['filename'], WY_command)
				end_time = resource.getrusage(resource.RUSAGE_CHILDREN)
				end_wall_time = time.time()
				WY_time = "real: {}, usr: {}, sys: {}, tot: {}".format(
					end_wall_time - start_wall_time,
					end_time.ru_utime - start_time.ru_utime,
					end_time.ru_stime - start_time.ru_stime,
					end_time.ru_utime - start_time.ru_utime + end_time.ru_stime - start_time.ru_stime
				)

			# Write output
			output.write("{} {} (algoT:[{}] pT:[{}])\n".format(algo['string'], in_file['filename'], algo_time, p_time))
			for line in open(out_file_conv): output.write(line)
			if algo['calc_p']:
				for line in open(scores_file): output.write(line)
			if algo['WY']:
				output.write("WYT:{}\n".format(WY_time))
				for line in open(WY_file): output.write(line)
				os.remove(WY_file)


			# Output sorted by p-values
			if algo['sort_p']:
				sort_file(scores_file)
				output.write("\n{} (P-sorted) {}\n".format(algo['string'], in_file['filename']))
				with open(scores_file, 'r') as f:
					for line in f: output.write(line)

			# Compute Sigspan values
			if algo['sigspan']:
				# Calculate sigspan p-values
				start_time = resource.getrusage(resource.RUSAGE_CHILDREN)
				start_wall_time = time.time()
				scores_file = get_p_scores(out_file_conv, in_file['filename'], sigspan_command)
				end_time = resource.getrusage(resource.RUSAGE_CHILDREN)
				end_wall_time = time.time()
				sig_time = "real: {}, usr: {}, sys: {}, tot: {}".format(
					end_wall_time - start_wall_time,
					end_time.ru_utime - start_time.ru_utime,
					end_time.ru_stime - start_time.ru_stime,
					end_time.ru_utime - start_time.ru_utime + end_time.ru_stime - start_time.ru_stime
				)

				sort_file(scores_file, reverse=False)
				output.write("\n{} (Sigspan) {} (T:{})\n".format(algo['string'], in_file['filename'], sig_time))
				with open(scores_file, 'r') as f:
					for line in f: output.write(line)

			output.write("\n")
			output.flush()

			# Clean up files
			if algo['numeric']:
				os.remove(out_file)
				os.remove(in_file_conv)
			os.remove(scores_file)
			os.remove(out_file_conv)
