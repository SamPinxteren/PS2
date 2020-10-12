#!/usr/bin/env python3

from helper_functions import *
from itertools import product
from collections import defaultdict, Counter
from random import choices, random, randint, sample
from math import log, ceil, factorial as fact
from scipy.stats import norm
import tempfile
import subprocess
import os
import resource
import time

# GoKrimp omitted since it doesn't find any of the planted patterns
algorithms = [
	{
		'command': "java -jar algorithms/spmf.jar run SKOPUS {in_file} {out_file} 1000 true 10 true 1",
		'numeric':  False,
		'sort':     False,
		'sort_p':   False,
		'sigspan':  False,
		'string':  "SKOPUS"
	},
	{
		'command': "java -jar algorithms/spmf.jar run CloFast {in_file} {out_file} {support}",
		'numeric':  True,
		'sort':     True,
		'sort_p':   True,
		'sigspan':  True,
		'string':  "Support"
	}
]

output_filename = "all_on_artificial"
p_command = "algorithms/p -P"
sigspan_command = "algorithms/p -I"

if __name__ == "__main__":
	with open(output_filename, 'w') as output:
		data = defaultdict(list)
		timing_data = defaultdict(list)
		for i in range(25):
			in_file = generate_datafile()
			for algo in algorithms:
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
					support=in_file['support']
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
				timing_data[(algo['string'], 'real')].append(end_wall_time - start_wall_time)
				timing_data[(algo['string'], 'usr')].append(end_time.ru_utime - start_time.ru_utime)
				timing_data[(algo['string'], 'sys')].append(end_time.ru_stime - start_time.ru_stime)
				timing_data[(algo['string'], 'tot')].append(end_time.ru_utime - start_time.ru_utime + end_time.ru_stime - start_time.ru_stime)

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

				# Write output
				output.write("{} {}\n".format(algo['string'], in_file['filename']))
				"""
				with open(out_file_conv, 'r') as f:
					for line in f: output.write(line)
				"""

				# Calculate p-values
				start_time = resource.getrusage(resource.RUSAGE_CHILDREN)
				start_wall_time = time.time()
				scores_file = get_p_scores(out_file_conv, in_file['filename'], p_command)
				end_time = resource.getrusage(resource.RUSAGE_CHILDREN)
				end_wall_time = time.time()
				algo_time = "real: {}, usr: {}, sys: {}, tot: {}".format(
					end_wall_time - start_wall_time,
					end_time.ru_utime - start_time.ru_utime,
					end_time.ru_stime - start_time.ru_stime,
					end_time.ru_utime - start_time.ru_utime + end_time.ru_stime - start_time.ru_stime
				)
				d = find_on_line(scores_file, in_file['patterns'])
				for pattern in in_file['patterns']:
					data[(algo['string'], pattern)].append(d[pattern])
				output.write("\n")
				with open(scores_file, 'r') as f:
					for line in f: output.write(line)

				# Output sorted by p-values
				if algo['sort_p']:
					sort_file(scores_file)
					output.write("\n{} (P-sorted) {}\n".format(algo['string'], in_file['filename']))
					d = find_on_line(scores_file, in_file['patterns'])
					for pattern in in_file['patterns']:
						data[("P-sorted", pattern)].append(d[pattern])
					output.write("\n")
					with open(scores_file, 'r') as f:
						for line in f: output.write(line)
					timing_data[('P-sorted', 'real')].append(end_wall_time - start_wall_time)
					timing_data[('P-sorted', 'usr')].append(end_time.ru_utime - start_time.ru_utime)
					timing_data[('P-sorted', 'sys')].append(end_time.ru_stime - start_time.ru_stime)
					timing_data[('P-sorted', 'tot')].append(end_time.ru_utime - start_time.ru_utime + end_time.ru_stime - start_time.ru_stime)

				output.write("\n")
				output.flush()
				os.remove(scores_file)

				if algo['sigspan']:
					# Calculate sigspan
					start_time = resource.getrusage(resource.RUSAGE_CHILDREN)
					start_wall_time = time.time()
					scores_file = get_p_scores(out_file_conv, in_file['filename'], sigspan_command)
					end_time = resource.getrusage(resource.RUSAGE_CHILDREN)
					end_wall_time = time.time()
					algo_time = "real: {}, usr: {}, sys: {}, tot: {}".format(
						end_wall_time - start_wall_time,
						end_time.ru_utime - start_time.ru_utime,
						end_time.ru_stime - start_time.ru_stime,
						end_time.ru_utime - start_time.ru_utime + end_time.ru_stime - start_time.ru_stime
					)

					# Output sorted by p-values
					sort_file(scores_file)
					output.write("\n{} (SigSpan) {}\n".format(algo['string'], in_file['filename']))
					d = find_on_line(scores_file, in_file['patterns'])
					for pattern in in_file['patterns']:
						data[("SigSpan", pattern)].append(d[pattern])
					output.write("\n")
					with open(scores_file, 'r') as f:
						for line in f: output.write(line)
					timing_data[('SigSpan', 'real')].append(end_wall_time - start_wall_time)
					timing_data[('SigSpan', 'usr')].append(end_time.ru_utime - start_time.ru_utime)
					timing_data[('SigSpan', 'sys')].append(end_time.ru_stime - start_time.ru_stime)
					timing_data[('SigSpan', 'tot')].append(end_time.ru_utime - start_time.ru_utime + end_time.ru_stime - start_time.ru_stime)

					output.write("\n")
					output.flush()
					os.rename(scores_file, "/tmp/sigspan{}".format(i))

				# Clean up files
				if algo['numeric']:
					os.remove(out_file)
					os.remove(in_file_conv)
				os.remove(out_file_conv)
			os.remove(in_file['filename'])
			for algo, pattern in sorted(data.keys()):
				print("{} & {} & ".format(algo, pattern) + " & ".join(map(str, data[(algo, pattern)])))

		print()

		for algo, pattern in sorted(data.keys()):
			l = data[(algo, pattern)]
			output.write("{} & {} & ".format(algo, pattern) + " & ".join(map(str, data[(algo, pattern)])))
			print("{} & {} & {} ({}/{})".format(algo, pattern, sum([int(j) for j in l if str(j).isnumeric()]) / len(l), sum([1 for j in l if j == "/"]), len(l)))

		for algo, stat in sorted(timing_data.keys()):
			l = timing_data[(algo, stat)]
			avg = sum(l) / len(l)
			print(algo, stat, avg)
