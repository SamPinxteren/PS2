from collections import defaultdict, Counter
from math import log, ceil, factorial as fact
from random import choices, random, randint, sample, shuffle
import tempfile
import os


def num_to_string(n, m=100):
	s = ""
	while n > 0:
		s = chr(65 + (n % 26)) + s
		n = n // 26
	return s.rjust(int(ceil(log(m) / log(26))), 'A')

g_dyn = dict()
def g(bars, elements):
	if (bars,elements) not in g_dyn:
		g_dyn[(bars, elements)] = fact(bars + elements) / (fact(bars) * fact(elements))
	return g_dyn[(bars, elements)]

def C_i(X, absolute=False):
	V = [0] * sum(X)
	l = 0 
	norm_term = 1 
	count = 0 

	if absolute:
		V[0] = 1 
	else:
		V[0] = 1.0 

	for i in range(1, len(X)):
		l += X[i-1]
		if not absolute: norm_term = g(l, X[i])
		prev = zip(range(0, l), V[:l])
		V[:l] = [0] * l 
		for j, f in prev:
			for p in range(j + 1, l + 1): 
				for t in range(X[i]):
					V[p + t] += f * g(j, t) * g(l-p, X[i]-t-1) / norm_term
					count += 1
	return sum(V)

def generate_random_dataset(L, S, P, pref):
	"""
	Generates an artificial dataset using the given parameters:
	 - L: The size of the alphabet.
	 - S: The number of sequences
	 - P: The length of the sequences
	Returns a list of sequences containing the pattern with an
	associated p-value of at most the given value.
	"""
	prefix = "{}{}".format(L, pref)
	pattern = [prefix + num_to_string(i, L) for i in range(L)]
	dataset = []
	total_p = 1
	for i in range(S):
		sequence = choices(pattern, k=P)
		dataset.append(sequence)

	return (dataset, " ".join(pattern))


def generate_artificial_dataset(L, P, occurs_prob=0.05, pref=""):
	"""
	Generates an artificial dataset using the given parameters:
	 - L: The length of the pattern.
	 - P: The minimum associated p-value.
	 - p: The probability for the sequence to occur.
	Returns a list of sequences containing the pattern with an
	associated p-value of at most the given value.
	"""
	prefix = "{}{}".format(L, pref)
	pattern = [prefix + num_to_string(i, L) for i in range(L)]
	dataset = []
	total_p = 1
	occur_sequences = 0
	no_occur_sequences = 0
	while total_p > P:
		if random() <= occurs_prob:
			sequence = choices(pattern, k=randint(0, L)) + pattern + choices(pattern, k=randint(0, L))
			p = C_i(list(Counter(sequence).values()))
			total_p *= p
			occur_sequences += 1
		else:
			symbols = sample(pattern, randint(1, L-1))
			sequence = choices(symbols, k=randint(1, 3*L))
			no_occur_sequences += 1
		dataset.append(sequence)

	return (dataset, " ".join(pattern))

def get_support_value(filename, n):
	"""
	For a given filename, get the number of percentage
	of lines which make up n lines.
	"""
	with open(filename, 'r') as f:
		line_count = sum([1 for line in f])
		f.close()
		return float(n) / line_count

numeric_dict = defaultdict(lambda: len(numeric_dict))
def file_to_numeric(filename):
	"""
	Convert a space seperated file to the numeric input
	expected by certain algorithms as implemented in
	SPMF.
	Creates a file that should be removed after use.
	"""
	global numeric_dict
	fd, filename_output = tempfile.mkstemp()
	try:
		with open(filename, 'r') as in_file:
			with os.fdopen(fd, 'w') as f:
				for line in in_file:
					for word in line.rstrip().split(" "):
						f.write("{} -1 ".format(numeric_dict[word]))
					f.write("-2\n")
	except:
		os.close(fd)
		os.remove(filename_output)
		raise
	return filename_output

def file_from_numeric(filename):
	"""
	Convert the numeric output of an SPMF algorithm
	into a string file.
	"""
	global numeric_dict
	reverse_dict = {v: k for k, v in numeric_dict.items()}
	try:
		with open(filename, 'r') as in_file:
			fd, filename_output = tempfile.mkstemp()
			with os.fdopen(fd, 'w') as f:
				for line in in_file:
					words = []
					data_part = False
					for word in line.rstrip().split(" "):
						if data_part or word.startswith("#"):
							data_part = True
							words.append(word)
							continue
						try:
							if int(word) in reverse_dict:
								words.append(reverse_dict[int(word)])
							else:
								words.append(word)
						except:
							words.append(word)
					f.write(" ".join(words) + "\n")
			f.close()
		in_file.close()
	except:
		os.remove(filename_output)
		os.close(filename_output)
		raise
	return filename_output

def get_p_scores(pattern_file, data_file, prob_command):
	"""
	Read an SPMF output file and calculate p values
	for each resulting pattern.
	 - Any pattern of length 1 will be skipped.
	 - Any duplicate words will be skipped
	"""
	# Convert output file to clean pattern file
	fd, pattern_file_conv = tempfile.mkstemp()
	try:
		with open(pattern_file, 'r') as in_file:
			with os.fdopen(fd, 'w') as f:
				for line in in_file:
					words = []
					seen = set()
					for word in line.rstrip().split(" "):
						if word.startswith("#"): break
						try:
							int(word)
						except:
							if word not in seen:
								words.append(word)
								seen.add(word)
					if len(words) <= 1: continue
					f.write(" ".join(words) + "\n")
			f.close()
		in_file.close()
	except:
		os.close(fd)
		os.remove(pattern_file_conv)
		raise

	# Calculate p-values
	fd, output_file = tempfile.mkstemp()
	os.close(fd)
	command = "{command} {data} {pattern} > {output}".format(
		command=prob_command,
		data=data_file,
		pattern=pattern_file_conv,
		output=output_file
	)
	print("Running: {}".format(command))
	os.system(command)
	os.remove(pattern_file_conv)

	return output_file

def sort_and_remove_singles(input_file):
	"""
	Sorts the contents of a file by the number after the
	first colon (:). Removes lines which only contain
	one -1. On the output of an SPMF algorithm this means
	removing any pattern of length 1.
	Creates a new file containing the result and returns
	the filename.
	"""
	fd, output_file = tempfile.mkstemp()
	os.close(fd)
	command = "sed -n '/.*-1.*-1.*/p' {in_file} | sort -n -t ':' -k 2 -r > {out_file}".format(
		in_file=input_file,
		out_file=output_file
	)
	os.system(command)
	return output_file

def sort_file(input_file, reverse=True):
	if reverse:
		r = "r"
	else:
		r = ""
	os.system("sort -{r}un -o {in_file} {in_file}".format(r=r,in_file=input_file))

def find_on_line(filename, patterns):
	"""
	Find the lines on which each of the
	given patterns are found.
	"""
	d = defaultdict(lambda: "/")
	with open(filename, 'r') as f:
		i = 1
		for line in f:
			pattern = " ".join(line.rstrip().split(" ")[1:])
			if pattern in patterns and d[pattern] == "/":
				d[pattern] = i
			i += 1
	return d

def generate_datafile_bursts(bursty_p, normal_p):
	fd, out_file = tempfile.mkstemp()
	patterns = []
	line_count = 0
	P = (1/fact(7)) ** 3
	with os.fdopen(fd, 'w') as f:
		datasets = []
		for L in range(3,7):
			for bursty in [True, False]:
				if bursty:
					p = bursty_p
					prefix = "B"
				else:
					p = normal_p
					prefix = "N"

				dataset, pattern = generate_artificial_dataset(L, P, p, prefix)
				patterns.append(pattern)
				for sequence in dataset:
					f.write(" ".join(sequence) + "\n")
					line_count += 1
	return {
		'filename': out_file,
		'support': "{:0.5f}%".format(1.5 * 100/line_count),
		'patterns': patterns
	}

def generate_datafile_randombursts(blocks=0):
	fd, out_file = tempfile.mkstemp()
	patterns = []
	line_count = 0
	P = (1/fact(6)) ** 3
	with os.fdopen(fd, 'w') as f:
		datasets = []
		# Generate data
		for L in range(3,7):
			dataset, pattern = generate_artificial_dataset(L, P, 0.15, "N")
			patterns.append(pattern)
			for sequence in dataset:
				f.write(" ".join(sequence) + "\n")
				line_count += 1

		# Generate random bursts
		for i in range(blocks):
			dataset, pattern = generate_random_dataset(4, 10, 20, "B|{}|".format(i))
			#patterns.append(pattern)
			for sequence in dataset:
				f.write(" ".join(sequence) + "\n")
				line_count += 1
	return {
		'filename': out_file,
		'support': "{:0.5f}%".format(1.5 * 100.0/line_count),
		'patterns': patterns
	}

def generate_datafile():
	fd, out_file = tempfile.mkstemp()
	patterns = []
	line_count = 0
	P = (1/fact(6)) ** 3
	with os.fdopen(fd, 'w') as f:
		datasets = []
		for L in range(2,7):
			dataset, pattern = generate_artificial_dataset(L, P)
			patterns.append(pattern)
			for sequence in dataset:
				f.write(" ".join(sequence) + "\n")
				line_count += 1
	return {
		'filename': out_file,
		'support': "{:0.5f}%".format(1.5 * 100/line_count),
		'patterns': patterns
	}

def randomize_sequence_order(filename):
	fd, out_file = tempfile.mkstemp()
	with os.fdopen(fd, 'w') as f:
		with open(filename, 'r') as in_f:
			for line in in_f:
				tokens = line.rstrip().split(" ")
				shuffle(tokens)
				f.write(" ".join(tokens) + "\n")
	return out_file
