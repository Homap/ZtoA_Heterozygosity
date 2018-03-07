#!/usr/bin/python
import sys
#test_var = 10 Try in ipython modulename.test_var
#sequence = open(sys.argv[1], "r")
#print "Enter window size", # Comma after the print statement in python 2 removes '\n'
#winSize = raw_input()
#print "Enter step size",
#step = raw_input()

# Function for sliding window
def slidingWindow(sequence, winSize, step):
	""" Returns a generator that will iterate through
		the defined chunks of input sequence. Input
		sequence must be iterable."""
		
	# Verify the inputs
	try: it = iter(sequence)
	except TypeError:
		raise Exception ("**ERROR** sequence must be iterable.")
	if not ((type(winSize) == type(0)) and (type(step) == type(0))):
		raise Exception("**ERROR** type(winSize) and type(step) must be int.")
	if step > winSize:
		raise Exception("**ERROR** step must not be larger than winSize.")
	if winSize > len(sequence):
#		print sequence, "NOPE"
#		print "NOPE"
		pass
#		raise Exception("**ERROR** winSize must not be larger than sequence length.")
		
	# Pre-compute number of chunks to emit
	numOfChunks = ((int(len(sequence)-winSize)/step))+1
	numOfChunks = int(numOfChunks)
	
	# Do the work
	for i in range(0, numOfChunks*step, step):
		yield sequence[i:i+winSize]
