"""
LMK_UTILS 
@author: Liam M. Kilcommons
My convenince functions and debugging tools.

"""
 
import pdb, traceback, sys
import numpy

def nNans(arr,column=0):
	"""
	Counts up the number of occurances of NaN in a row of a 2-d array
	"""
	if len(arr.shape) != 2:
		raise RuntimeError('nanRows works only on 2-d arrays, this array had shape %s' % (str(arr.shape)))

	nans = arr[:,column] == numpy.nan
	return len(numpy.nonzero(nans)[0])

def asColumn(arr):
	"""
	Reshapes a 1 x n numpy array to be n x 1. Throws an error if 
	input is not 1 x n or n x 1.
	"""

	if not isinstance(arr,numpy.ndarray):
		raise TypeError('asColumn expects as numpy array, got type %s') % (type(arr))

	#Reshape things to be column vector if nessecary
	if len(arr.shape)==1: #its currently 1-d, make it 2-d and a column array
		arr = numpy.resize(arr,(arr.size,1))
	elif arr.shape[1] > 1 and arr.shape[0] == 1:
			arr = numpy.resize(arr,(max(arr.shape),1))
	elif arr.shape[0] > 1 and arr.shape[1] == 1:
		pass #do nothing
	else:
		raise RuntimeError('asColumn was passed a non 1-d numpy array. Actual dimensions %d x %d' % (arr.shape[0],arr.shape[1]))
	return arr		

def nan(shape):
	"""
	Shorthand for new nan-filled numpy array
	"""

	if len(shape) != 2 or not isinstance(shape,tuple):
		raise ValueError('Bad shape input, must be tuple of length 2')

	arr = numpy.empty(shape)
	arr.fill(numpy.nan)
	return arr

def deduplicate(arr,column=0,terse=False):
	"""
	Show and remove any duplicate values in a numpy array
	"""
	if arr.shape[1] > arr.shape[0]:
		raise RuntimeError('This function only works on arrays which have one timestamp per row (spreadsheet style)')

	unqvals,unqinv = numpy.unique(arr[:,column],return_inverse=True)
	unqvals,unqind = numpy.unique(arr[:,column],return_index=True)

	#Count up the number of occurances of each unique value
	n_occurances = numpy.zeros_like(unqvals)
	for inv in unqinv:
		n_occurances[inv] += 1

	#Find the indices of any value with duplicates
	multiple_offenders = numpy.nonzero(n_occurances>1)

	mo_vals = unqvals[multiple_offenders]
	mo_counts = n_occurances[multiple_offenders]

	if not terse:
		for k in range(len(mo_vals)):
			print("Duplicate value %s occured %d times in column %d" % (mo_vals[k],mo_counts[k],column)) 
			print(str(arr[arr[:,column]==mo_vals[k],:]))
	else:
		print("Of %d input values, %d were unique, %d had more than one occurance." % (len(arr[:,column]),len(unqvals),len(mo_vals)))

	return arr[unqind,:]

def arraycheck(inobj,returnList=False):
	"""
	I made this to help input debugging.
	Check to see if an input is a numpy array,
	-if its an array and returnList is false do nothing
	-if it's a list, make it a column vector
	-if it's a single value, make it a list and then a column vector
	-if it's neither throw a TypeError
	-if returnList make the output a list
	-If returning an array, make sure it's always at least 2-d
	"""
	if isinstance(inobj,numpy.ndarray):
		pass
	else:
		#Okay it's not a numpy array
		#Does it act like a list and have a length? 
		try: 
			nelements = len(inobj)
		except TypeError: #if not, it's probably a single number,
						  #so make it a list (Duck Typing)
			nelements = 1 
			inobj = [inobj]

		inobj = numpy.array(inobj,ndmin=2) #Make sure its a column array

	if returnList:
		return inobj.tolist()
	else:
		return inobj


	

