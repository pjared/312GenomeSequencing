#!/usr/bin/python3

from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
	from PyQt4.QtCore import QLineF, QPointF
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import math
import time
import random

# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1

class GeneSequencing:

	def __init__( self ):
		self.matrix = ''
		pass


	def getCost(self, char1, char2):
		#Match
		if (char1 == char2):
			return -3
		#Substitution
		elif(char1 != char2):
			return 1
		else:
			return 5

	def doAlignment(self, seq1, seq2):
		seq1 = ' ' + seq1
		seq2 = ' ' + seq2
		self.matrix = [[0 for j in range(0, len(seq1))] for i in range(0, len(seq2))]

		#Set up base case
		incrementer = 0
		for i in range(0, len(seq2)):
			self.matrix[i][0] = incrementer
			incrementer += 5
		incrementer = 0
		for i in range(0, len(seq1)):
			self.matrix[0][i] = incrementer
			incrementer += 5
		#TODO: FIX THIS
		for i in range(1,  len(seq2)):
			for j in range(1, len(seq1)):
				temp = seq2[i - 1]
				temp2 = seq1[j - 1]
				totalDist = self.getCost(seq1[j], seq2[i])  # Subtract one since range is starting at 1

				bestDistance = self.matrix[i - 1][j]  # TOP
				topLeft = self.matrix[i - 1][j - 1]
				if (topLeft < bestDistance):
					bestDistance = topLeft
				left = self.matrix[i][j - 1]
				if (left < bestDistance):
					bestDistance = left
				self.matrix[i][j] = totalDist + bestDistance

		#TODO: Make sure that matrix  works, along with getCost()

# This is the method called by the GUI.  _seq1_ and _seq2_ are two sequences to be aligned, _banded_ is a boolean that tells
# you whether you should compute a banded alignment or full alignment, and _align_length_ tells you 
# how many base pairs to use in computing the alignment

	def align(self, seq1, seq2, banded, align_length):
		self.banded = banded
		self.MaxCharactersToAlign = align_length

###################################################################################################
# your code should replace these three statements and populate the three variables: score, alignment1 and alignment2
		score = random.random()*100;
		alignment1 = 'abc-easy  DEBUG:({} chars,align_len={}{})'.format(
			len(seq1), align_length, ',BANDED' if banded else '')
		alignment2 = 'as-123--  DEBUG:({} chars,align_len={}{})'.format(
			len(seq2), align_length, ',BANDED' if banded else '')
###################################################################################################					

		self.doAlignment(seq1, seq2)
		#now i need to get the sequence from the matrix
		self.doAlignment(seq2, seq1)
		return {'align_cost':score, 'seqi_first100':alignment1, 'seqj_first100':alignment2}