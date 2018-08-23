"""This module contains many list-based functions to calculate descriptive statistics."""
from math import sqrt 
from collections import Counter

def mode(arr):
   """get the most frequent value"""
   return max(set(arr),key=arr.count)

def average(arr):
  """average of the values, must have more than 0 entries.

  :param arr: list of numbers
  :type arr: number[] a number array
  :return: average
  :rtype: float

  """
  if len(arr) == 0:
    sys.stderr.write("ERROR: no content in array to take average\n")
    sys.exit()
  if len(arr) == 1:  return arr[0]
  return float(sum(arr))/float(len(arr))

def median(arr):
  """median of the values, must have more than 0 entries.

  :param arr: list of numbers
  :type arr: number[] a number array
  :return: median
  :rtype: float

  """
  if len(arr) == 0:
    sys.stderr.write("ERROR: no content in array to take average\n")
    sys.exit()
  if len(arr) == 1: return arr[0]
  quot = int(len(arr)/2)
  rem = len(arr)%2
  if rem != 0:
    return sorted(arr)[quot]
  return float(sum(sorted(arr)[quot-1:quot+1]))/float(2)

def standard_deviation(arr):
  """standard deviation of the values, must have 2 or more entries.

  :param arr: list of numbers
  :type arr: number[] a number array
  :return: standard deviation
  :rtype: float

  """
  return sqrt(variance(arr))

def variance(arr):
  """variance of the values, must have 2 or more entries.

  :param arr: list of numbers
  :type arr: number[] a number array
  :return: variance
  :rtype: float

  """
  avg = average(arr)
  return sum([(float(x)-avg)**2 for x in arr])/float(len(arr)-1)

def N50(arr):
  """N50 often used in assessing denovo assembly.

  :param arr: list of numbers
  :type arr: number[] a number array
  :return: N50
  :rtype: float

  """
  if len(arr) == 0:
    sys.stderr.write("ERROR: no content in array to take N50\n")
    sys.exit()
  tot = sum(arr)
  half = float(tot)/float(2)
  cummulative = 0
  for l in sorted(arr):
    cummulative += l
    if float(cummulative) > half: 
      return l
  sys.stderr.write("ERROR: problem finding M50\n")
  sys.exit()
