from read_iterdb_x import *
from read_EFIT import *
from finite_differences_x import *
from interp import *
import matplotlib.pyplot as plt
import numpy as np
import sys

iterdbFileName = sys.argv[1]
efitFileName = sys.argv[2]
EFITdict = read_EFIT(efitFileName)
ITERDBdict = read_iterdb_x(iterdbFileName)

