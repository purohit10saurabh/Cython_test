from libcpp.string cimport string
from libcpp.pair cimport pair
from libcpp cimport bool
from libcpp.vector cimport vector
from scipy.sparse import csc_matrix
import numpy as np
import time
from p_mat cimport *
from c_mat cimport *

cdef extern from "cpp_header.h":
	SMatF* my_func( SMatF* score_mat, SMatF* Y_token, float alpha )
	SMatF* my_GM( SMatF* score_mat, SMatF* Y_token, float alpha )