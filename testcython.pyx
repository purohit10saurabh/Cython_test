from libcpp.string cimport string
from libcpp.pair cimport pair
from libcpp cimport bool
from libcpp.vector cimport vector
from scipy.sparse import *
import numpy as np
import time
import random
from c_mat cimport *
from p_mat cimport *
from testcython cimport *

def my_c_func( csc_score_mat, csc_Y_token, p_alpha):
	cdef P_SMatF p_score_mat = SMatF_csc_to_cpp( csc_score_mat )
	cdef P_SMatF p_Y_token = SMatF_csc_to_cpp( csc_Y_token )
	cdef SMatF* score_mat = p_score_mat.smat;
	cdef SMatF* Y_token = p_Y_token.smat;
	cdef float alpha = p_alpha;
	cdef SMatF* lbl_score_smat = my_func(score_mat,Y_token,alpha)
	lbl_score_p_smat = P_SMatF.wrap(lbl_score_smat)
	lbl_score_csc_smat = SMatF_cpp_to_csc(lbl_score_p_smat)
	return lbl_score_csc_smat

def my_GM_func( csc_score_mat, csc_Y_token, p_alpha):
	cdef P_SMatF p_score_mat = SMatF_csc_to_cpp( csc_score_mat )
	cdef P_SMatF p_Y_token = SMatF_csc_to_cpp( csc_Y_token )
	cdef SMatF* score_mat = p_score_mat.smat;
	cdef SMatF* Y_token = p_Y_token.smat;
	cdef float alpha = p_alpha;
	cdef SMatF* lbl_score_smat = my_GM(score_mat,Y_token,alpha)
	lbl_score_p_smat = P_SMatF.wrap(lbl_score_smat)
	lbl_score_csc_smat = SMatF_cpp_to_csc(lbl_score_p_smat)
	return lbl_score_csc_smat
