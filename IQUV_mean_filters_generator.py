#!/usr/bin/env python

"""
Generate <I>,<Q>,<U>,<V> from the medians already calcaluated.

Author: Juan Camilo Guevara Gomez
Latest version: 17 August 2018
"""
import multiprocessing as mp
import csv
import astropy.units as u
import matplotlib
import matplotlib.pyplot as plt
import urllib.request
import numpy as np
from matplotlib.patches import Circle
from PIL import Image
import cv2
import rawpy
from IPython.display import clear_output
from IPython import display
import time,os,sys
import datetime
from os import listdir
from os.path import isfile, join
import numpy.ma as ma
import sunpy 
import sunpy.map
from progress.bar import Bar
import glob

st_path = './Stokes_IQUV_averages/'
i_path = './Stokes_IQUV_medians/'

files = glob.glob(i_path+'hmi.S_90s.20170910_*_TAI.3.I0*')
list_ii = [i.split('.I0')[0].split('/')[-1] for i in files]


bar = Bar('Processing', max=len(list_ii))

for idx,inter in enumerate(list_ii):
	bar.next()
	"Make <I>"
	fname0 = i_path+'%s.I0_recentered.medians_4pics.fits'%inter
	fname1 = i_path+'%s.I1_recentered.medians_4pics.fits'%inter
	fname2 = i_path+'%s.I2_recentered.medians_4pics.fits'%inter
	fname3 = i_path+'%s.I3_recentered.medians_4pics.fits'%inter
	fname4 = i_path+'%s.I4_recentered.medians_4pics.fits'%inter
	fname5 = i_path+'%s.I5_recentered.medians_4pics.fits'%inter
	i_chunk = np.array([sunpy.map.Map(fname0).data,sunpy.map.Map(fname1).data,sunpy.map.Map(fname2).data,sunpy.map.Map(fname3).data,sunpy.map.Map(fname4).data,sunpy.map.Map(fname5).data])
	i_avg = np.nanmean(i_chunk, axis=0)
	i_head = sunpy.map.Map(fname0).meta
	i_map = sunpy.map.Map(i_avg,i_head)
	#rot_ang = i_map.meta['crota2']
	#i_rot = i_map.rotate(angle = rot_ang * u.deg)
	i_map.save(st_path+'%sI_recentered.medians_4pics.meanvalue.fits'%inter)
	del fname0,fname1,fname2,fname3,fname4,fname5,i_chunk,i_avg,i_head,i_map
	"Make <Q>"
	fname0 = i_path+'%s.Q0_recentered.medians_4pics.fits'%inter
	fname1 = i_path+'%s.Q1_recentered.medians_4pics.fits'%inter
	fname2 = i_path+'%s.Q2_recentered.medians_4pics.fits'%inter
	fname3 = i_path+'%s.Q3_recentered.medians_4pics.fits'%inter
	fname4 = i_path+'%s.Q4_recentered.medians_4pics.fits'%inter
	fname5 = i_path+'%s.Q5_recentered.medians_4pics.fits'%inter
	q_chunk = np.array([sunpy.map.Map(fname0).data,sunpy.map.Map(fname1).data,sunpy.map.Map(fname2).data,sunpy.map.Map(fname3).data,sunpy.map.Map(fname4).data,sunpy.map.Map(fname5).data])
	q_avg = np.nanmean(q_chunk, axis=0)
	q_head = sunpy.map.Map(fname0).meta
	q_map = sunpy.map.Map(q_avg,q_head)
#	rot_ang = q_map.meta['crota2']
#	q_rot = q_map.rotate(angle = rot_ang * u.deg)
	q_map.save(st_path+'%sQ_recentered.medians_4pics.meanvalue.fits'%inter)
	del fname0,fname1,fname2,fname3,fname4,fname5,q_chunk,q_avg,q_head,q_map
	"Make <V>"
	fname0 = i_path+'%s.V0_recentered.medians_4pics.fits'%inter
	fname1 = i_path+'%s.V1_recentered.medians_4pics.fits'%inter
	fname2 = i_path+'%s.V2_recentered.medians_4pics.fits'%inter
	fname3 = i_path+'%s.V3_recentered.medians_4pics.fits'%inter
	fname4 = i_path+'%s.V4_recentered.medians_4pics.fits'%inter
	fname5 = i_path+'%s.V5_recentered.medians_4pics.fits'%inter
	v_chunk = np.array([sunpy.map.Map(fname0).data,sunpy.map.Map(fname1).data,sunpy.map.Map(fname2).data,sunpy.map.Map(fname3).data,sunpy.map.Map(fname4).data,sunpy.map.Map(fname5).data])
	v_avg = np.nanmean(v_chunk, axis=0)
	v_head = sunpy.map.Map(fname0).meta
	v_map = sunpy.map.Map(v_avg,v_head)
#	rot_ang = v_map.meta['crota2']
#	v_rot = v_map.rotate(angle = rot_ang * u.deg)
	v_map.save(st_path+'%sV_recentered.medians_4pics.meanvalue.fits'%inter)
	del fname0,fname1,fname2,fname3,fname4,fname5,v_chunk,v_avg,v_head,v_map
	"Make <U>"
	fname0 = i_path+'%s.U0_recentered.medians_4pics.fits'%inter
	fname1 = i_path+'%s.U1_recentered.medians_4pics.fits'%inter
	fname2 = i_path+'%s.U2_recentered.medians_4pics.fits'%inter
	fname3 = i_path+'%s.U3_recentered.medians_4pics.fits'%inter
	fname4 = i_path+'%s.U4_recentered.medians_4pics.fits'%inter
	fname5 = i_path+'%s.U5_recentered.medians_4pics.fits'%inter
	u_chunk = np.array([sunpy.map.Map(fname0).data,sunpy.map.Map(fname1).data,sunpy.map.Map(fname2).data,sunpy.map.Map(fname3).data,sunpy.map.Map(fname4).data,sunpy.map.Map(fname5).data])
	u_avg = np.nanmean(u_chunk, axis=0)
	u_head = sunpy.map.Map(fname0).meta
	u_map = sunpy.map.Map(u_avg,u_head)
	#rot_ang = u_map.meta['crota2']
	#u_rot = u_map.rotate(angle = rot_ang * u.deg)
	u_map.save(st_path+'%sU_recentered.medians_4pics.meanvalue.fits'%inter)
	del fname0,fname1,fname2,fname3,fname4,fname5,u_chunk,u_avg,u_head,u_map
bar.finish()
