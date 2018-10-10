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
import glob
from progress.bar import Bar

#rebin function
def rebin8x8(mydata):
	#8x8
	newdata8 = np.zeros(np.shape(mydata),dtype=np.float32)
	for i in range(0,len(mydata),8):
		for j in range(0,len(newdata8),8):
			(a) = (mydata[i,j]+mydata[i,j+1]+mydata[i,j+2]+mydata[i,j+3]+mydata[i,j+4]+mydata[i,j+5]+mydata[i,j+6]+mydata[i,j+7]
				+mydata[i+1,j]+mydata[i+1,j+1]+mydata[i+1,j+2]+mydata[i+1,j+3]+mydata[i+1,j+4]+mydata[i+1,j+5]+mydata[i+1,j+6]+mydata[i+1,j+7]
				+mydata[i+2,j]+mydata[i+2,j+1]+mydata[i+2,j+2]+mydata[i+2,j+3]+mydata[i+2,j+4]+mydata[i+2,j+5]+mydata[i+2,j+6]+mydata[i+2,j+7]
				+mydata[i+3,j]+mydata[i+3,j+1]+mydata[i+3,j+2]+mydata[i+3,j+3]+mydata[i+3,j+4]+mydata[i+3,j+5]+mydata[i+3,j+6]+mydata[i+3,j+7]
				+mydata[i+4,j]+mydata[i+4,j+1]+mydata[i+4,j+2]+mydata[i+4,j+3]+mydata[i+4,j+4]+mydata[i+4,j+5]+mydata[i+4,j+6]+mydata[i+4,j+7]
				+mydata[i+5,j]+mydata[i+5,j+1]+mydata[i+5,j+2]+mydata[i+5,j+3]+mydata[i+5,j+4]+mydata[i+5,j+5]+mydata[i+5,j+6]+mydata[i+5,j+7]
				+mydata[i+6,j]+mydata[i+6,j+1]+mydata[i+6,j+2]+mydata[i+6,j+3]+mydata[i+6,j+4]+mydata[i+6,j+5]+mydata[i+6,j+6]+mydata[i+6,j+7]
				+mydata[i+7,j]+mydata[i+7,j+1]+mydata[i+7,j+2]+mydata[i+7,j+3]+mydata[i+7,j+4]+mydata[i+7,j+5]+mydata[i+7,j+6]+mydata[i+7,j+7])
			b = a
			newdata8[i,j] = b
			newdata8[i,j+1] = b
			newdata8[i,j+2] = b
			newdata8[i,j+3] = b
			newdata8[i,j+4] = b
			newdata8[i,j+5] = b
			newdata8[i,j+6] = b
			newdata8[i,j+7] = b
			newdata8[i+1,j] = b
			newdata8[i+1,j+1] = b
			newdata8[i+1,j+2] = b
			newdata8[i+1,j+3] = b
			newdata8[i+1,j+4] = b
			newdata8[i+1,j+5] = b
			newdata8[i+1,j+6] = b
			newdata8[i+1,j+7] = b
			newdata8[i+2,j] = b
			newdata8[i+2,j+1] = b
			newdata8[i+2,j+2] = b
			newdata8[i+2,j+3] = b
			newdata8[i+2,j+4] = b
			newdata8[i+2,j+5] = b
			newdata8[i+2,j+6] = b
			newdata8[i+2,j+7] = b
			newdata8[i+3,j] = b
			newdata8[i+3,j+1] = b
			newdata8[i+3,j+2] = b
			newdata8[i+3,j+3] = b
			newdata8[i+3,j+4] = b
			newdata8[i+3,j+5] = b
			newdata8[i+3,j+6] = b
			newdata8[i+3,j+7] = b
			newdata8[i+4,j] = b
			newdata8[i+4,j+1] = b
			newdata8[i+4,j+2] = b
			newdata8[i+4,j+3] = b
			newdata8[i+4,j+4] = b
			newdata8[i+4,j+5] = b
			newdata8[i+4,j+6] = b
			newdata8[i+4,j+7] = b
			newdata8[i+5,j] = b
			newdata8[i+5,j+1] = b
			newdata8[i+5,j+2] = b
			newdata8[i+5,j+3] = b
			newdata8[i+5,j+4] = b
			newdata8[i+5,j+5] = b
			newdata8[i+5,j+6] = b
			newdata8[i+5,j+7] = b
			newdata8[i+6,j] = b
			newdata8[i+6,j+1] = b
			newdata8[i+6,j+2] = b
			newdata8[i+6,j+3] = b
			newdata8[i+6,j+4] = b
			newdata8[i+6,j+5] = b
			newdata8[i+6,j+6] = b
			newdata8[i+6,j+7] = b
			newdata8[i+7,j] = b
			newdata8[i+7,j+1] = b
			newdata8[i+7,j+2] = b
			newdata8[i+7,j+3] = b
			newdata8[i+7,j+4] = b
			newdata8[i+7,j+5] = b
			newdata8[i+7,j+6] = b
			newdata8[i+7,j+7] = b
	return newdata8

rt_path = './Stokes_IQUV_rebin/'
st_path = './Stokes_IQUV_averages/'
i_path = './Stokes_IQUV_medians/'

files = glob.glob(st_path+'hmi.S_90s.20170910_*_TAI.3I_*')
list_ii = [i.split('I_')[0].split('/')[-1] for i in files]

bar = Bar('Processing', max=len(list_ii))

for idx,inter in enumerate(list_ii):
	bar.next()
	"Make <I>"
	fname0 = st_path+'%sI_recentered.medians_4pics.meanvalue.fits'%inter
	i_rebin = rebin8x8(sunpy.map.Map(fname0).data)
	i_head = sunpy.map.Map(fname0).meta
	i_map = sunpy.map.Map(i_rebin,i_head)
	rot_ang = i_map.meta['crota2']
	i_rot = i_map.rotate(angle = rot_ang * u.deg)
	i_rot.save(rt_path+'%sI_recentered.medians_4pics.meanvalue.rebin8x8.fits'%inter)
	del fname0,i_head,i_map,i_rebin,i_rot,rot_ang
	"Make <Q>"
	fname0 = st_path+'%sQ_recentered.medians_4pics.meanvalue.fits'%inter
	q_rebin = rebin8x8(sunpy.map.Map(fname0).data)
	q_head = sunpy.map.Map(fname0).meta
	q_map = sunpy.map.Map(q_rebin,q_head)
	rot_ang = q_map.meta['crota2']
	q_rot = q_map.rotate(angle = rot_ang * u.deg)
	q_map.save(rt_path+'%sQ_recentered.medians_4pics.meanvalue.rebin8x8.fits'%inter)
	del fname0,q_head,q_map,q_rebin,q_rot,rot_ang
	"Make <V>"
	fname0 = st_path+'%sV_recentered.medians_4pics.meanvalue.fits'%inter
	v_rebin = rebin8x8(sunpy.map.Map(fname0).data)
	v_head = sunpy.map.Map(fname0).meta
	v_map = sunpy.map.Map(v_rebin,v_head)
	rot_ang = v_map.meta['crota2']
	v_rot = v_map.rotate(angle = rot_ang * u.deg)
	v_map.save(rt_path+'%sV_recentered.medians_4pics.meanvalue.rebin8x8.fits'%inter)
	del fname0,v_head,v_map,v_rebin,v_rot,rot_ang
	"Make <U>"
	fname0 = st_path+'%sU_recentered.medians_4pics.meanvalue.fits'%inter
	u_rebin = rebin8x8(sunpy.map.Map(fname0).data)
	u_head = sunpy.map.Map(fname0).meta
	u_map = sunpy.map.Map(u_rebin,u_head)
	rot_ang = u_map.meta['crota2']
	u_rot = u_map.rotate(angle = rot_ang * u.deg)
	u_map.save(rt_path+'%sU_recentered.medians_4pics.meanvalue.rebin8x8.fits'%inter)
	del fname0,u_head,u_map,u_rebin,u_rot,rot_ang
bar.finish()
