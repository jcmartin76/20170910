#!/usr/bin/env python
"""
This code read the QUV/I HMI Fits and make medians of the images each 6 mins
Also 
Author: Juan Camilo Guevara Gomez
Latest version:  July 2018
"""
import numpy as np
import sunpy
import sunpy.map
import os 
from os import listdir
from os.path import isfile, join
#Reading name of files
mypath = 'Stokes_recentered/'
filenames = [f for f in listdir(mypath) if isfile(join(mypath, f))]
filenames = sorted(filenames)
fnames = sorted(filenames,key = lambda x: x.split('.')[-2])
#Function for divide fnames
def chunkIt(seq, num):
	avg = len(seq) / float(num)
	out = []
	last = 0.0
	while last < len(seq):
		out.append(seq[int(last):int(last + avg)])
		last += avg
	return out
#Grouping QUV
N=24
iquv = chunkIt(fnames,N)
#Calculating medians each 6 minutes or 4 images
for item in iquv:
	s_i = chunkIt(item,90)   #Se divide en 90 porque entonces quedan grupos de 4 imagenes correspondientes a 6 minutos
	print('Calculating medians for %s'%s_i[0][0][-18:-5])
	for jitem in s_i:
		for k in range(len(jitem)):
			locals()["im"+str(k)] = sunpy.map.Map(mypath+'/%s'%jitem[k]).data
		images = np.dstack((im0,im1,im2,im3))
		mydata = np.nanmedian(images, axis=2)
		myheader = sunpy.map.Map(mypath+'%s'%jitem[0]).meta
		mymap = sunpy.map.Map(mydata,myheader)
		mymap.save('Stokes_IQUV_medians/%s.medians_4pics.fits'%jitem[0][:-5],filetype='auto')
print('Code has finished')