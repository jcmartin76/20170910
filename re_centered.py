#!/usr/bin/env python

"""
This code scale and center HMI images
Author: Juan Camilo Guevara Gomez
Version:  July 2018
"""

import numpy as np
import csv
import cv2
import os, os.path
import sys,string
import exifread # to read information from the image
import urllib.request
import astropy.units as u
import matplotlib.pyplot as plt
import pandas as pd
from collections import defaultdict
import matplotlib
from PIL import Image
from PIL import ImageChops
import glob
from os import listdir
from os.path import isfile, join
import os.path
import os
import sunpy
import sunpy.map
import numpy as np
import numpy.ma as ma
import matplotlib.patches as patches

#Reading path with original files
mypath = '/Volumes/Mnemosyne/work_20170910_flare/Stokes'
filenames = ['/Volumes/Mnemosyne/work_20170910_flare/Stokes/'+f for f in listdir(mypath) if isfile(join(mypath, f))]
filenames = sorted(filenames)
#Grouping by Stokes
I0 = []
I1 = []
I2 = []
I3 = []
I4 = []
I5 = []
Q0 = []
Q1 = []
Q2 = []
Q3 = []
Q4 = []
Q5 = []
U0 = []
U1 = []
U2 = []
U3 = []
U4 = []
U5 = []
V0 = []
V1 = []
V2 = []
V3 = []
V4 = []
V5 = []
for item in filenames:
    if item[-7:-5] == 'I0':
        I0.append(item)
    if item[-7:-5] == 'I1':
        I1.append(item)
    if item[-7:-5] == 'I2':
        I2.append(item)
    if item[-7:-5] == 'I3':
        I3.append(item)
    if item[-7:-5] == 'I4':
        I4.append(item)
    if item[-7:-5] == 'I5':
        I5.append(item)
    if item[-7:-5] == 'Q0':
        Q0.append(item)
    if item[-7:-5] == 'Q1':
        Q1.append(item)
    if item[-7:-5] == 'Q2':
        Q2.append(item)
    if item[-7:-5] == 'Q3':
        Q3.append(item)
    if item[-7:-5] == 'Q4':
        Q4.append(item)
    if item[-7:-5] == 'Q5':
        Q5.append(item)
    if item[-7:-5] == 'U0':
        U0.append(item)
    if item[-7:-5] == 'U1':
        U1.append(item)
    if item[-7:-5] == 'U2':
        U2.append(item)
    if item[-7:-5] == 'U3':
        U3.append(item)
    if item[-7:-5] == 'U4':
        U4.append(item)
    if item[-7:-5] == 'U5':
        U5.append(item)
    if item[-7:-5] == 'V0':
        V0.append(item)
    if item[-7:-5] == 'V1':
        V1.append(item)
    if item[-7:-5] == 'V2':
        V2.append(item)
    if item[-7:-5] == 'V3':
        V3.append(item)
    if item[-7:-5] == 'V4':
        V4.append(item)
    if item[-7:-5] == 'V5':
        V5.append(item)
#Concatenating all arrays with files paths
ALL = I0+I1+I2+I3+I4+I5+Q0+Q1+Q2+Q3+Q4+Q5+U0+U1+U2+U3+U4+U5+V0+V1+V2+V3+V4+V5  
total_fits = len(ALL)
#Recentering process
a_fits = 0
for item in ALL:
    a_fits = a_fits + 1
    print('Processing fits %i'%a_fits)
    map_i = sunpy.map.Map(item)
    data_test = map_i.data
    h,w = np.shape(data_test)[:2]
    #read original radii
    rix = map_i.meta['rsun_obs']/map_i.meta['cdelt1']
    riy = map_i.meta['rsun_obs']/map_i.meta['cdelt2']
    #desired Sun's radii for all images
    rfx = 1900.
    rfy = 1900.
    #scale constant
    kxf = rfx/rix
    kyf = rfy/riy
    #scaling image according to scales constants
    res1 = cv2.resize(data_test,(int(w*kxf),int(h*kyf)), interpolation = cv2.INTER_CUBIC)
    hm,wm = res1.shape[:2]
    #Define new image size to9000x6000 (New radii will be kept)
    nh = 4096
    nw = 4096
    #Read original centers in X and Y
    xi = map_i.meta['crpix1']
    yi = map_i.meta['crpix2']
    #Centering the Sun's image in New image with size 4096x4096
    kx = nw/wm
    xn = kxf*xi
    ky = nh/hm
    yn = kyf*yi
    xc = 2048
    yc = 2048
    tx = xc - xn
    ty = yc - yn
    #Creating final image with new conditions
    M = np.float32([[1,0,tx],[0,1,ty]])
    new_data = cv2.warpAffine(res1,M,(nw,nh))
    new_header = map_i.meta
    new_header['cdelt1'] = map_i.meta['rsun_obs']/rfx
    new_header['cdelt2'] = map_i.meta['rsun_obs']/rfy
    new_header['crpix1'] = xc
    new_header['crpix2'] = yc
    new_fits = sunpy.map.Map(new_data,new_header)
#    mypath_fits = '/Volumes/data/Megamovie/guevara-pipeline/Stokes_recentered'
#    print(mypath_fits[0])
#    files_fits = [f for f in listdir(mypath) if isfile(join(mypath_fits, f))]
    files_fits_spl = item.split('/')
    files_fits_spl2 = files_fits_spl[-1].split('.')
    out_file = '%s.%s.%s.%s.%s'%(files_fits_spl2[0],files_fits_spl2[1],files_fits_spl2[2],files_fits_spl2[3],files_fits_spl2[4])
#    print(out_file)
    if os.path.exists('Stokes_recentered/%s_recentered.fits'%out_file):
        print("FITS already exists")
    else:
        print("Saving FITS")
        new_fits.save('/Volumes/Mnemosyne/work_20170910_flare/Stokes_recentered/%s_recentered.fits'%out_file,filetype='auto')
    del new_fits, map_i, new_header,data_test
#    print('Sending new fits to server')
#    os.system("sshpass -p '12345' scp %s_recentered.fits megamovie@swaves-TM1:/Volumes/data/Megamovie/guevara-pipeline/Stokes_recentered"%item[61:112])
#    os.remove('%s_recentered.fits'%item[61:112])
print('The code has finished')
