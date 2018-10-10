#!/usr/bin/env python

"""
Plotting <Q>, <I>, <U/I>, <Q/I> Stored in SWAVES server

Author: Juan Camilo Guevara Gomez
Latest version: 25 July 2018
"""
from __future__ import print_function, division

import astropy.units as u
import matplotlib
import matplotlib.pyplot as plt
import urllib.request
import numpy as np
import peakutils
import math
from matplotlib.patches import Circle
import matplotlib.patches as patches
from PIL import Image
#import cv2
from IPython.display import clear_output
from IPython import display
import time,os,sys
import datetime
import matplotlib.dates as mdates
from os import listdir
from os.path import isfile, join
import numpy.ma as ma

import matplotlib.dates as mdates

import sunpy
import sunpy.map
from sunpy.timeseries import TimeSeries
from sunpy.time import TimeRange, parse_time
from sunpy.net import hek, Fido, attrs as a
import astropy.units as u
from astropy.coordinates import SkyCoord

import pickle

import warnings
warnings.filterwarnings('ignore')

def mad(arr):
    """ Median Absolute Deviation: a "Robust" version of standard deviation.
        Indices variabililty of the sample.
        https://en.wikipedia.org/wiki/Median_absolute_deviation
    """
    arr = np.ma.array(arr).compressed() # should be faster to not use masked arrays.
    med = np.nanmedian(arr)
    return np.nanmedian(np.abs(arr - med))


# Centro del area y tamano
bxc = [2330,80]
#bxc = [2324,65]
bsz = 40
rsun = 953.3

if bxc[1] > 2048:
	nxc = bxc[1] - 2048
else:
	nxc = 2048 -bxc[1]

if bxc[0] > 2048:
	nyc = bxc[0] - 2048
else:
	nyc = 2048 -bxc[0]

nbx = [nyc, nxc]
upx = bsz*0.5*0.5016345536842105

rpx = np.sqrt(nbx[0]**2 + nbx[1]**2)
rpa = rpx*0.5016345536842105
dpa = (rpa - rsun)*0.725*100000000
upa = upx*0.725*100000000

# Pascal plot datapoints
psh_points = np.genfromtxt('hmi_055_psh.csv', delimiter=',')
#plt.plot(my_data[1:,0],my_data[1:,1])
ne2d  = np.interp(rpa/rsun,psh_points[:,0],psh_points[:,1])
ne2dn = np.interp((rpa-upx)/rsun,psh_points[:,0],psh_points[:,1])
ne2dm = np.interp((rpa+upx)/rsun,psh_points[:,0],psh_points[:,1])

#Datos de GOES
tr = TimeRange(['2017-09-10 12:30', '2017-09-10 21:30'])
results = Fido.search(a.Time(tr), a.Instrument('XRS'))

files = Fido.fetch(results)
goes = TimeSeries(files)

tci = datetime.datetime.strptime('2017-09-10 12:28:41.40','%Y-%m-%d %H:%M:%S.%f')
tcf = datetime.datetime.strptime('2017-09-10 21:22:41.30','%Y-%m-%d %H:%M:%S.%f')
xrsa = goes.data['xrsb']

client = hek.HEKClient()
flares_hek = client.search(hek.attrs.Time(tr.start, tr.end),
                           hek.attrs.FL, hek.attrs.FRM.Name == 'SWPC')

xrgoes = np.array([xrsa[i] for i in range(len(xrsa)) if tci<=xrsa.index[i]<=tcf])
xrtiempo = np.array([xrsa.index[i] for i in range(len(xrsa)) if tci<=xrsa.index[i]<=tcf])

####  HMI

mypathI = 'Stokes_IQUV_averages/'
#mypathI = 'Stokes_IQUV_rebin/'

filenamesI = [mypathI+f for f in listdir(mypathI) if isfile(join(mypathI, f))]

I_files = []
for item in filenamesI:
	if item.split('_')[5] == 'TAI.3I':
		I_files.append(item)
Q_files = []
for item in filenamesI:
	if item.split('_')[5] == 'TAI.3Q':
#		if item.split('_')[-1].split('.')[2] == 'rebin8x8':
		Q_files.append(item)
U_files = []
for item in filenamesI:
	if item.split('_')[5] == 'TAI.3U':
#		if item.split('_')[-1].split('.')[2] == 'rebin8x8':
		U_files.append(item)
V_files = []
for item in filenamesI:
	if item.split('_')[5] == 'TAI.3V':
#		if item.split('_')[-1].split('.')[2] == 'rebin8x8':
		V_files.append(item)

I_files.sort()
Q_files.sort()
U_files.sort()
V_files.sort()

#rot_pol_ang = np.radians(90. - 8.81220944)
rot_pol_ang = np.radians(90. - 8.81220944)

Int_i = []
Int_i0 = []
tiempo = []
for item in I_files:
	print('Making interval %s I'%item)
	q_i_rot = sunpy.map.Map(item)
	q_i_head = q_i_rot.meta
	#new_array = np.array(q_i_rot.data[1750:1790,3990:4030])
	new_array0 = np.array(q_i_rot.data[2028:2068,2028:2068])
	new_array = np.array(q_i_rot.data[int(bxc[0]-bsz/2):int(bxc[0]+bsz/2),int(bxc[1]-bsz/2):int(bxc[1]+bsz/2)])
	Int = np.nanmean(new_array)
	Int0 = np.nanmean(new_array0)
	Int_i.append(Int)
	Int_i0.append(Int0)
	tiempo.append(q_i_head['date-obs'])
	del q_i_rot,q_i_head,new_array,Int,Int0

Int_q = []
for item in Q_files:
	print('Making interval %s Q'%item)
	q_i_rot = sunpy.map.Map(item)
	q_i_head = q_i_rot.meta
	#new_array = np.array(q_i_rot.data[1750:1790,3990:4030])
	new_array = np.array(q_i_rot.data[int(bxc[0]-bsz/2):int(bxc[0]+bsz/2),int(bxc[1]-bsz/2):int(bxc[1]+bsz/2)])
	Int = np.nanmean(new_array)
	Int_q.append(Int)
	del q_i_rot,q_i_head,new_array,Int

Int_v = []
for item in V_files:
	print('Making interval %s V'%item)
	q_i_rot = sunpy.map.Map(item)
	q_i_head = q_i_rot.meta
	#new_array = np.array(q_i_rot.data[1750:1790,3990:4030])
	new_array = np.array(q_i_rot.data[int(bxc[0]-bsz/2):int(bxc[0]+bsz/2),int(bxc[1]-bsz/2):int(bxc[1]+bsz/2)])
	Int = np.nanmean(new_array)
	Int_v.append(Int)
	del q_i_rot,q_i_head,new_array,Int

Int_u = []
for item in U_files:
	print('Making interval %s U'%item)
	q_i_rot = sunpy.map.Map(item)
	q_i_head = q_i_rot.meta
	#new_array = np.array(q_i_rot.data[1750:1790,3990:4030])
	new_array = np.array(q_i_rot.data[int(bxc[0]-bsz/2):int(bxc[0]+bsz/2),int(bxc[1]-bsz/2):int(bxc[1]+bsz/2)])
	Int = np.nanmean(new_array)
	Int_u.append(Int)
	del q_i_rot,q_i_head,new_array,Int

#Int_q = np.array(Int_q) / np.array(Int_i)
#Int_u = np.array(Int_u) / np.array(Int_i)

Int_q_prima = np.array(Int_q)*np.cos(2*rot_pol_ang) + np.array(Int_u)*np.sin(2*rot_pol_ang)
Int_u_prima = -1.0*np.array(Int_q)*np.sin(2*rot_pol_ang) + np.array(Int_u)*np.cos(2*rot_pol_ang)

Int_q_old = Int_q
Int_u_old = Int_u

Int_q = Int_q_prima
Int_u = Int_u_prima

from scipy import stats

rmsq = np.nanmean(Int_q[0:20]) #np.sqrt(np.mean(np.square(np.array(Int_q[0:20]))))
stdq = np.std(np.array(Int_q[0:20]))
blq = np.ones(len(Int_q))*rmsq

rmsu = -1.0*np.sqrt(np.mean(np.square(np.array(Int_u[0:50]))))
stdu = np.std(np.array(Int_u[0:50]))
blu = np.ones(len(Int_u))*rmsu

rmsv = np.sqrt(np.mean(np.square(np.array(Int_v[0:6]))))
stdv = np.std(np.array(Int_v[0:6]))
new_v = []
for item in Int_v:
	if item < -1*rmsv-1*stdv:
		new_v.append(-1*rmsv-1*stdv)
	else:
		new_v.append(item)
blv = peakutils.baseline(np.array(new_v),3)


rmsi = np.sqrt(np.mean(np.square(np.array(Int_i[0:30]))))
stdi = np.std(np.array(Int_i[0:6]))
bli = np.ones(len(Int_i))*rmsi

di = []
dq = []
du = []
dv = []
for i in range(len(bli)):
	if Int_i[i]> 0: #rmsi+2*stdi:
		di.append(math.hypot(0,Int_i[i]-bli[i]))
	else:
		di.append(np.nan)
	if Int_q[i] < 100:# -1*rmsq-1*stdq:
		dq.append(math.hypot(0,Int_q[i]-blq[i]))
		du.append(math.hypot(0,Int_u[i]-blu[i]))
		dv.append(math.hypot(0,Int_v[i]-blv[i]))
	else:
		dq.append(np.nan)
		du.append(np.nan)
		dv.append(np.nan)


dq_di = np.array(dq)/np.array(di)
du_di = np.array(du)/np.array(di)
dv_di = np.array(dv)/np.array(di)

dq_di[0:34]  = np.nan
dq_di[56:] = np.nan
du_di[0:34]  = np.nan
du_di[56:] = np.nan
dv_di[0:34]  = np.nan
dv_di[56:] = np.nan


pol = np.sqrt(dq_di**2+du_di**2+dv_di**2)
n_e = ((np.array(di)/np.array(Int_i0))/ne2d)/dpa #/ 0.354 #27e8#/2.2e8
n_emin = ((np.array(di)/np.array(Int_i0))/ne2dn)/(dpa-upa)#/ 0.354
n_emax = ((np.array(di)/np.array(Int_i0))/ne2dm)/(dpa+upa)#/ 0.354

#n_e = (dq_di/ne2d*0.3)/dpa #27e8#/2.2e8
#n_emin = (dq_di/ne2dn)/(dpa-upa)
#n_emax = (dq_di/ne2dm)/(dpa+upa)


n_e[0:34]  = np.nan
n_e[56:] = np.nan

n_emin[0:34]  = np.nan
n_emin[56:] = np.nan

n_emax[0:34]  = np.nan
n_emax[56:] = np.nan

hfmt = mdates.DateFormatter('%H:%M:%S')

tiempo = [datetime.datetime.strptime(i,'%Y-%m-%dT%H:%M:%S.%f') for i in tiempo]
#dates = matplotlib.dates.date2num(list_of_datetimes)


f, axarr = plt.subplots(7, sharex=True, figsize=(6,8))

axarr[0].plot(tiempo,bli,color='grey',drawstyle='steps-mid')
axarr[0].plot(tiempo,Int_i,color='k',drawstyle='steps-mid',label='<I>')
axarr[0].set_ylabel('Stokes I')
axarr[0].yaxis.set_major_locator(plt.MaxNLocator(4))
#axarr[0].set_yticks([170,180,190,200,210])

ax2 = axarr[0].twinx()  # instantiate a second axes that shares the same x-axis

ax2.set_ylabel('GOES flux', color='saddlebrown')  # we already handled the x-label with ax1

ax2.plot(xrtiempo, xrgoes, color='saddlebrown')
ax2.set_yscale('log')
ax2.tick_params(axis='y', labelcolor='saddlebrown')
ax2.axvline(parse_time(flares_hek[2].get('event_peaktime')))
ax2.axvspan(parse_time(flares_hek[2].get('event_starttime')),
            parse_time(flares_hek[2].get('event_endtime')),
            alpha=0.2)
ax2.set_yticks([10e-7,10e-6,10e-5,10e-4])

axarr[1].plot(tiempo,blq,color='grey',drawstyle='steps-mid',linestyle='-')
axarr[1].plot(tiempo,Int_q,color='red',drawstyle='steps-mid',label='<Q>')
axarr[1].plot(tiempo,-1.5 + np.array(Int_u),color='darkblue',drawstyle='steps-mid',label='<U>')
axarr[1].plot(tiempo,-2.5 + np.array(Int_v),color='darkgreen',drawstyle='steps-mid',label='<V>')
axarr[1].errorbar(tiempo[1], 1, yerr=stdq,color='darkred',capsize=3)
axarr[1].text(.025,.7,"Q'",horizontalalignment='center',color='red',transform=axarr[1].transAxes)
axarr[1].text(.025,.3,"U'",horizontalalignment='center',color='darkblue',transform=axarr[1].transAxes)
axarr[1].text(.025,.1,"V",horizontalalignment='center',color='darkgreen',transform=axarr[1].transAxes)
axarr[1].set_ylabel("Stokes Q',U',V")
axarr[1].set_ylim(-3.25,2)

axarr[2].plot(tiempo,dq_di,color='red',drawstyle='steps-mid',label='<Q/I>')
axarr[2].set_ylabel(r'<$\Delta Q/ \Delta I$>')
axarr[2].set_ylim(-.025,0.3)
axarr[2].fill_between(tiempo, dq_di- mad(dq_di), dq_di+mad(dq_di),color='red',alpha=0.2,step='mid')

axarr[3].plot(tiempo,du_di,color='darkblue',drawstyle='steps-mid',label='<U/I>')
axarr[3].set_ylabel(r'<$\Delta U/ \Delta I$>')
axarr[3].set_ylim(-.025,0.3)
axarr[3].fill_between(tiempo, du_di-mad(du_di), du_di+mad(du_di),color='darkblue',alpha=0.2,step='mid')

axarr[4].plot(tiempo,dv_di,color='darkgreen',drawstyle='steps-mid',label='<V/I>')
axarr[4].set_ylabel(r'<$\Delta V/ \Delta I$>')
axarr[4].set_ylim(-.025,0.3)
axarr[4].fill_between(tiempo, dv_di-mad(dv_di), dv_di+mad(dv_di),color='darkgreen',alpha=0.2,step='mid')

axarr[5].plot(tiempo,pol,color='darkorange',drawstyle='steps-mid',label='P')
axarr[5].set_ylabel('Polarization')
axarr[5].set_ylim(0,0.3)

axarr[6].plot(tiempo,n_e,color='indigo',drawstyle='steps-mid',label='P')
axarr[6].set_ylabel(r'n$_e$[cm$^{-3}$]')
axarr[6].set_yscale('log')#axarr[6].set_ylim(0,0.3)
axarr[6].fill_between(tiempo, n_emax, n_emin,color='indigo',alpha=0.2,step='mid')
axarr[6].set_xlabel('Time (UTC)')
axarr[6].set_ylim(1e9,1e13)
#axarr[6].yaxis.set_major_locator(plt.MaxNLocator(3))

f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.1, hspace=0.1)

plt.gcf().autofmt_xdate()
myFmt = mdates.DateFormatter('%H:%M')
plt.gca().xaxis.set_major_formatter(myFmt)
plt.xticks(rotation=0)
plt.tight_layout()
plt.savefig('graph_r_%3d.pdf'%(100*(rpa/rsun)),format='pdf',dpi=600)
#plt.show()


hmi_data_avg_area = { "time": tiempo, "Int_i": Int_i, "Int_q": Int_q, "Int_u": Int_u, "Int_v": Int_v, "dq_di": dq_di, "du_di": du_di, "dv_di": dv_di, "pol":pol}
f=open( "save.p", "wb" )
pickle.dump( hmi_data_avg_area, f)
f.close()

plt.close()

for i in range(len(tiempo)):
    print(i, tiempo[i], ((np.array(di[i])/np.array(Int_i0[i]))/ne2d)/dpa, ne2d/(np.array(di[i])/np.array(Int_i0[i])), dpa)

#plt.plot(tiempo,Int_q,'r',tiempo,-1*np.array(Int_q_old),'b',tiempo,Int_u,'g',tiempo,-1*np.array(Int_u_old),'k')
#plt.show()


idx = 38

inmap = sunpy.map.Map(I_files[idx])

bottom_left = SkyCoord(0.5*nxc*u.arcsec - 40*u.arcsec, -0.5*nyc*u.arcsec - 30*u.arcsec,frame=inmap.coordinate_frame)
top_right = SkyCoord(0.5*nxc*u.arcsec + 30*u.arcsec, -0.5*nyc*u.arcsec + 30*u.arcsec,frame=inmap.coordinate_frame)

submap = inmap.submap(bottom_left, top_right)

# Create a new matplotlib figure, larger than default.
fig = plt.figure(figsize=(12,6))

# Add a first Axis, using the WCS from the map.
ax1 = fig.add_subplot(1,2,1, projection=inmap)

# Plot the Map on the axes with default settings.
inmap.plot()
inmap.draw_rectangle(bottom_left, 60* u.arcsec, 60* u.arcsec)

# Draw a box on the image
ax2 = fig.add_subplot(1,2,2, projection=submap)

bottom_left = SkyCoord(0.5*nxc*u.arcsec - 10*u.arcsec, -0.5*nyc*u.arcsec - 10*u.arcsec,frame=submap.coordinate_frame)
submap.plot(vmax=500)
submap.draw_rectangle(bottom_left, 20* u.arcsec, 20* u.arcsec)

plt.show()

print("frame\t         Date Obs\t \t n_e\t\t N_e\t\t L")

for i in range(len(tiempo)):
    print("%3d\t%s\t%.2e\t%.2e\t%.2e"%(i, tiempo[i], ((np.array(di[i])/np.array(Int_i0[i]))/ne2d)/dpa,
    ne2d/(np.array(di[i])/np.array(Int_i0[i])), dpa))

#plt.plot(tiempo,Int_q,'r',tiempo,-1*np.array(Int_q_old),'b',tiempo,Int_u,'g',tiempo,-1*np.array(Int_u_old),'k')
#plt.show()
