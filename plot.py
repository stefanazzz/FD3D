# read and plot results from triffy FD simulation
import pandas as pd
import os as os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import subprocess
import sys
from scipy.special import lambertw
import time
import matplotlib.colors as colors
from glob import glob

#%%#############################
#### DEFINE WORKING FIRST FOLDER
##############################
os.chdir('/home/lcbt87/FD3D/SRC_NEW/')
# read parameters of simulation from parameterisation file fpar:
f=open(file='fpar',mode='r')
fpar = f.readlines();
dx = fpar[1].split(" ");dx = float(dx[0])
dt = fpar[2].split(" ");dt = float(dt[0])
nr=len(fpar)-20
ntime = fpar[3].split(" ");ntime = int(ntime[0])
maxv=0;
#%%
def seismic_wiggle(section, gpos, spacing, dt,  normalize=False, scale=1., 
                   ranges = None, fillcolor='k', AGC=-1.,redvel=0.0, 
                   clip=False, colo=['k','r','k'], fill=True):
    """
    Plot a seismic section (numpy 2D array matrix) as wiggles.

    Parameters:

    * section :  2D array
        matrix of traces (first dimension time, second dimension traces)
    * dt : float
        sample rate in seconds (default 4 ms)
    * ranges : (x1, x2)
        min and max horizontal values (default trace number)
    * scale : float
        scale factor multiplied by the section values before plotting
    * color : tuple of strings
        Color for filling the wiggle, positive  and negative lobes.
    * normalize :
        True to normalizes each trace in the section using individual trace max/min
        data will be in the range (-0.5, 0.5) zero centered

    .. warning::
        Slow for more than 200 traces, in this case decimate your
        data or use ``seismic_image``.

    """
    import numpy as np
    import matplotlib.pyplot as pyplot
    
    def max_rolling1(a, window,axis =1):
        shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
        strides = a.strides + (a.strides[-1],)
        rolling = np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)
        return np.max(rolling,axis=axis)
    def rms_rolling(a, window):
        asq=a**2
        rrms=np.zeros(len(a))
        for i in range(window//2,len(a)-(window-1)//2):
            rva=np.sum(asq[i-window//2:i+(window-1)//2])
            rva=np.sqrt(rva/window)
            rrms[i]=rva
        return rrms
    
    npts, ntraces = section.shape  # time/traces
    if ntraces < 1:
        raise IndexError("Nothing to plot")
    if npts < 1:
        raise IndexError("Nothing to plot")
    t = np.linspace(0, dt*npts, npts)
    amp = 1.  # normalization factor
    gmin = 0.  # global minimum
    toffset = 0.  # offset in time to make 0 centered
    pyplot.ylim(max(t), 0)
    if ranges is None:
        ranges = (0, ntraces)
    x0, x1 = ranges
    # horizontal increment
    dx = spacing*float((x1-x0)/ntraces)
    slowness=0.0;tlab='t (s)'
    if redvel != 0.0: slowness = 1./redvel; tlab='t - x/%d (s)'% redvel
    pyplot.xlim(x0, x1)
    if AGC > 0:
        nlen=int(AGC/dt)
        if nlen < 1:
            print('ERROR: length of AGC window is too small (<dt)!')
            sys.exit(1)
    col=colo[0]
    for i, trace in enumerate(section.transpose()):
        if len(colo) == 2 and i >= ntraces/2: col=colo[1]
        if len(colo) == 3 and i >= ntraces/3: col=colo[1]
        if len(colo) == 3 and i >= 2*ntraces/3: col=colo[2]
        if AGC > 0:
            trace_sq=trace**2
            trace_rrms=np.convolve(trace_sq, np.ones(nlen), 'same') / float(nlen)
            tr1=np.where(trace_rrms > 0,np.sqrt(trace_rrms), 1.0)
            tr0=trace/tr1
        else:
            tr0=trace
        if normalize:
            gmax = max(tr0.min(), tr0.max(), key=abs)
            gmin = tr0.min()
            amp=abs(gmax-gmin)
            if (amp==0): 
                if (gmax == 0): amp=1
                else: amp=abs(gmax)
            toffset = 0.0
        # tr = (((trace-gmin)/amp)-toffset)*scale*dx
        tr = ((tr0/amp)-toffset)*scale*dx
        if clip:
            for j in range(len(tr)):
                if abs(tr[j]) > dx/2: tr[j]=np.sign(tr[j])*dx/2
#         if mog[0] == i+1: 
#             x = mog[1]
#         else: 
#             x = x0+float(i)*dx  # x positon for this trace
#            
        x=gpos[i]
        t0=t-x*slowness
        pyplot.plot(x+tr*0.0, t0, 'lightgrey')
        pyplot.plot(x+tr, t0, col)
        if fill==True: pyplot.fill_betweenx(t0, x, x+tr, tr > 0 , color=fillcolor)
    pyplot.xlim(left=0,right=x+gpos[0])
    pyplot.ylim(top=-20*dt)
    pyplot.ylabel(tlab)
    pyplot.xlabel('position (dx units)')

#%%################################
files = glob('./RES_plain/tracex*')
traces=[];pos=[]
for j, file in enumerate(files):
    trace=pd.read_csv(file).to_numpy().flatten()
    traces.append(trace)
    pos.append( float((fpar[20+j].split()[0]))/dx )
traces_plain=np.array(traces).T
files = glob('./RES_sphere/tracex*')
traces=[];pos=[]
for j, file in enumerate(files):
    trace=pd.read_csv(file).to_numpy().flatten()
    traces.append(trace)
    pos.append( float((fpar[20+j].split()[0]))/dx )
traces_sphere=np.array(traces).T
#%%
traces_both=np.concatenate((traces_plain,traces_sphere),axis=1)
pos=np.array(pos)
pos_both=np.concatenate((pos,pos),axis=0)
fig,ax=plt.subplots(num=2,clear=True)
seismic_wiggle(traces_both, pos_both, dx, dt, normalize=True, scale=.1, AGC=0.00, 
               redvel=0, clip=False, colo=['k','b'], fill=True)
ax.set_title('sphere and plain')
fig.savefig('wiggles_both.pdf')
#%%
pos=np.array(pos);pos_all=np.concatenate((pos,pos,pos),axis=0)
traces_all=np.concatenate((traces_plain,traces_sphere,traces_sphere-traces_plain),axis=1);
fig,ax=plt.subplots(num=3,clear=True)
seismic_wiggle(traces_all, pos_all, dx, dt, normalize=True, scale=.1, AGC=0.00, 
               redvel=0, clip=False, colo=['lightgray','lightblue','r'], fill = False)
ax.set_title('difference')
fig.savefig('wiggles_diff.pdf')

#%%#######################################
#### READ AND PLOT MOVIE
########################################
bobo = open('./iter','r')
iter= int(bobo.read().split()[2]) - 1;
iter = int(iter / 100.0) * 100
fig,ax=plt.subplots(num=1,clear=True)
maxfinal=1.0
for i in range(100,iter+100,100):
    print(i)
    inf = './RES_sphere/rdi1'+str(i).zfill(5) # current or final iteration
    fld = pd.read_csv(inf,header=None,sep='\s+', skipinitialspace=True)
    fld = fld.to_numpy()
    maxv=np.amax(np.abs(fld)); 
    if i == 600: maxfinal=maxv
    if i > 600: maxv = maxfinal
    ax.cla()
    ax.imshow(fld, norm=colors.SymLogNorm(
        linthresh=0.03, linscale=1.0,vmin=-maxv, vmax=maxv))
    #ax.imshow(fld,vmin=-maxv,vmax=maxv)
    ax.set_title(str(i).zfill(5)+' dt')
    fig.canvas.draw()
    fig.savefig('./'+str(i).zfill(5)+'.png')
    time.sleep(0.4)
    fig.canvas.flush_events()