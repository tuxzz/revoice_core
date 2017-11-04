import numpy as np
import pylab as pl
from revoice import *
from revoice.common import *
import pyrevoice as p
import gc

#w, sr = loadWav("voices/renri_i_A3.wav")
w, sr = loadWav("voices/yuri_orig.wav")

nX = len(w)

print("Without prefilter...")
x = w
rtpyinProc = rtpyin.Processor(sr, prefilter = False)
nHop = getNFrame(nX, rtpyinProc.hopSize)
obsProbList = []
iInHop = 0
iOutHop = 0
while(True):
    data = x[iInHop * rtpyinProc.hopSize:(iInHop + 1) * rtpyinProc.hopSize]
    if(len(data) == 0):
        data = None
    out = rtpyinProc(data)
    if(out is not None):
        obsProbList.append(out)
        iOutHop += 1
    elif(data is None):
        break
    iInHop += 1
if(iOutHop != nHop):
    print("nHop mismatch(expected %d, got %d)" % (nHop, iOutHop))
    exit(1)
del x
f0List = p.pyin.extractF0(obsProbList)
del rtpyinProc

print("With prefilter...")
x = w
rtpyinProc = rtpyin.Processor(sr, prefilter = True)
nHop = getNFrame(nX, rtpyinProc.hopSize)
obsProbList_pf = []
iInHop = 0
iOutHop = 0
cBuffer = None
while(True):
    data = x[iInHop * rtpyinProc.hopSize:(iInHop + 1) * rtpyinProc.hopSize]
    if(len(data) == 0):
        data = None
    out = rtpyinProc(data)
    if(out is not None):
        if(cBuffer is None):
            cBuffer = rtpyinProc.buffer
        obsProbList_pf.append(out)
        iOutHop += 1
    elif(data is None):
        break
    iInHop += 1
if(iOutHop != nHop):
    print("nHop mismatch(expected %d, got %d)" % (nHop, iOutHop))
    exit(1)
del x
f0List_pf = p.pyin.extractF0(obsProbList_pf)

print("Python without prefilter...")
x = w
rtpyinProc = p.rtpyin.Processor(sr, prefilter = False)
nHop = getNFrame(nX, rtpyinProc.hopSize)
obsProbList_o = []
iInHop = 0
iOutHop = 0
while(True):
    data = x[iInHop * rtpyinProc.hopSize:(iInHop + 1) * rtpyinProc.hopSize]
    if(len(data) == 0):
        data = None
    out = rtpyinProc(data)
    if(out is not None):
        obsProbList_o.append(out)
        iOutHop += 1
    elif(data is None):
        break
    iInHop += 1
if(iOutHop != nHop):
    print("nHop mismatch(expected %d, got %d)" % (nHop, iOutHop))
    exit(1)
del x
f0List_o = p.pyin.extractF0(obsProbList_o)

print("Python with prefilter...")
x = w
rtpyinProc = p.rtpyin.Processor(sr, prefilter = True)
nHop = getNFrame(nX, rtpyinProc.hopSize)
obsProbList_pf_o = []
iInHop = 0
iOutHop = 0
pyBuffer = None
while(True):
    data = x[iInHop * rtpyinProc.hopSize:(iInHop + 1) * rtpyinProc.hopSize]
    if(len(data) == 0):
        data = None
    out = rtpyinProc(data)
    if(out is not None):
        if(pyBuffer is None):
            pyBuffer = rtpyinProc.buffer
        obsProbList_pf_o.append(out)
        iOutHop += 1
    elif(data is None):
        break
    iInHop += 1
if(iOutHop != nHop):
    print("nHop mismatch(expected %d, got %d)" % (nHop, iOutHop))
    exit(1)
del x
f0List_pf_o = p.pyin.extractF0(obsProbList_pf_o)

if((np.abs(f0List - f0List_o) > 1e-5).any()):
    print("Test failed without prefilter, max diff = %lf" % np.max(np.abs(f0List - f0List_o)))
if((np.abs(f0List_pf - f0List_pf_o) > 1e-5).any()):
    print("Test failed with prefilter, max diff = %lf" % np.max(np.abs(f0List_pf - f0List_pf_o)))

for i, obsProb in enumerate(obsProbList):
    if(obsProbList_o[i].shape != obsProb.shape or (np.abs(obsProbList_o[i] - obsProb) > 1e-5).any()):
        print("Test failed without prefilter @ obsProb", i)
        print("  Diff:", obsProb - obsProbList_o[i])

for i, obsProb in enumerate(obsProbList_pf):
    if(obsProbList_pf_o[i].shape != obsProb.shape or (np.abs(obsProbList_pf_o[i] - obsProb) > 1e-5).any()):
        print("Test failed with prefilter @ obsProb", i)
        print("  Diff:", obsProb - obsProbList_pf_o[i])

gc.collect()
rvExitCheck()

t = np.arange(len(f0List)) * rtpyinProc.hopSize / sr
#t = np.arange(len(f0List))
pl.figure()
pl.plot(t, f0List, label = "Direct")
pl.plot(t, f0List_o, label = "Direct Py")
pl.legend()
pl.figure()
pl.plot(t, f0List_pf, label = "Prefiltered")
pl.plot(t, f0List_pf_o, label = "Prefiltered Py")
pl.legend()
pl.show()