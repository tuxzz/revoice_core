import numpy as np
import pylab as pl
from revoice import *
from revoice.common import *
import pyrevoice as p
import gc

#w, sr = loadWav("voices/renri_i_A3.wav")
#w, sr = loadWav("voices/chihaya_01.wav")
w, sr = loadWav("voices/yuri_orig.wav")
prefilter = True

print("Py pYIN...")
pyinProc = p.pyin.Processor(sr, prefilter = prefilter)
obsProbList = pyinProc(w)

print("C rtPYIN...")
x = w
nX = len(x)
rtpyinProc = rtpyin.Processor(sr, prefilter = prefilter)
nHop = getNFrame(nX, rtpyinProc.hopSize)
obsProbList_c = []
iInHop = 0
iOutHop = 0
while(True):
    data = x[iInHop * rtpyinProc.hopSize:(iInHop + 1) * rtpyinProc.hopSize]
    if(len(data) == 0):
        data = None
    out = rtpyinProc(data)
    if(out is not None):
        obsProbList_c.append(out)
        iOutHop += 1
    elif(data is None):
        break
    iInHop += 1
if(iOutHop != nHop):
    print("nHop mismatch(expected %d, got %d)" % (nHop, iOutHop))
    exit(1)
del x, rtpyinProc

print("C...")
x = w
nX = len(x)
nHop = getNFrame(nX, pyinProc.hopSize)

rtmonopitchProc = rtmonopitch.Processor(*rtmonopitch.parameterFromPYin(pyinProc))
f0List = np.zeros(nHop)
for iHop in range(nHop):
    frame = getFrame(x, iHop * rtmonopitchProc.hopSize, 2 * rtmonopitchProc.hopSize)
    out = rtmonopitchProc(frame, obsProbList_c[iHop])
    if(out is None):
        continue

    nOut = len(out)
    f0List[iHop - len(out) + 1:iHop + 1] = out

del rtmonopitchProc
gc.collect()
rvExitCheck()

print("Python...")
monopitchProc = p.monopitch.Processor(*p.monopitch.parameterFromPYin(pyinProc))
f0List_o = monopitchProc(w, obsProbList)

t = np.arange(len(f0List)) * pyinProc.hopSize / sr
#t = np.arange(len(f0List))
pl.plot(t, f0List, label = "c")
pl.plot(t, f0List_o, label = "py")
pl.legend()
pl.show()