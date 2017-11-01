import numpy as np
import pylab as pl
from revoice import *
from revoice.common import *
import pyrevoice as p
import gc

#w, sr = loadWav("voices/renri_i_A3.wav")
#w, sr = loadWav("voices/chihaya_01.wav")
w, sr = loadWav("voices/yuri_orig.wav")

print("pYIN...")
pyinProc = p.pyin.Processor(sr)
obsProbList = pyinProc(w)

print("C...")
x = w
nX = len(x)
nHop = getNFrame(nX, pyinProc.hopSize)

rtmonopitchProc = rtmonopitch.Processor(*rtmonopitch.parameterFromPYin(pyinProc))
f0List = np.zeros(nHop)
for iHop in range(nHop):
    frame = getFrame(x, iHop * rtmonopitchProc.hopSize, 2 * rtmonopitchProc.hopSize)
    out = rtmonopitchProc(frame, obsProbList[iHop])
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

zF0List = f0List.copy()
zF0List[f0List < 0.0] = 0.0

zF0List_o = f0List_o.copy()
zF0List_o[f0List_o < 0.0] = 0.0

t = np.arange(len(f0List)) * pyinProc.hopSize / sr
#t = np.arange(len(f0List))
pl.plot(t, f0List, label = "c")
pl.plot(t, f0List_o, label = "py")
pl.legend()
pl.show()