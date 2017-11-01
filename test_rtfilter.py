import numpy as np
import scipy.signal as sp
import pylab as pl
from revoice import *
from revoice.common import *
import pyrevoice as p
import gc

for kernelSize in (3, 5, 7, 9, 111, 129, 199, 255):
    cKernel = rtfilter.firwinSingleBand(kernelSize, 0.0, 2500.0, "blackman", 22050.0)
    pyKernel = sp.firwin(kernelSize, 2500.0, window = "blackman", nyq = 22050.0)
    if((np.abs(cKernel - pyKernel) > 1e-5).any()):
        print("firwinSingleBand test failed at kernel size %d" % kernelSize)
        exit(1)
    kernel = np.random.uniform(0.0, 1.0, kernelSize)
    rtfilterProc = rtfilter.Procressor(kernel, 4097)
    pyfilterProc = p.rtfilter.Procressor(kernel)
    for dataSize in (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 128, 777, 999, 1024, 4096, 4097):
        data = np.random.uniform(0.0, 1.0, dataSize)
        for m in range(1, 10):
            rtOut = np.zeros(0)
            i = 0
            delayed = 0
            while(i < dataSize - dataSize % m):
                out = rtfilterProc(data[i:i + m])
                if(out is not None):
                    rtOut = np.concatenate((rtOut, out))
                else:
                    delayed += 1
                del out
                i += m
            if(m == 1 and delayed != rtfilterProc.delay and dataSize >= rtfilterProc.delay):
                print("Delay mismatch expected %d, real %d @ kernelSize = %d, dataSize = %d, m = %d" % (rtfilterProc.delay, delayed, kernelSize, dataSize, m))
                exit(1)
            if(i < dataSize):
                out = rtfilterProc(data[i:])
                if(out is not None):
                    rtOut = np.concatenate((rtOut, out))
            delayed = rtfilterProc(None)
            rtOut = np.concatenate((rtOut, delayed))

            pyOut = np.zeros(0)
            i = 0
            delayed = 0
            while(i < dataSize - dataSize % m):
                out = rtfilterProc(data[i:i + m])
                if(out is not None):
                    pyOut = np.concatenate((pyOut, out))
                else:
                    delayed += 1
                del out
                i += m
            if(m == 1 and delayed != rtfilterProc.delay and dataSize >= rtfilterProc.delay):
                print("Delay mismatch expected %d, real %d @ kernelSize = %d, dataSize = %d, m = %d" % (rtfilterProc.delay, delayed, kernelSize, dataSize, m))
                exit(1)
            if(i < dataSize):
                out = rtfilterProc(data[i:])
                if(out is not None):
                    pyOut = np.concatenate((pyOut, out))
            delayed = rtfilterProc(None)
            pyOut = np.concatenate((pyOut, delayed))
            
            if(rtOut.shape != pyOut.shape):
                print("Shape mismatch rtOut = %s, out = %s @ kernelSize = %d, dataSize = %d, m = %d" % (str(rtOut.shape), str(pyOut.shape), kernelSize, dataSize, m))
                exit(1)
            if((np.abs(rtOut - pyOut) < 1e-10).all()):
                print("Test passed @ kernelSize = %d, dataSize = %d, m = %d" % (kernelSize, dataSize, m))
            else:
                print("Wrong value @ kernelSize = %d, dataSize = %d, m = %d" % (kernelSize, dataSize, m))
                pl.plot(rtOut - pyOut)
                pl.show()
                exit(1)
    del rtfilterProc
gc.collect()
rvExitCheck()
print("Everything passed")