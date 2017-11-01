import numpy as np
import scipy.signal as sp

from .common import *

class Procressor:
    def __init__(self, kernel):
        self.kernel = kernel

        kernelSize = len(self.kernel)
        self.delay = kernelSize // 2
        self.delayed = 0
        if(kernelSize < 3):
            raise ValueError("length of kernel cannot be less than 3")
        if(kernelSize % 2 == 0):
            raise ValueError("length of kernel must be odd")
        
        self.buffer = np.zeros(kernelSize - 1, dtype = np.float64)
    
    def __call__(self, x):
        kernelSize = len(self.kernel)
        if(x is not None):
            x = np.asarray(x, dtype = np.float64)
            if(len(x) + kernelSize - 1 < 128):
                conv = np.convolve(x, self.kernel)
            else:
                conv = sp.fftconvolve(x, self.kernel)
            conv[:kernelSize - 1] += self.buffer
            self.buffer = conv[-kernelSize + 1:].copy()
            if(self.delayed < self.delay):
                conv = conv[self.delay - self.delayed:len(x)]
                self.delayed += min(len(x), self.delay - self.delayed)
                if(len(conv) > 0):
                    return conv
            else:
                return conv[:len(x)]
        else:
            if(self.delayed == 0):
                return None
            self.buffer, tmp = np.zeros(kernelSize - 1, dtype = np.float64), self.buffer[self.delay - self.delayed:kernelSize // 2]
            self.delayed = 0
            return tmp