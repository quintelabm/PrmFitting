import numpy as np
import chaospy as cp
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

plt.show()

v = cp.Anglit(scale=0.02,shift=9)

#v = cp.Bradford(1, 0.01, 0.2)
v = cp.Beta(1.8, 2, 0.01, 0.2)

i = 0
num = 1000
samp = np.zeros(num)
while i<num:
    samp[i] = v.sample()
    i = i +1
    

num_bins = 50
plt.figure(0)
plt.title("pdf-delta")
plt.hist(samp, num_bins, facecolor='blue', alpha=0.5)


#v = cp.Bradford(10, 0.9995, 0.9999)
v = cp.GeneralizedHalfLogistic(shape=1, scale=0.0003, shift=0.9996)
i = 0
num = 5000
samp = np.zeros(num)
while i<num:
    samp[i] = v.sample()
    i = i +1
    

num_bins = 50
plt.figure(1)
plt.title("pdf-epsilon")
plt.hist(samp, num_bins, facecolor='green', alpha=0.5)


v = cp.Bradford(6, 3, 5)
i = 0
num = 5000
samp = np.zeros(num)
while i<num:
    samp[i] = v.sample()
    i = i +1
    

num_bins = 50
plt.figure(2)
plt.title("pdf-rho")
plt.hist(samp, num_bins, facecolor='purple', alpha=0.5)

#C
v = cp.GeneralizedHalfLogistic(shape=1, scale=5, shift=5)
i = 0
num = 5000
samp = np.zeros(num)
while i<num:
    samp[i] = v.sample()
    i = i +1
    

num_bins = 50
plt.figure(3)
plt.title("pdf-c")
plt.hist(samp, num_bins, facecolor='red', alpha=0.5)

#K
v = cp.GeneralizedHalfLogistic(shape=0.8, scale=4, shift=2)
i = 0
num = 5000
samp = np.zeros(num)
while i<num:
    samp[i] = v.sample()
    i = i +1
    

num_bins = 50
plt.figure(4)
plt.title("pdf-K")
plt.hist(samp, num_bins, facecolor='black', alpha=0.5)
