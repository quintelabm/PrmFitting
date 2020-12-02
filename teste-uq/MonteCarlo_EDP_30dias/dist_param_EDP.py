import numpy as np
import chaospy as cp
import matplotlib.pyplot as plt

c = cp.Normal(19,1)

i = 0
num = 1000
samp = np.zeros(num)
while i<num:
    samp[i] = c.sample()
    i = i +1
    

num_bins = 50
plt.figure(0)
plt.title("c-pdf")
plt.hist(samp, num_bins, facecolor='blue', alpha=0.5)


c = cp.Normal(4,1)

i = 0
num = 1000
samp = np.zeros(num)
while i<num:
    samp[i] = c.sample()
    i = i +1
    

num_bins = 50
plt.figure()
plt.title("sigma-pdf")
plt.hist(samp, num_bins, facecolor='blue', alpha=0.5)


c = cp.Normal(0.8,0.1)

i = 0
num = 1000
samp = np.zeros(num)
while i<num:
    samp[i] = c.sample()
    i = i +1
    

num_bins = 50
plt.figure()
plt.title("mu_t-pdf")
plt.hist(samp, num_bins, facecolor='blue', alpha=0.5)


c = cp.Normal(2.8,0.5)

i = 0
num = 1000
samp = np.zeros(num)
while i<num:
    samp[i] = c.sample()
    i = i +1
    

num_bins = 50
plt.figure()
plt.title("mu_c-pdf")
plt.hist(samp, num_bins, facecolor='blue', alpha=0.5)


v = cp.GeneralizedHalfLogistic(shape=1, scale=0.009, shift=0.990)
i = 0
num = 1000
samp = np.zeros(num)
while i<num:
    samp[i] = v.sample()
    i = i +1
    

num_bins = 50
plt.figure()
plt.title("pdf-epsilon-s")
plt.hist(samp, num_bins, facecolor='green', alpha=0.5)


v = cp.Bradford(2, 0.9, 0.96)
i = 0
num = 1000
samp = np.zeros(num)
while i<num:
    samp[i] = v.sample()
    i = i +1
    

num_bins = 50
plt.figure()
plt.title("pdf-epsilon-a")
plt.hist(samp, num_bins, facecolor='purple', alpha=0.5)


c = cp.Normal(0.3,0.03)

i = 0
num = 1000
samp = np.zeros(num)
while i<num:
    samp[i] = c.sample()
    i = i +1
    

num_bins = 50
plt.figure()
plt.title("pdf-epsilon-r")
plt.hist(samp, num_bins, facecolor='blue', alpha=0.5)


v = cp.Bradford(2, 0.01, 0.3)
i = 0
num = 1000
samp = np.zeros(num)
while i<num:
    samp[i] = v.sample()
    i = i +1
    

num_bins = 50
plt.figure()
plt.title("pdf-delta")
plt.hist(samp, num_bins, facecolor='purple', alpha=0.5)


v = cp.GeneralizedHalfLogistic(shape=1, scale=9, shift=15)
i = 0
num = 1000
samp = np.zeros(num)
while i<num:
    samp[i] = v.sample()
    i = i +1
    

num_bins = 50
plt.figure()
plt.title("pdf-alpha")
plt.hist(samp, num_bins, facecolor='green', alpha=0.5)


v = cp.GeneralizedHalfLogistic(shape=1.3, scale=6, shift=5)
i = 0
num = 1000
samp = np.zeros(num)
while i<num:
    samp[i] = v.sample()
    i = i +1
    

num_bins = 50
plt.figure()
plt.title("pdf-r")
plt.hist(samp, num_bins, facecolor='green', alpha=0.5)



v = cp.GeneralizedHalfLogistic(shape=1, scale=3, shift=9.7)
i = 0
num = 1000
samp = np.zeros(num)
while i<num:
    samp[i] = v.sample()
    i = i +1
    

num_bins = 50
plt.figure()
plt.title("pdf-rho")
plt.hist(samp, num_bins, facecolor='green', alpha=0.5)