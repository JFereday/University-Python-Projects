import numpy as np
import matplotlib.pyplot as plt


#setup arrays for x and y based on data given in coursework document
x = np.linspace(0,(15*np.pi/8),16)
x_spacing = x[1]-x[0]
y = np.array([-0.2,-0.1,0.3,0.2,0.4,0.5,0.0,-0.4,-0.4,-0.2,0.1,0.2,0.2,0.1,0.1,-0.1])



#Coefficients calculated using numpy module
a = np.fft.fft(y)


#Set up array for continuous formula by creating an array with large number of points
array_length = 100
xf = np.linspace(0,2*np.pi,array_length)
yf = np.zeros(array_length)

n = len(y)
m = n//2

# loop over all x values
for i in range(0,array_length):
    #first term of Fourier transform
    yf[i] = 0.5 * a[0].real
    #summation part of Fourier transform
    for j in range(1, m - 1):
        ak = (a[j] + a[-j]).real
        bk = (a[-j] - a[j]).imag
        yf[i] = yf[i] + ak * np.cos(j * xf[i]) + bk * np.sin(j * xf[i])
    #Last term of Fourier transform
    yf[i] = yf[i] + 0.5 * a[-m].real * np.cos(m * xf[i])

#divide through by n as the fft uses a formula that produces a value n times too large
yf = yf/n

#plot figure
plt.plot(x,y,'bo',label="Numerical values")
plt.plot(xf,yf,label="FFT transformed formula")
plt.legend()
plt.title("Fast Fourier transform")
plt.xlabel("x")
plt.savefig("Fast Fourier Transform.png")
