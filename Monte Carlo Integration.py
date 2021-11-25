# This project required using Monte Carlo integration in order to calculate the volume of the unit n-ball.
# The unit n-ball is a 'ball' that is contained within 1 unit distance of the origin
# i.e. in 2 dimension, the 'volume' of the n-ball is the area of a circle with a radius of 1

import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sp

#Setup an array of the dimensions between 2 and 14
arrayLen = 14
n_array = np.linspace(2,15,arrayLen)

#solve analytic solution for plotting on graph
analytic_sol = np.pi**(n_array/2) / sp.gamma((n_array/2)+1)

#Function for Monte-Carlo solution
def MonteSol(n,npts):
    #reset count for each repetition of the function
    count = 0
    #Set up a 2d array of npts (testing points) by n (dimension number). Square all the elements of the array and take the sum of all coloumns
    array = np.random.random((npts,n))
    square_array = np.square(array)
    sum_array = np.sum(square_array,axis=1)
    #Loop over all elements of the sum array to check if the value is within one unit distance of the origin
    #If yes, increase count by 1
    for i in range(0,npts):
        if (sum_array[i]<1.0):
            count+=1
    #take the count value and multiply by 2 to the power of the dimension n. Divide this through by npts
    #npts in this case is the length of the sum array
    value = (count*(2**n))/npts
    return value

#Create array of Monte-Carlo solutions and fill it using values from the function
#n has to be a very large number in order to minimize the error at higher values of the dimension
#between the analytical solution and Monte-Carlo solution
Monte_sol = np.zeros(arrayLen)
for j in range(0,arrayLen): 
    Monte_sol[j] = MonteSol(n=int(n_array[j]),npts=1000000)

#Draw a plot of the analytical solution and Monte-Carlo solution
#Adding a legend and a title to the graph
#Save graph to file
plt.plot(n_array,analytic_sol,label='Analytical solution')
plt.plot(n_array,Monte_sol,label='Monte-Carlo solution',linestyle="dashed")
plt.title("Volume of the unit n-ball")
plt.legend()
plt.xlabel("spatial dimension")
plt.savefig("unit n-ball graph.png")
