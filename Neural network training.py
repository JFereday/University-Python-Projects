import numpy as np
import matplotlib.pyplot as plt

# I1 ---W14
#   ----W15---   H4
#   ----W24--
# I2                           O6
#   ---W25---    H5
#   ----W34
# I3-------W35

#Weight table = [W14,W15,W24,W25,W34,W35,W46,W56] (len:8)

weight = 2.0*np.random.random(8)-1.0
new_weight = []


#Truth table = [[I1],[I2],[I3]]
inode = np.array([[0,0,0,0,1,1,1,1],[0,0,1,1,0,0,1,1],[0,1,0,1,0,1,0,1]])

desiredoutput = np.array([0,1,1,0,1,0,0,1])

alpha = 0.05
k = 1.0



def calculation(w):
    o = np.zeros(np.shape(inode)[1],dtype=float)
    #calculate output values for inputs
    for i in range(0,np.shape(inode)[1]):
        H4 = (inode[0][i] * w[0]) + (inode[1][i] * w[2]) + (inode[2][i] * w[4])
        H5 = (inode[0][i] * w[1]) + (inode[1][i] * w[3]) + (inode[2][i] * w[5])
        o[i] = (H4 * w[6]) + (H5 * w[7])
        #print(value)
        
    err = 0.5*sum((act(o)-act(desiredoutput))**2)
    #print(err)
    return o, err


def total_error(checkarray,output):
    check = 0.5*np.sum((act(checkarray)-act(output))**2)           
    return check

def weightmodify(w):
    rnd = np.random.randint(0,8,2)
    w[rnd[0]], w[rnd[1]] = w[rnd[1]], w[rnd[0]]
    #w[rnd[0]] = w[rnd[0]]+np.random.random(1)-1
    return w

def act(x):
    return 1.0/(1+np.exp(-k*x))

tmax = 1000
beta = 0.01

error_array = np.zeros(tmax)


for t in range(0,tmax):
    beta = beta + 1.0/tmax

    new_weight = weightmodify(w=weight)
    new_output, new_error = calculation(w=new_weight) 
    
    output, error = calculation(w=weight)
    print(np.round(new_output),np.round(output),new_error,error)
    #print(new_error,error)
    
    #want to minimize error so do sum(expect-calc)     
    if(error>new_error):
        print("shorter, accepting") 
        weight=new_weight
        error_array[t] = new_error
    else:
        prob = np.exp(beta*(act(new_error)-act(error)))
        print(prob)
        rn = np.random.random(1)
        if(rn<prob):
            print("longer, but accepting ")
            weight=new_weight
            error_array[t] = new_error
        else:
            print("longer, not accepting ")
            error_array[t] = error
            
    
plt.plot(error_array)

print(np.round(output))