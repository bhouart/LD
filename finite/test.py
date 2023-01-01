import math
import numpy as np
import matplotlib.pyplot as plt



class Finite():
    def __init__(self, Z, M, N, F, c, beta) -> None:
        self.Z = Z        
        self.M = M
        self.N = N
        self.F = F
        self.c = c
        self.beta = beta

    def pi_d(self, k):
        if (k-self.M) >= 0:
            return k*self.F*self.c/self.N
        else:
            return 0
    
    def pi_c(self, k):
        return self.pi_d(k) - self.c

    # Average fitness
    def f_c(self, k):
        res = 0
        for j in range(self.N):
            if k > 0 and self.Z-k > 0:
                res += math.comb(k-1, j)*math.comb(self.Z-k, self.N-j-1)*self.pi_c(j+1)

        return (1/math.comb(self.Z-1, self.N-1))*res 

    def f_d(self, k):
        res = 0
        for j in range(self.N):
            if k > 0 and self.Z-k-1 > 0:
                res += math.comb(k,j)*math.comb(self.Z-k-1,self.N-j-1)*self.pi_d(j)
        return (1/math.comb(self.Z-1,self.N-1))*res


    def transition_probabilities(self, k, coeff = 1):
        return (k/self.Z) * ((self.Z-k)/self.Z) * (1/(1 + math.exp(coeff*self.beta*(self.f_c(k) - self.f_d(k)))))


    def g1(self, k):
        #return self.transition_probabilities(k, -1) - self.transition_probabilities(k, 1)
        return (k/self.Z)*((self.Z-k)/self.Z)*math.tanh((self.beta/2)*(self.f_c(k) - self.f_d(k)))


    def g2(self, k):
        return (k/self.Z)*((self.Z-k)/self.Z)*math.tanh((self.beta/2)*(self.f_c_minus_f_d(k)))


    def fermi(self, f_a, f_b):
        return 1/(1+math.exp(-self.beta*(f_a - f_b)))


    def f_c_minus_f_d(self, k):
        return self.c*((self.F/self.N)*(1-((self.N-1)/(self.Z-1)))-1)



def Fig2_plot(z):
    finite = Finite(z, 0, 10, 12, 1, 0.9)
    steps = np.array([i/z for i in range(0, z+1)])
    y_value = np.zeros(z+1)

    for k in range(len(steps)):
        y_value[k] = finite.g2(k)

    plt.plot(steps, y_value, label = str(z))
    

def Figure2():
    
    z_array = np.array([20, 55, 640, 20000])
    for z in z_array:
        Fig2_plot(z)
    plt.legend()
    plt.show()


def Fig3_plot(z, f):
    finite = Finite(z, 5, 10, f, 1, 0.9)
    steps = np.array([i/z for i in range(0, z+1)])
    y_value = np.zeros(z+1)

    for k in range(len(steps)):
        y_value[k] = finite.g1(k)

    plt.plot(steps, y_value, label = str(z))

def Figure3():
    f_array = [12, 8]
    z_array = [10,20,40,100,500,100000]
    for f in f_array:
        for z in z_array:
            Fig3_plot(z, f)
        plt.legend()
        plt.show()

Figure2()
Figure3()
    