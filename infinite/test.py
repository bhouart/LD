import math
import numpy as np
import matplotlib.pyplot as plt

M = 10
N = 20
F = 12
c = 1

def R(x):
    part1 = 0
    for k in range(M, N, 1):
        part1 += (math.comb(N-1, k)) * x**(k-M+1) * (1-x)**(N-1-k)
        
    part2 = M * math.comb(N-1, M-1) * (1-x)**(N-M)
    
    return x**(M-1) * (part1 + part2)

def Q(x):
    return -c * (1-F/N*R(x))

def pi_d(k):
    if (k-M) >= 0:
        return k*F*c/N
    else:
        return 0
    
def pi_c(k):
    return pi_d(k) - c

def fitness_d(x):
    ret = 0
    for k in range(0, N):
        ret += math.comb(N-1, k) * x**k * (1-x)**(N-1-k) * pi_d(k)
    return ret

def fitness_c(x):
    ret = 0
    for k in range(0, N):
        ret += math.comb(N-1, k) * x**k * (1-x)**(N-1-k) * pi_c(k+1)
    return ret
    

F = 1/(R(M/N))*N

steps = np.array([i/1000 for i in range(0, 1001)])
xprime = np.zeros(1001)
xprime2 = np.zeros(1001)
fit_diff = np.zeros(1001)

for i in range(len(steps)):
    x = steps[i]
    xprime[i] = x*(1-x)*Q(x)
    xprime2[i] = x*(1-x)*(fitness_c(x) - fitness_d(x))
    fit_diff[i] = fitness_c(x) - fitness_d(x)

#plt.plot(steps, xprime)
plt.plot(steps, xprime2)
#plt.plot(steps, fit_diff)
plt.plot(steps, np.zeros(1001))
plt.show()
    