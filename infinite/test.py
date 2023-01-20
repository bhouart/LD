import math
import numpy as np
import matplotlib.pyplot as plt

M = 10
N = 20
F = 12
c = 1
steps = np.array([i/1000 for i in range(0, 1001)])
nbr = 1


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
     
        


def replicator_equation_multiples():
    
    xprime = np.zeros(1001)
    for i in range(len(steps)):
        x = steps[i]
        #xprime[i] = x*(1-x)*Q(x)
        xprime[i] = x*(1-x)*(fitness_c(x) - fitness_d(x))
    return xprime

"""
Cases
0 : F/N < l*
1 : F/N = l*
2 : F/N > l* & F/N < 1
3 : F/N > 1
"""
def fit_diff_cases(case):
    fit_diff = np.zeros(1001)    
    roots = []
    for i in range(len(steps)):
        x = steps[i]
        fit_diff[i] = Q(x)
        
        if i > 0 and fit_diff[i] * fit_diff[i-1] <= 0:
            roots.append([x, 0])
            
    plt.plot(steps, fit_diff, zorder=1)
    plt.plot(steps, np.zeros(1001), linestyle="dotted", zorder=1)


    plt.scatter(0,0, c="black", zorder=2)
    if case < 3:
        plt.scatter(1,0, c="white", edgecolor="black", zorder=2)
    
    
    if case  == 1:
        plt.scatter(roots[0][0], roots[0][1], color="grey", edgecolors="black", zorder=2)
    
    elif  case == 2:
        plt.scatter(roots[0][0], roots[0][1], color="white", edgecolors="black", zorder=2)
        plt.scatter(roots[1][0], roots[1][1], color="black", edgecolors="black", zorder=2)
    
    elif case == 3:
        plt.scatter(roots[0][0], roots[0][1], color="white", edgecolors="black", zorder=2)
        plt.scatter(1,0, c="black", edgecolor="black", zorder=2)
    
    
    plt.xlabel("Proportion of cooperators")
    plt.ylabel("Fitness difference")
    plt.show()  


def rx_multipleF():
    f_values = (8, 8.84, 10, 25)
    r = np.zeros(1001)
    for i in range(len(steps)):
        r[i] = R(steps[i])
    plt.plot(steps, r, label=f'M = {M}, N = {N}')
    
    for f in f_values:
        plt.plot(steps, np.array([20/f for i in range(1001)]), linestyle="dotted", zorder=1, label=f'F/N = {f/N}')
    
    plt.scatter(M/N, R(M/N), color="black", edgecolors="black", zorder=2)
    plt.xlabel("Proportion of cooperators")
    plt.ylabel("R(x)")
    plt.legend(loc="lower right")
    plt.show()



if nbr == 1:
    combos = [(6, 14), (10, 12), (14,8)]
    curves = []
    for M, F in combos:
        curves.append(replicator_equation_multiples())
        plt.plot(steps, curves[-1], label=f'M = {M}, F = {F}')
    
    plt.plot(steps, np.zeros(1001), linestyle="dotted", zorder=1)
    plt.xlabel("Proportion of cooperators")
    plt.ylabel("Replicator equation")
    plt.legend(loc="lower right")
    plt.show()
    

elif nbr == 2:
    lambda_star = 1/R(M/N)
    f_values = (6, lambda_star*N, 12, 25) # 6 8.842 12 25
    for i in range(len(f_values)):
        F = f_values[i]
        fit_diff_cases(i)

elif nbr == 3:
    rx_multipleF()
