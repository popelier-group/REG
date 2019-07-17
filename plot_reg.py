"""
plot_reg.py v0.0
L. J. Duarte, XXXXXXX, P. L. A. Popelier 

Library to easily buid plots for REG analysis. 
Check for updates at github.com/ljduarte
For details about the method, please see XXXXXXX

Please, report bugs and issues to leo.j.duarte@hotmail.com.br
coded by L. J. Duarte
"""

import matplotlib.pyplot as plt

def standardize(A):
    avgA = sum(A)/len(A)
    sqA = [(a-avgA)**2 for a in A]
    stdevA = (sum(sqA)/(len(A)-1))**(0.5)
    A = [(a-avgA)/stdevA for a in A]
    
    return A


def plot_energies(X, wfn_energies, seg=False, seg_index=[], save = False, terms=[], terms_header=[]):
    plt.figure(figsize=(16,8))
    plt.title("REG Analysis")
    plt.xlabel("Control Coordinate")
    plt.ylabel("Energy")
    plt.plot(X, wfn_energies, color = 'lightgrey', linewidth = 6 )
    plt.xlim(X[0], X[-1])
    for i in seg_index:
        plt.axvline(X[i], linestyle = 'dotted', color = 'black' )
    
    for contribution in terms:
        plt.plot(X, contribution, linewidth = 2, linestyle = '-.' )
        
    
    plt.show()