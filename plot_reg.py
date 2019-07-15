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

def plot_energies(X, wfn_energies, IQA_contributions=[], seg=False, seg_index=[], save = False):
    plt.figure(figsize=(16,8))
    plt.title("REG Analysis")
    plt.xlabel("Control Coordinate")
    plt.ylabel("Energy")
    plt.plot(X, wfn_energies)
    for i in seg_index:
        plt.axvline(X[i])
    
    plt.show()