import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch


def plot_spectrum(irrep_file_names):
    colors = [ 'grey','r', 'g', 'b','yellow','orange']
    dic_colors = {0:'grey',1:'r',2:'g',3:'b',4:'tab:olive',10:'orange'}
    dic_alpha = {-1:1,1:0.5}
    dic_irreps = {"A1mP.txt":"$A1^{-+}$","A1pP.txt":"$A1^{++}$","T1mP.txt":"$T1^{-+}$","EpP.txt":"$E^{++}$","T2pP.txt":"$T2^{++}$","T1pP.txt":"$T1^{++}$"}
    omega_baryon_mass = 0.353
    axis = plt.gca()
    axis.set_ylim(2,3)
    axis.spines['top'].set_visible(False)
    axis.spines['right'].set_visible(False)
    
    for irrep_file_name in irrep_file_names:

        irrep = np.loadtxt(irrep_file_name)
        energies = irrep[:, 0]/0.353
        error_energy = irrep[:, 1]/0.353
        angular_mom = irrep[:, 2]
        parity = irrep[:, 3]
        x = [dic_irreps[irrep_file_name] for i in range(len(energies))]
        for i in range(len(energies)):
            bar_width = 0.3
            
    
            plt.bar(x[i], 2 * error_energy[i], bottom=energies[i] - error_energy[i], color=dic_colors[angular_mom[i]], alpha=dic_alpha[parity[i]], width=bar_width, edgecolor=dic_colors[angular_mom[i]])

            plt.errorbar(x[i], energies[i], yerr=error_energy[i], fmt='o', color=dic_colors[angular_mom[i]], alpha=dic_alpha[parity[i]],marker = 's')
    plt.ylabel("$\\frac{m}{m_{\Omega}}$")
    plt.xlabel("$\Lambda^{PC}$")
    plt.title("$c\\bar c$ spectrum: $16^3\\times 128$,n128,ml0743")
    plt.show()
plot_spectrum(["A1pP.txt","EpP.txt","T2pP.txt","T1pP.txt"])
        
