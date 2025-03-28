import matplotlib.pyplot as plt
import numpy as np

def plot_amplitude(path_data):
    x_values,y_values,err_values = [],[],[]
    with open(path_data, 'r') as f:
        data = f.readlines()
        for i,d in enumerate(data):
            if i != 0 and d != '\n':
                f = d.split('||')
                x = float(f[0])
                ys = f[1].strip().split('  ')
                x = float(x)
                
                ys = ys.strip().split('  ')
                y = float(ys[0])
                err  = float(ys[1])
                x_values.append(x)
                y_values.append(y)
                err_values.append(err)
    plt.errorbar(x_values, y_values, yerr=err_values, fmt='o')
    plt.plot(x_values, y_values, 'r-',linewidth=0.5)
    plt.xlabel('Ecm')
    plt.ylabel('Amplitude')
    plt.show()
def plot_amplitudes(paths_data,label,J,P):
    fig,ax = plt.subplots()
    colors = ['red','blue','green','black','orange','purple','brown','pink','gray','cyan']
    for i,path in enumerate(paths_data):
        print(i)
        x_values,y_values,err_values = [],[],[]
        with open(path, 'r') as f:
            data = f.readlines()
            for j,d in enumerate(data):
                if j != 0 and d != '\n':
                    
                    f = d.split('||')
                    x = float(f[0])
                    ys = f[1].strip().split('  ')
                    x = float(x)
                    
                    y = float(ys[0])
                    err  = float(ys[1])
                    x_values.append(x)
                    y_values.append(y)
                    err_values.append(err)
        #ax.plot(x_values, y_values,'o',label=label[i],color=colors[i],)
        ax.plot(x_values, y_values, '-',linewidth=0.5,label = label[i],color=colors[i])
        #ax.errorbar(x_values, y_values, yerr=err_values, fmt='o',label=label[i],color=colors[i])
        ax.fill_between(x_values, np.array(y_values)-np.array(err_values), np.array(y_values)+np.array(err_values), alpha=0.3,color = colors[i])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.legend()
    p = '+' if P == 1 else '-'
    j = str(J)
    plt.title(f'$J^P = {J}^{p}$')
    
    plt.xlabel('$E_{cm}$')
    plt.ylabel('$\\rho_i\\rho_j |t_{ij}|^2$')
    plt.show()
#plot_amplitudes(['Amplitudes_plot_data/plot_DDbar.dat','Amplitudes_plot_data/plot_psieta.dat','Amplitudes_plot_data/plot_psieta_DDbar.dat'],["$D\\bar{D}\\to D\\bar{D}$","$\\psi\\eta\\to\\psi\\eta$","$\\psi\\eta\\to D\\bar{D}$"],1,-1)
#paths = ['Amplitudes_plot_data/Dbar{D}-1^F_3--Dbar{D}-1^F_3.dat',"Amplitudes_plot_data/Dbar{D}-1^F_3--psieta-3^F_3.dat","Amplitudes_plot_data/psieta-3^F_3--psieta-3^F_3.dat"]
labels = ["$D\\bar{D}\\to D\\bar{D}$"]
paths = ['Amplitudes_plot_data/phase_shift.dat']
#paths = ['Data/plot.dat']
#labels = ['']
def plot_amplitudes_phase_shift(paths_data,label,J,P,ax):
    colors = ['red','blue','green','black','orange','purple','brown','pink','gray','cyan']
    for i,path in enumerate(paths_data):
        print(i)
        x_values,y_values,err_values = [],[],[]
        with open(path, 'r') as f:
            data = f.readlines()
            for j,d in enumerate(data):
                if j != 0 and d != '\n':
                    
                    f = d.split('||')
                    x = float(f[0])
                    ys = f[1].strip().split('  ')
                    x = float(x)
                    
                    y = float(ys[0])
                    err  = float(ys[1])
                    x_values.append(x)
                    y_values.append(y)
                    err_values.append(err)
        x_values = np.array(x_values)
        y_values = np.array(y_values)
        y_values = y_values*180/np.pi
        err_values = np.array(err_values)
        err_values = err_values*180/np.pi
        
        
        #ax.plot(x_values, y_values,'o',label=label[i],color=colors[i],)
        ax.plot(x_values, y_values, '-',linewidth=0.5,label = label[i],color=colors[i])
        #ax.errorbar(x_values, y_values, yerr=err_values, fmt='o',label=label[i],color=colors[i])
        ax.fill_between(x_values, np.array(y_values)-np.array(err_values), np.array(y_values)+np.array(err_values), alpha=0.3,color = colors[i])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend()
    p = '+' if P == 1 else '-'
    j = str(J)
    #plt.title(f'$J^P = {J}^{p}$')
    
    #plt.xlabel('$E_{cm}$')
    #plt.ylabel('$\\rho_i\\rho_j |t_{ij}|^2$')
