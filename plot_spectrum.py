import numpy as np
import os
import matplotlib.pyplot as plt
from colors import color_coding_file
from colors import color_coding_file_simple
from colors import color_coding_dict
import matplotlib.patches as mpatches
from colors import save_color_code_state
import sys
from colors import operator_identification
from colors import operator_identification_plus
from particle import read_particles
import grouptheory as gt
import matplotlib
import particle as p
import no_int_E_levels as no_int
from matplotlib.patches import FancyBboxPatch


def plot_spectrum(irrep_file_names):
    colors = [ 'grey','r', 'g', 'b','yellow','orange']
    dic_colors = {0:'grey',1:'r',2:'g',3:'b',4:'tab:olive',10:'orange'}
    dic_alpha = {-1:1,1:0.5}
    dic_irreps = {"GroupTheory/A1mP.txt":"$A1^{-+}$","GroupTheory/A1pP.txt":"$A1^{++}$","GroupTheory/T1mP.txt":"$T1^{-+}$","GroupTheory/EpP.txt":"$E^{++}$","GroupTheory/T2pP.txt":"$T2^{++}$","GroupTheory/T1pP.txt":"$T1^{++}$"}
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
#plot_spectrum(["GroupTheory/A1pP.txt","GroupTheory/EpP.txt","GroupTheory/T2pP.txt","GroupTheory/T1pP.txt"])
        
def plot_irrep_mass_with_get_finite_ax(spectrums,irrep_name,Ls,Ps,max_state,channel,th,path_to_file,ax):
    
        colors = ['red','blue','green','orange','purple','pink','black','brown','teal','cyan','magenta','grey','lime','olive','yellow','navy','maroon','aqua','fuchsia','silver','red','blue','green','orange','purple','pink','black','brown','yellow','cyan','magenta','grey','lime','olive','teal','navy','maroon','aqua','fuchsia','silver']
        c = channel.copy()
        C_p = '+' if channel['C_parity'] == 1 else '-'
        threshold = 0.74
        channelsDs = p.Ds(channel,threshold)
        chh = p.channels(channel,threshold)[0]
        for key in channelsDs.keys():
            chh[key] = channelsDs[key]
        dict = {}
        channels = chh.keys()
        for i,chan in enumerate(channels):
            dict[chan] = colors[i]
        dict = color_coding_dict()


        

        


        dict2 = {}
        for i,chan in enumerate(channels):
            dict2[i] = chan
        ch = list(channels)
        counts = {}
        for i,chan in enumerate(channels):
            counts[chan] = 0
        labels = []
        av_c = []
        E_levels = no_int.get_E_levels_in_flight(channel,irrep_name,Ps,Ls)
        for channel in E_levels.keys():
            Es,multiplcities = E_levels[channel][0],E_levels[channel][1]
            padding = np.array([0.001 for l in Es])
            label = "$" + channel[0] + " " + channel[1] + "$"
            chan = (channel[0],channel[1])
            for j in range(multiplcities):
                if label not in labels and min(Es) < th:
                    ax.plot(Ls,Es+j*padding,color = dict[chan],label = label)
                    labels.append(label)
                    av_c.append(chan) 
                else:
                    ax.plot(Ls,Es+j*padding,color = dict[chan])
      
            
        all_particles = p.read_particles('Particles/particles_unfl.txt')+p.read_particles('Particles/charmonium.txt')+p.read_particles('Particles/Ds.txt')

        for i,channel in enumerate(av_c):
            for particle in all_particles:
                if particle.name == channel[0]:
                    mass_1 = particle.Mass
                if particle.name == channel[1]:
                    mass_2 = particle.Mass
            mass = mass_1 + mass_2
            ax.plot(Ls,[mass for l in Ls],color = dict[channel],linestyle = '--',linewidth = 0.7)
                    
        
        ax.set_xlim(14,26)
        ax.set_ylim(0.6,th)
        
        
        
        title_name = ""
        flag = True
        for i,letter in enumerate(irrep_name):
            if i == 0:
                title_name += letter

            elif letter.isnumeric():
                title_name +=  "_"+letter + "^{"
                flag = False
            elif i == 1:
                if letter == "M":
                    let = "-"
                elif letter == "P":
                    let = "+"
                elif letter == "m":
                    let = "-"
                elif letter == "p":
                    let = "+"
                title_name += "^{"+let
                flag = False
            elif flag:
                if letter == "M":
                    let = "-"
                elif letter == "P":
                    let = "+"
                elif letter == "m":
                    let = "-"
                elif letter == "p":
                    let = "+"
                title_name += "{" + let
                flag = False
            elif i != len(irrep_name)-1:
                if letter == "M":
                    let = "-"
                elif letter == "P":
                    let = "+"
                elif letter == "m":
                    let = "-"
                elif letter == "p":
                    let = "+"
                title_name += let
            else:
                if letter == "M":
                    let = "-"
                elif letter == "P":
                    let = "+"
                elif letter == "m":
                    let = "-"
                elif letter == "p":
                    let = "+"
                title_name += let + "}"

            

        title = 'Irrep: $ '+title_name+'$' + " " + "$[" +Ps+"]$" 
        #plt.title(title)
        volumes = []
        t0s = []
        skpss = []
        max_states = []
        irreps = []
        for spectrum in spectrums:
            irreps.append(spectrum.irrep)
            volumes.append(spectrum.volume)
            t0s.append(spectrum.t0)
            skpss.append(spectrum.skips)
            max_states.append(max_state)
        irr = irreps[0]
        for i in range(1,len(irreps)):
            pass
            if irr != irreps[i]:
                pass

        for i,volume in enumerate(volumes):
            irr = irreps[i]
            t0 = t0s[i]
            max_state = max_states[i]
            skps = skpss[i]

            mass,err = spectrums[i].get_masses()
            mass_,err_ = [],[]
            for i in range(max_state):
                mass_.append(mass[i])
                err_.append(err[i])
            mass = mass_
            err = err_

            i = range(len(mass))
            xs = [volume+0.001*mas for mas in i]
            nams = [i for i in range(54)]
            for skip in skps:
                if skip in nams:
                    nams.remove(skip)
            colors = ['black' for l in range(max_state)]
            dcol = color_coding_file_simple()
            for j in range(max_state):
                name = f"{irr}/Volume_{volume}/t0{t0}/StateColorFiles/state{j}.txt"
                if os.path.exists(name):
                    with open(name) as f:
                        lines = f.readlines()
                        number = int(lines[0])
                
                    colors[j] = dcol[number]

            

            
            for i in range(len(mass)):
                if mass[i] <= th:
                    if Ps == "100":
                        m = mass[i]**2-((2*np.pi)/(xs[i]*3.444))**2
                        m = np.sqrt(m)
                    else:
                        m = mass[i]

                    ax.errorbar(xs[i],m,yerr=err[i],fmt='o',c=colors[i],markersize=3)
                    ax.text(xs[i],m,nams[i],fontsize = 8)
        xs = []
        errs = []
        ys = []
        with open(path_to_file,'r') as f:
            lines = f.readlines()
            for line in lines:
                elems = line.split("   ")
                volume = float(elems[0])
                elm2 = elems[1].split()
                energy = float(elm2[0])
                error = float(elm2[1])
                
                errs.append(error)
                xs.append(volume)
                ys.append(energy)
       
        #plt.errorbar(xs,ys,yerr=errs,color = 'gold',markersize = 3,marker = '*',linestyle = '')
        initial_L = xs[0]
        level = []
        levels = []
        levels = {}
        count = 0
        for i,L in enumerate(xs):
            if L != initial_L:
                level = ys[i]
                err = errs[i]
                initial_L = L
                count = 0
                if count not in levels.keys():
                    levels[count] = [(level,L,err)]
                else:
                    levels[count] += [(level,L,err)]
                count += 1
            else:
                level = ys[i]
                err = errs[i]
                if count not in levels.keys():
                    levels[count] = [(level,L,err)]
                else:
                    levels[count] += [(level,L,err)]
                count += 1
        for key in levels.keys():
            xs = []
            ys = []
            errs = []
            for elem in levels[key]:
                ys.append(elem[0])
                xs.append(elem[1])
                errs.append(elem[2])
                
            ax.plot(xs,ys,marker = '',linestyle = '--',color = 'gold',markersize = 3)
            #ax.errorbar(xs,ys,yerr=errs,color = 'gold',markersize = 3,marker = '*',linestyle = '')
            ax.fill_between(xs,np.array(ys)-np.array(errs),np.array(ys)+np.array(errs),color = 'gold',alpha = 0.3)
            


            


        #plt.fill_between(xs,np.array(ys)-np.array(errs),np.array(ys)+np.array(errs),color = 'gold',alpha = 0.3)

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        #plt.xlabel('Volume')
        #plt.ylabel('$a_tE_{cm}$',rotation = 0)
        irrep_name = irr 
        name = f"{irrep_name}/Spectrum_{irrep_name}.pdf"
        ax.legend(ncol = 2,loc = 'lower left',fontsize = 6)
        #plt.savefig(name)
        #plt.show()"
def plot_irrep_mass_no_sigma_grey_out(spectrums,irrep_name,Ls,Ps,max_state,channel,th,included,ax):
    
        colors = ['red','blue','green','orange','purple','pink','black','brown','teal','cyan','magenta','grey','lime','olive','yellow','navy','maroon','aqua','fuchsia','silver','red','blue','green','orange','purple','pink','black','brown','yellow','cyan','magenta','grey','lime','olive','teal','navy','maroon','aqua','fuchsia','silver']
        c = channel.copy()
        C_p = '+' if channel['C_parity'] == 1 else '-'
        threshold = 0.74
        channelsDs = p.Ds(channel,threshold)
        chh = p.channels(channel,threshold)[0]
        for key in channelsDs.keys():
            chh[key] = channelsDs[key]
        dict = {}
        channels = chh.keys()
        for i,chan in enumerate(channels):
            dict[chan] = colors[i]
        dict = color_coding_dict()


        

        


        dict2 = {}
        for i,chan in enumerate(channels):
            dict2[i] = chan
        ch = list(channels)
        counts = {}
        for i,chan in enumerate(channels):
            counts[chan] = 0
        labels = []
        av_c = []
        E_levels = no_int.get_E_levels_in_flight(channel,irrep_name,Ps,Ls)
        for channel in E_levels.keys():
            Es,multiplcities = E_levels[channel][0],E_levels[channel][1]
            padding = np.array([0.001 for l in Es])
            label = "$" + channel[0] + " " + channel[1] + "$"
            chan = (channel[0],channel[1])
            for j in range(multiplcities):
                if label not in labels and min(Es) < th:
                    ax.plot(Ls,Es+j*padding,color = dict[chan],label = label)
                    labels.append(label)
                    av_c.append(chan) 
                else:
                    ax.plot(Ls,Es+j*padding,color = dict[chan])
        
        
            
        all_particles = p.read_particles('Particles/particles_unfl.txt')+p.read_particles('Particles/charmonium.txt')+p.read_particles('Particles/Ds.txt')

        for i,channel in enumerate(av_c):
            for particle in all_particles:
                if particle.name == channel[0]:
                    mass_1 = particle.Mass
                if particle.name == channel[1]:
                    mass_2 = particle.Mass
            mass = mass_1 + mass_2
            ax.plot(Ls,[mass for l in Ls],color = dict[channel],linestyle = '--',linewidth = 0.7)
                    

        plt.xlim(14,26)
        plt.ylim(0.6,th)
        title_name = ""
        flag = True
        for i,letter in enumerate(irrep_name):
            if i == 0:
                title_name += letter

            elif letter.isnumeric():
                title_name +=  "_"+letter + "^{"
                flag = False
            elif i == 1:
                if letter == "M":
                    let = "-"
                elif letter == "P":
                    let = "+"
                elif letter == "m":
                    let = "-"
                elif letter == "p":
                    let = "+"
                title_name += "^{"+let
                flag = False
            elif flag:
                if letter == "M":
                    let = "-"
                elif letter == "P":
                    let = "+"
                elif letter == "m":
                    let = "-"
                elif letter == "p":
                    let = "+"
                title_name += "{" + let
                flag = False
            elif i != len(irrep_name)-1:
                if letter == "M":
                    let = "-"
                elif letter == "P":
                    let = "+"
                elif letter == "m":
                    let = "-"
                elif letter == "p":
                    let = "+"
                title_name += let
            else:
                if letter == "M":
                    let = "-"
                elif letter == "P":
                    let = "+"
                elif letter == "m":
                    let = "-"
                elif letter == "p":
                    let = "+"
                title_name += let + "}"

            

        title = 'Irrep: $ '+title_name+'$' + " " + "$[" +Ps+"]$" 
        #plt.title(title)
        volumes = []
        t0s = []
        skpss = []
        max_states = []
        irreps = []
        for spectrum in spectrums:
            irreps.append(spectrum.irrep)
            volumes.append(spectrum.volume)
            t0s.append(spectrum.t0)
            skpss.append(spectrum.skips)
            max_states.append(max_state)
        irr = irreps[0]
        for i in range(1,len(irreps)):
            pass
            if irr != irreps[i]:
                raise ValueError("Different irreps")

        for i,volume in enumerate(volumes):
            irr = irreps[i]
            t0 = t0s[i]
            max_state = max_states[i]
            skps = skpss[i]
            include = included[i]

            mass,err = spectrums[i].get_masses()
            mass_,err_ = [],[]
            for i in range(max_state):
                mass_.append(mass[i])
                err_.append(err[i])
            mass = mass_
            err = err_

            i = range(len(mass))
            xs = [volume+0.001*mas for mas in i]
            nams = [i for i in range(54)]
            for skip in skps:
                if skip in nams:
                    nams.remove(skip)
            colors = ['black' for l in range(max_state)]
            dcol = color_coding_file_simple()
            for j in range(max_state):
                name = f"{irr}/Volume_{volume}/t0{t0}/StateColorFiles/state{j}.txt"
                if os.path.exists(name):
                    with open(name) as f:
                        lines = f.readlines()
                        number = int(lines[0])
                
                    colors[j] = dcol[number]

            

            
            for i in range(len(mass)):
                if mass[i] <= th:
                    
                    if i not in include:
                        ax.errorbar(xs[i],mass[i],yerr=err[i],fmt='o',c='grey',markersize=3)
                    else:
                        ax.errorbar(xs[i],mass[i],yerr=err[i],fmt='o',c=colors[i],markersize=3)
                    #ax.text(xs[i],mass[i],nams[i],fontsize = 8)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        #ax.set_xlabel('Volume')
        #ax.set_ylabel('$a_tE_{cm}$',rotation = 0)
        irrep_name = irr 
        name = f"{irrep_name}/Spectrum_{irrep_name}.pdf"
        ax.legend(ncol = 2,loc = 'lower left',fontsize = 6)
        #plt.savefig(name)
        #plt.show()
def plot_irrep_mass(spectrums,irrep_name,Ls,Ps,max_state,channel,th):
    
    colors = ['red','blue','green','orange','purple','pink','black','brown','teal','cyan','magenta','grey','lime','olive','yellow','navy','maroon','aqua','fuchsia','silver','red','blue','green','orange','purple','pink','black','brown','yellow','cyan','magenta','grey','lime','olive','teal','navy','maroon','aqua','fuchsia','silver']
    c = channel.copy()
    C_p = '+' if channel['C_parity'] == 1 else '-'
    threshold = 0.74
    channelsDs = p.Ds(channel,threshold)
    chh = p.channels(channel,threshold)[0]
    for key in channelsDs.keys():
        chh[key] = channelsDs[key]
    dict = {}
    channels = chh.keys()
    for i,chan in enumerate(channels):
        dict[chan] = colors[i]
    dict = color_coding_dict()


    

    


    dict2 = {}
    for i,chan in enumerate(channels):
        dict2[i] = chan
    ch = list(channels)
    counts = {}
    for i,chan in enumerate(channels):
        counts[chan] = 0
    fig, ax = plt.subplots()
    labels = []
    av_c = []
    E_levels = no_int.get_E_levels_in_flight(channel,irrep_name,Ps,Ls)
    for channel in E_levels.keys():
        Es,multiplcities = E_levels[channel][0],E_levels[channel][1]
        padding = np.array([0.001 for l in Es])
        label = "$" + channel[0] + " " + channel[1] + "$"
        chan = (channel[0],channel[1])
        for j in range(multiplcities):
            if label not in labels and min(Es) < th:
                ax.plot(Ls,Es+j*padding,color = dict[chan],label = label)
                labels.append(label)
                av_c.append(chan) 
            else:
                ax.plot(Ls,Es+j*padding,color = dict[chan])
    
    E_levels_sigma = no_int.get_E_levels_irrep_in_flight_sigma(c,irrep_name,Ps)
    for channel in E_levels_sigma.keys():
        if channel[1] == 'f_0':
            ch2 = "\sigma"
            ch1 = channel[0]
        else:
            ch1 = "\sigma"
            ch2 = channel[0]
        chan = (ch1,ch2)

        Vs,Es,Errs,mults = E_levels_sigma[channel][0],E_levels_sigma[channel][1],E_levels_sigma[channel][2],E_levels_sigma[channel][3]
        if Vs == []:
            print("missing sigma information for ",channel)
            continue
        color = dict[chan]
        label = "$" + ch1 + " " + ch2 + "$"
        for i,V in enumerate(Vs):
            Ess = Es[i]
            Vss = [V for l in Ess]
            Erss = Errs[i]
            for j in range(mults):
                if label not in labels and min(Ess) < th:
                    padding = np.array([0.1 for l in Ess])

                    ax.errorbar(Vss+padding*j,Ess,yerr=Erss,color = color,alpha = 0.3,fmt='o',markersize = 3,label = label)
                    labels.append(label)
                    av_c.append((ch1,ch2))
                else:
                    padding = np.array([0.1 for l in Ess])
                    ax.errorbar(Vss+padding*j,Ess,yerr=Erss,color = color,alpha = 0.3,fmt='o',markersize = 3)

        
    all_particles = p.read_particles('Particles/particles_unfl.txt')+p.read_particles('Particles/charmonium.txt')+p.read_particles('Particles/Ds.txt')

    for i,channel in enumerate(av_c):
        for particle in all_particles:
            if particle.name == channel[0]:
                mass_1 = particle.Mass
            if particle.name == channel[1]:
                mass_2 = particle.Mass
        mass = mass_1 + mass_2
        ax.plot(Ls,[mass for l in Ls],color = dict[channel],linestyle = '--',linewidth = 0.7)
                

    plt.xlim(14,26)
    plt.ylim(0.6,th)
    title_name = ""
    flag = True
    for i,letter in enumerate(irrep_name):
        if i == 0:
            title_name += letter

        elif letter.isnumeric():
            title_name +=  "_"+letter + "^{"
            flag = False
        elif i == 1:
            if letter == "M":
                let = "-"
            elif letter == "P":
                let = "+"
            elif letter == "m":
                let = "-"
            elif letter == "p":
                let = "+"
            title_name += "^{"+let
            flag = False
        elif flag:
            if letter == "M":
                let = "-"
            elif letter == "P":
                let = "+"
            elif letter == "m":
                let = "-"
            elif letter == "p":
                let = "+"
            title_name += "{" + let
            flag = False
        elif i != len(irrep_name)-1:
            if letter == "M":
                let = "-"
            elif letter == "P":
                let = "+"
            elif letter == "m":
                let = "-"
            elif letter == "p":
                let = "+"
            title_name += let
        else:
            if letter == "M":
                let = "-"
            elif letter == "P":
                let = "+"
            elif letter == "m":
                let = "-"
            elif letter == "p":
                let = "+"
            title_name += let + "}"

        

    title = 'Irrep: $ '+title_name+'$' + " " + "$[" +Ps+"]$" 
    plt.title(title)
    volumes = []
    t0s = []
    skpss = []
    max_states = []
    irreps = []
    
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.xlabel('Volume')
    plt.ylabel('$a_tE_{cm}$',rotation = 0)
    #irrep_name = irr
    name = f"{irrep_name}/Spectrum_{irrep_name}.pdf"
    plt.legend(fontsize = 6,ncols = 2)
    #plt.savefig(name)
    plt.show()
def plot_elastic_phase_shift_discrete(discrete_path,ax):
    with open(file=discrete_path,mode='r') as f:
        lines = f.readlines()
        for line in lines:
            elems = line.split("  ")

            print(elems)
            ecm = elems[1].split('ecm = ')[1]
            ecm_val = float(ecm.split('+/-')[0])
            ecm_err = float(ecm.split('+/-')[1])
            print(f"ecm = {ecm_val} +/- {ecm_err}")
            phase = elems[3].split('real[delta](deg) = ')[1]
            phase_val = float(phase.split('+/-')[0])
            phase_err = float(phase.split('+/-')[1])
            irrep = elems[0]
            if "T1m" in irrep:
                color = "blue"
            if "A1" in irrep:
                color = "red"
            if "A2m" in irrep:
                color = "green"

            

            #L = int(elems[0])
            #phase = float(elems[1])
            ax.errorbar(ecm_val,phase_val,xerr = ecm_err,yerr = phase_err,fmt = '.',color = color,markersize = 1)

def plot_irrep_mass_with_get_finite_E_vs_L_ax(spectrums,irrep_name,Ls,Ps,max_state,channel,th,path_to_folder,ax):
    
        colors = ['red','blue','green','orange','purple','pink','black','brown','teal','cyan','magenta','grey','lime','olive','yellow','navy','maroon','aqua','fuchsia','silver','red','blue','green','orange','purple','pink','black','brown','yellow','cyan','magenta','grey','lime','olive','teal','navy','maroon','aqua','fuchsia','silver']
        c = channel.copy()
        C_p = '+' if channel['C_parity'] == 1 else '-'
        threshold = 0.74
        channelsDs = p.Ds(channel,threshold)
        chh = p.channels(channel,threshold)[0]
        for key in channelsDs.keys():
            chh[key] = channelsDs[key]
        dict = {}
        channels = chh.keys()
        for i,chan in enumerate(channels):
            dict[chan] = colors[i]
        dict = color_coding_dict()


        

        


        dict2 = {}
        for i,chan in enumerate(channels):
            dict2[i] = chan
        ch = list(channels)
        counts = {}
        for i,chan in enumerate(channels):
            counts[chan] = 0
        labels = []
        av_c = []
        E_levels = no_int.get_E_levels_in_flight(channel,irrep_name,Ps,Ls)
        for channel in E_levels.keys():
            Es,multiplcities = E_levels[channel][0],E_levels[channel][1]
            padding = np.array([0.001 for l in Es])
            label = "$" + channel[0] + " " + channel[1] + "$"
            chan = (channel[0],channel[1])
            for j in range(multiplcities):
                if label not in labels and min(Es) < th:
                    ax.plot(Ls,Es+j*padding,color = dict[chan],label = label)
                    labels.append(label)
                    av_c.append(chan) 
                else:
                    ax.plot(Ls,Es+j*padding,color = dict[chan])
      
            
        all_particles = p.read_particles('Particles/particles_unfl.txt')+p.read_particles('Particles/charmonium.txt')+p.read_particles('Particles/Ds.txt')

        for i,channel in enumerate(av_c):
            for particle in all_particles:
                if particle.name == channel[0]:
                    mass_1 = particle.Mass
                if particle.name == channel[1]:
                    mass_2 = particle.Mass
            mass = mass_1 + mass_2
            ax.plot(Ls,[mass for l in Ls],color = dict[channel],linestyle = '--',linewidth = 0.7)
                    
        
        ax.set_xlim(14,26)
        ax.set_ylim(0.6,th)
        
        
        
        title_name = ""
        flag = True
        for i,letter in enumerate(irrep_name):
            if i == 0:
                title_name += letter

            elif letter.isnumeric():
                title_name +=  "_"+letter + "^{"
                flag = False
            elif i == 1:
                if letter == "M":
                    let = "-"
                elif letter == "P":
                    let = "+"
                elif letter == "m":
                    let = "-"
                elif letter == "p":
                    let = "+"
                title_name += "^{"+let
                flag = False
            elif flag:
                if letter == "M":
                    let = "-"
                elif letter == "P":
                    let = "+"
                elif letter == "m":
                    let = "-"
                elif letter == "p":
                    let = "+"
                title_name += "{" + let
                flag = False
            elif i != len(irrep_name)-1:
                if letter == "M":
                    let = "-"
                elif letter == "P":
                    let = "+"
                elif letter == "m":
                    let = "-"
                elif letter == "p":
                    let = "+"
                title_name += let
            else:
                if letter == "M":
                    let = "-"
                elif letter == "P":
                    let = "+"
                elif letter == "m":
                    let = "-"
                elif letter == "p":
                    let = "+"
                title_name += let + "}"

            

        title = 'Irrep: $ '+title_name+'$' + " " + "$[" +Ps+"]$" 
        #plt.title(title)
        volumes = []
        t0s = []
        skpss = []
        max_states = []
        irreps = []
        for spectrum in spectrums:
            irreps.append(spectrum.irrep)
            volumes.append(spectrum.volume)
            t0s.append(spectrum.t0)
            skpss.append(spectrum.skips)
            max_states.append(max_state)
        irr = irreps[0]
        for i in range(1,len(irreps)):
            pass
            if irr != irreps[i]:
                pass

        for i,volume in enumerate(volumes):
            irr = irreps[i]
            t0 = t0s[i]
            max_state = max_states[i]
            skps = skpss[i]

            mass,err = spectrums[i].get_masses()
            mass_,err_ = [],[]
            for i in range(max_state):
                mass_.append(mass[i])
                err_.append(err[i])
            mass = mass_
            err = err_

            i = range(len(mass))
            xs = [volume+0.001*mas for mas in i]
            nams = [i for i in range(54)]
            for skip in skps:
                if skip in nams:
                    nams.remove(skip)
            colors = ['black' for l in range(max_state)]
            dcol = color_coding_file_simple()
            for j in range(max_state):
                name = f"{irr}/Volume_{volume}/t0{t0}/StateColorFiles/state{j}.txt"
                if os.path.exists(name):
                    with open(name) as f:
                        lines = f.readlines()
                        number = int(lines[0])
                
                    colors[j] = dcol[number]

            

            
            for i in range(len(mass)):
                if mass[i] <= th:
                    if Ps == "100":
                        m = mass[i]**2-((2*np.pi)/(xs[i]*3.444))**2
                        m = np.sqrt(m)
                    else:
                        m = mass[i]

                    ax.errorbar(xs[i],m,yerr=err[i],fmt='o',c=colors[i],markersize=3)
                    ax.text(xs[i],m,nams[i],fontsize = 8)
        xs = []
        errs = []
        ys = []
        files = os.listdir(path_to_folder)
        print(files)
        path_files = []

        for file in files:
            path_files.append(f"{path_to_folder}/{file}")
        for path_to_file in path_files:
            with open(path_to_file,'r') as f:
                lines = f.readlines()
                level_E = []
                level_L = []
                level_err = []
                for line in lines:
                    elems = line.split("   ")
                    volume = float(elems[0])
                    elm2 = elems[1].split()
                    if "v=" in elm2[1]:
                        energy = float(elm2[0].strip())
                        error = 0.0
                    else:
                        energy = float(elm2[0].strip())
                        error = float(elm2[1].strip())

                        
                    level_err.append(error)
                    level_E.append(energy)
                    level_L.append(volume)
                ax.plot(level_L,level_E,marker = '',linestyle = '--',color = 'gold',markersize = 3)
                ax.fill_between(level_L,np.array(level_E)-np.array(level_err),np.array(level_E)+np.array(level_err),color = 'gold',alpha = 0.3)
        """with open(path_to_file,'r') as f:
            lines = f.readlines()
            for line in lines:
                elems = line.split("   ")
                volume = float(elems[0])
                elm2 = elems[1].split()
                energy = float(elm2[0])
                error = float(elm2[1])
                
                errs.append(error)
                xs.append(volume)
                ys.append(energy)
       
        #plt.errorbar(xs,ys,yerr=errs,color = 'gold',markersize = 3,marker = '*',linestyle = '')
        initial_L = xs[0]
        level = []
        levels = []
        levels = {}
        count = 0
        for i,L in enumerate(xs):
            if L != initial_L:
                level = ys[i]
                err = errs[i]
                initial_L = L
                count = 0
                if count not in levels.keys():
                    levels[count] = [(level,L,err)]
                else:
                    levels[count] += [(level,L,err)]
                count += 1
            else:
                level = ys[i]
                err = errs[i]
                if count not in levels.keys():
                    levels[count] = [(level,L,err)]
                else:
                    levels[count] += [(level,L,err)]
                count += 1
        for key in levels.keys():
            xs = []
            ys = []
            errs = []
            for elem in levels[key]:
                ys.append(elem[0])
                xs.append(elem[1])
                errs.append(elem[2])
                
            ax.plot(xs,ys,marker = '',linestyle = '--',color = 'gold',markersize = 3)
            #ax.errorbar(xs,ys,yerr=errs,color = 'gold',markersize = 3,marker = '*',linestyle = '')
            ax.fill_between(xs,np.array(ys)-np.array(errs),np.array(ys)+np.array(errs),color = 'gold',alpha = 0.3)
            

"""
            


        #plt.fill_between(xs,np.array(ys)-np.array(errs),np.array(ys)+np.array(errs),color = 'gold',alpha = 0.3)

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        #plt.xlabel('Volume')
        #plt.ylabel('$a_tE_{cm}$',rotation = 0)
        irrep_name = irr 
        name = f"{irrep_name}/Spectrum_{irrep_name}.pdf"
        ax.legend(ncol = 2,loc = 'lower left',fontsize = 6)
        #plt.savefig(name)
        #plt.show()"