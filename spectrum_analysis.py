import numpy as np
import os
import matplotlib.pyplot as plt
from colors import color_coding_file
from colors import color_coding_file_simple
from colors import color_coding_dict
from colors import save_color_code_state
import sys
from colors import operator_identification
import grouptheory as gt
import particle as p
import no_int_E_levels as no_int
def read_ax_file(filename):
    ## Check if file empty

    blue_initial_1 = []
    blue_initial_2 = []
    blue_initial_3 = []
    with open(filename, 'r') as f:
        lines = f.readlines()
        line1 = lines[0]
        max_x = float(line1.split(' ')[2])
        min_x = float(line1.split(' ')[1])
        block = []
        for i in range(1,len(lines)):
            if '#e0' in lines[i] or '#c0' in lines[i] or '#m' in lines[i]:
                continue
            block += [lines[i]]
        #return block
    block1 = []
    dic = {0:[], 1:[], 2:[], 3:[],4:[],5:[],6:[],7:[],8:[],9:[],10:[],11:[]}
    index = 0
    flag = True

    for i in range(len(block)):
        if block[i] == '\n' and flag:
            flag = False
            index += 1
            continue
        elif block[i] == '\n' and not flag:
            continue
        elif block[i] == '#\n':
            flag = True
            continue
        elif "#y" in block[i]:
            min_y = float(block[i].split(' ')[1])
            max_y = float(block[i].split(' ')[2])
        elif block[i][0] == '#':
            continue

        dic[index] += [block[i]]

    return dic,max_x,min_x,max_y,min_y       
   

def list_to_x_y(list):
    x = []
    y = []
    for i in range(len(list)):
        x += [float(list[i].split('  ')[0])]
        y += [float(list[i].split('  ')[1])]
    return x, y
def list_to_x_y_err(list):
    x = []
    y = []
    err = []
    for i in range(len(list)):
        x += [float(list[i].split('  ')[0])]
        y += [float(list[i].split('  ')[1])]
        err += [float(list[i].split('  ')[2])]
    return x, y,err
def plot_ax_file(filename,show,save,save_path):
    dic,max_x,min_x,max_y,min_y = read_ax_file(filename)
    blue_1_fit_x = []
    blue_2_fit_x = []
    blue_3_fit_x = []
    blue_1_fit_y = []
    blue_2_fit_y = []
    blue_3_fit_y = []
    blue_1_fit_x, blue_1_fit_y = list_to_x_y(dic[0])
    blue_2_fit_x, blue_2_fit_y = list_to_x_y(dic[1])
    blue_3_fit_x, blue_3_fit_y = list_to_x_y(dic[2])
    red_1_fit_x, red_1_fit_y = list_to_x_y(dic[3])
    red_2_fit_x, red_2_fit_y = list_to_x_y(dic[4])
    red_3_fit_x, red_3_fit_y = list_to_x_y(dic[5])
    blue_4_fit_x, blue_4_fit_y = list_to_x_y(dic[6])
    blue_5_fit_x, blue_5_fit_y = list_to_x_y(dic[7])
    blue_6_fit_x, blue_6_fit_y = list_to_x_y(dic[8])
    not_in_color = 'teal'
    in_color = 'darkred'
    fig,ax = plt.subplots()
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.xlim(min_x,max_x)
    plt.ylim(min_y,max_y)

    plt.plot(blue_1_fit_x, blue_1_fit_y, label='blue_1',color = not_in_color,linewidth = 0.5)
    plt.plot(blue_2_fit_x, blue_2_fit_y, label='blue_2',color = not_in_color,linewidth = 0.5)
    plt.plot(blue_3_fit_x, blue_3_fit_y, label='blue_3',color = not_in_color,linewidth = 0.5)
    plt.plot(red_1_fit_x, red_1_fit_y, label='red_1',color = in_color,linewidth = 0.5)
    plt.plot(red_2_fit_x, red_2_fit_y, label='red_2',color = in_color,linewidth = 0.5)
    plt.plot(red_3_fit_x, red_3_fit_y, label='red_3',color = in_color,linewidth = 0.5)
    plt.plot(blue_4_fit_x, blue_4_fit_y, label='blue_4',color = not_in_color,linewidth = 0.5)
    plt.plot(blue_5_fit_x, blue_5_fit_y, label='blue_5',color = not_in_color,linewidth = 0.5)
    plt.plot(blue_6_fit_x, blue_6_fit_y, label='blue_6',color = not_in_color,linewidth = 0.5)

    data_in_x,data_in_y,data_in_err = list_to_x_y_err(dic[9])
    plt.errorbar(data_in_x,data_in_y,yerr=data_in_err,fmt='o',color=in_color,label='data',markersize=1.5)
    data_not_in_x,data_not_in_y,data_not_in_err = list_to_x_y_err(dic[10])
    plt.errorbar(data_not_in_x,data_not_in_y,yerr=data_not_in_err,fmt='o',color=not_in_color,label='data_not_in',markersize=1.5)
    #print(dic[11])
    title = dic[11][1]

    m = title.split('; ')[1]
    m_val = m.split('\+-')[0].split(')=')[1]
    m_err = m.split('\+-')[1].split('"\n')[0]
    mbefore = m.split(')')[0]
    m_title = f"${mbefore}) = {float(m_val)} \pm {float(m_err)}$"
    

    #print(m_title)
    before = title.split('; ')[0]
    #print(before)
    val = before.split('=')[1]
    chi_title = "$\\chi^2/N_{Dof} = $" + val
    title_def = chi_title + '\n' + m_title
    ax.set_xlabel('t')
    ax.set_ylabel('$\\lambda_n \cdot e^{m_n t}$')
    plt.title(title_def)
    if show:
        plt.show()
    if save:
        plt.savefig(save_path)
        plt.close()


    
def calc(file):
    data = np.loadtxt(file,skiprows=1)
    n = len(data)
    avg = 0
    err = 0
    for i in range(n):
        avg += data[i][1]/n
    for i in range(n):
        err += (data[i][1]-avg)**2
    err = np.sqrt(err/(n*(n-1)))
    return avg,err

class Spectrum(object):
    def __init__(self,volume,t0,irrep,create_files = False,saveAllPlots = False):
        self.irrep = irrep
        self.volume = volume
        self.t0 = t0
        self.states = []
        self.operators = []
        self.skips = []
        self.path_home = f"{irrep}\\Volume_{volume}\\t0{t0}"
        states = self.get_states()
        self.operators = self.get_operators()
        self.states = max(states)+1
        #self.states = max(states)
        for j in range(self.states):
            if j not in states:
                self.skips.append(j)
        if create_files:
            print("Creating files")
            self.obtain_Z_t0()
            self.obtain_mass()
            [self.renormalize(op) for op in range(self.operators)]
        if saveAllPlots:
            print("Saving all plots")
            for state in range(self.states):
                if state not in self.skips:
                    self.plot_correlator(state,False,True)
        

    def get_operators(self):
        ## Get the operators from the file ops.txt
        path = f"{self.irrep}\\Volume_{self.volume}\\ops.txt"
        with open(path, "r") as f:
            line = f.readlines()
            operators = len(line)
        return operators

    def get_states(self):
        ## Get the states from the files in MassJackFiles
        path = f"{self.path_home}\\MassJackFiles"
        files = os.listdir(path)
        states = []
        for file in files:
            split1 = file.split('state')[1]
            number = int(split1.split('.')[0])
            states.append(number)
        return states

    def create_names(self):
        names = []
        operators = [i for i in range(self.operators)]
        states = [i for i in range(self.states)]
        new_states = []
        for state in states:
            if state not in self.skips:
                new_states.append(state)
        for operator in operators:
            for state in new_states:
                names.append(f"Z_t0_{self.t0}_reorder_state{state}_op{operator}.jack")
        return names
    def create_names_operators(self,operator):
        names = []
        states = [i for i in range(self.states)]
        new_states = []
        for state in states:
            if state not in self.skips:
                new_states.append(state)
        for state in new_states:
                names.append(f"Z_t0_{self.t0}_reorder_state{state}_op{operator}.jack")
        return names
    def create_names_states(self,state,):
        names = []
        operators = [i for i in range(self.operators)]
        for operator in operators:
            names.append(f"Z_t0_{self.t0}_reorder_state{state}_op{operator}.jack")

        return names
    
    def obtain_Z_t0(self,):
        names = self.create_names()

        Z_t0 = []
        for name in names:
            path = f"{self.irrep}\\Volume_{self.volume}\\t0{self.t0}\\ZJackFiles\\{name}"
            path2 = f"{self.irrep}\\Volume_{self.volume}\\t0{self.t0}\\Zvalues\\{name}"
            with open(path2, "w") as f:
                f.write(f"{abs(calc(path)[0])} {calc(path)[1]}")
    def renormalize(self,operator):
        names = self.create_names_operators(operator)
        values = np.zeros(len(names))
        errors = np.zeros(len(names))
        for name in names:
            path = f"{self.irrep}\\Volume_{self.volume}\\t0{self.t0}\\Zvalues\\{name}"
            a = np.loadtxt(path)
            val,err = a[0],a[1]
            values[names.index(name)] = val
            errors[names.index(name)] = err
        maximum = np.max(values)
        for i in range(len(values)):
            values[i] = values[i]/maximum
            errors[i] = errors[i]/maximum
        
        for i,name in enumerate(names):
            path = f"{self.irrep}\\Volume_{self.volume}\\t0{self.t0}\\ZvaluesRenormalized\\{name}"
            with open(path, "w") as f:
                f.write(f"{values[i]} {errors[i]}")

    def plot_histogram_state(self,state):
        file = f"{self.irrep}\\Volume_{self.volume}\\ops.txt"
        names = self.create_names_states(state)
        values = np.zeros(len(names))
        errors = np.zeros(len(names))
        for name in names:
            path = f"{self.irrep}\\Volume_{self.volume}\\t0{self.t0}\\ZvaluesRenormalized\\{name}"
            a = np.loadtxt(path)
            val,err = a[0],a[1]
            values[names.index(name)] = val
            errors[names.index(name)] = err
        x = [f"op{i}" for i in range(self.operators)]
        dict_operators,color_code_dict = color_coding_file(file)
        colors = [color_code_dict[dict_operators[name]] for name in x]
        labels = []
        for op in x:

            if dict_operators[op] not in labels:
                labels.append(dict_operators[op])
            else:
                labels.append('_'+dict_operators[op])
        fig,ax = plt.subplots()
        ax.bar(x,values,yerr=errors,color=colors,label=labels)
        name2 =  f"mass_t0_{self.t0}_reorder_state{state}.jack"
        v = np.loadtxt(f"{self.irrep}\\Volume_{self.volume}\\t0{self.t0}\\MassValues\\{name2}")
        value,error = round(v[0],3),round(v[1],3)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        plt.title(f"State {state} $m= {value} \pm {error}$")

        plt.legend()
        plt.xticks([])
        plt.show()
    
    def plot_correlator(self,state,show,save = False):

        if state in self.skips or state >= self.states:
            raise ValueError("State not found")
        name = f"prin_corr_fit_t0{self.t0}_reorder_state{state}.ax"
        path =  f"{self.irrep}\\Volume_{self.volume}\\t0{self.t0}\\PrinCorrPlots\\{name}"
        name = f"prin_corr_fit_t0{self.t0}_reorder_state{state}.pdf"
        save_path = f"{self.irrep}\\Volume_{self.volume}\\t0{self.t0}\\CorrPlots\\{name}"
        ## Check if file empty
        if os.stat(path).st_size == 0:
            return 
        plot_ax_file(path,show,save,save_path)


        
    def create_names_mass(self,):
        names = []
        states = [i for i in range(self.states)]
        new_states = []
        for state in states:
            if state not in self.skips:
                new_states.append(state)
        for state in new_states:
                names.append(f"mass_t0_{self.t0}_reorder_state{state}.jack")
        return names

    def obtain_mass(self):
        names = self.create_names_mass()
        #print(names)
        mass = []
        for name in names:
            path = f"{self.irrep}\\Volume_{self.volume}\\t0{self.t0}\\MassJackFiles\\{name}"
            path2 = f"{self.irrep}\\Volume_{self.volume}\\t0{self.t0}\\MassValues\\{name}"
            with open(path2, "w") as f:
                f.write(f"{abs(calc(path)[0])} {calc(path)[1]}")

    def get_masses(self):
        names = self.create_names_mass()
        
        masses = []
        errors = []
        for name in names:
            path = f"{self.irrep}\\Volume_{self.volume}\\t0{self.t0}\\MassValues\\{name}"
            a = np.loadtxt(path)
            val,err = a[0],a[1]
            masses.append(val)
            errors.append(err)
        return masses,errors
    @staticmethod
    def plot_irrep_mass(spectrums,irrep_name,Ls,Ps,max_state,channel):
        colors = ['red','blue','green','orange','purple','pink','black','brown','teal','cyan','magenta','grey','lime','olive','yellow','navy','maroon','aqua','fuchsia','silver','red','blue','green','orange','purple','pink','black','brown','yellow','cyan','magenta','grey','lime','olive','teal','navy','maroon','aqua','fuchsia','silver']
        C_p = '+' if channel['C_parity'] == 1 else '-'
        threshold = 0.72
        channelsDs = p.Ds(channel,threshold)
        chh = p.channels(channel,threshold)[0]
        for key in channelsDs.keys():
            chh[key] = channelsDs[key]
        dict = {}
        channels = chh.keys()
        for i,chan in enumerate(channels):
            dict[chan] = colors[i]
        dict = color_coding_dict()
        save_sigmas = []

        

        


        dict2 = {}
        for i,chan in enumerate(channels):
            dict2[i] = chan
        ch = list(channels)
        counts = {}
        for i,chan in enumerate(channels):
            counts[chan] = 0
        labels = []
        fig, ax = plt.subplots()

        for P in Ps:

            Es,indeces,multiplcities = no_int.E_levels_in_irrep_P(irrep_name,channels,Ls,P,channel,4,threshold)
            # Assign a color to each channel
            #colors = plt.cm.get_cmap('hsv', len(ch))

            for i,E in enumerate(Es):
                counts[ch[indeces[i]]] += 1

                for j in range(multiplcities[i]):
                    
                    padding = np.array([0.001 for l in E])
                    label = '$'+ ch[indeces[i]][0] + ' ' +ch[indeces[i]][1] +'$'
                    if ch[indeces[i]][1] != '\sigma':
                        if label not in labels:
                            ax.plot(Ls,E+j*padding,color =dict[dict2[indeces[i]]],linewidth = 0.7,label = label)

                            labels.append(label)
                        else:
                            ax.plot(Ls,E+j*padding,color = dict[dict2[indeces[i]]],linewidth = 0.7)
                    else:
                        if ch[indeces[i]][0] not in save_sigmas:
                            save_sigmas.append(ch[indeces[i]][0])
        print(save_sigmas)
        #print(save_sigmas)

        plt.xlim(14,26)
        plt.ylim(0.6,0.74)
        plt.xticks([16,20,24])
        plt.yticks([0.62,0.64,0.66,0.68,0.70,0.72,0.74])

        for i,ch in enumerate(channels):
            if counts[ch] != 0:
                E = no_int.non_int_elevel_at_rest(ch[0],ch[1],np.array([0,0,0]),3.444,Ls)

                plt.plot(Ls,E,color = dict[dict2[i]],linewidth = 0.7,linestyle = '--')

        plt.legend()
        new_Ls = [16,20,24]
        for particle in save_sigmas:

            Es, multiplcities, Lss = plot_sigma_E_levels(particle, new_Ls,Ps,irrep_name)
            print(Es,particle,Lss)
            for i,E in enumerate(Es):
                for j in range(multiplcities[i]):
                    #padding = np.array([0.001 for l in Lss[i][0]])
                    chh = (particle,"\sigma")
                    #print(E,Lss[i])
                    #print(Lss[i][0])
                    plt.errorbar(Lss[i][0]+j*0.1,E[0],yerr=0.01,color = 'brown',fmt='o',markersize = 3)

        title = 'Irrep: $ '+irrep_name.split('^')[0] + "^{"+irrep_name.split('^')[1] +C_p+'}$'
        plt.title(title)
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
            if irr != irreps[i]:
                raise ValueError("Different irreps")

        for i,volume in enumerate(volumes):
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
                name = f"{irr}\\Volume_{volume}\\t0{t0}\\StateColorFiles\\state{j}.txt"
                if os.path.exists(name):
                    with open(name) as f:
                        lines = f.readlines()
                        number = int(lines[0])
                
                    colors[j] = dcol[number]

            

            
            for i in range(len(mass)):
                ax.errorbar(xs[i],mass[i],yerr=err[i],fmt='o',c=colors[i],markersize=3)
                ax.text(xs[i],mass[i],nams[i],fontsize = 8)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        plt.xlabel('Volume')
        plt.ylabel('$a_tE_{cm}$',rotation = 0)
        plt.show()

    @staticmethod
    def plot_spectrum_multiple(spectrums,max_state):
        
        

       

        volumes = []
        t0s = []
        skpss = []
        max_states = []
        irreps = []
        names = []
        for spectrum in spectrums:
            irreps.append(spectrum.irrep)
            volumes.append(spectrum.volume)
            t0s.append(spectrum.t0)
            skpss.append(spectrum.skips)
            max_states.append(max_state)
            names.append(spectrum.irrep + f" $Volume = {spectrum.volume}t_0 = {spectrum.t0}$")
        
        fig, ax = plt.subplots()
        plt.ylim(0.6,0.74)
        plt.yticks([0.62,0.64,0.66,0.68,0.70,0.72,0.74])


        for i,volume in enumerate(volumes):
            t0 = t0s[i]
            max_state = max_states[i]
            skps = skpss[i]
            irr = irreps[i]

            mass,err = spectrums[i].get_masses()
            mass_,err_ = [],[]
            for k in range(max_state):
                mass_.append(mass[k])
                err_.append(err[k])
            mass = mass_
            err = err_


            xs = [names[i] for j in range(len(mass))]
            nams = [i for i in range(54)]
            for skip in skps:
                if skip in nams:
                    nams.remove(skip)
            colors = ['black' for l in range(max_state)]
            dcol = color_coding_file_simple()
            for j in range(max_state):
                name = f"{irr}\\Volume_{volume}\\t0{t0}\\StateColorFiles\\state{j}.txt"
                if os.path.exists(name):
                    with open(name) as f:
                        lines = f.readlines()
                        number = int(lines[0])
                    colors[j] = dcol[number]
            

            for i in range(len(mass)):
                ax.errorbar(xs[i],mass[i],yerr=err[i],fmt='o',c=colors[i],markersize=3)
                ax.text(xs[i],mass[i],nams[i],fontsize = 8)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        plt.xlabel('Spectra')
        plt.ylabel('$a_tE_{cm}$',rotation = 0)
        plt.show()
    def automatic_coloring(self,max_state):
        num_to_color = color_coding_file_simple()
        color_to_num = {}
        for key in num_to_color.keys():
            color_to_num[num_to_color[key]] = key
        new_dict,color_code_dict_new = color_coding_file(f"{self.irrep}\\Volume_{self.volume}\\ops.txt")
        dict_color_ops = {}
        for op in new_dict.keys():
            dict_color_ops[op] = color_code_dict_new[new_dict[op]]
        

        for state in range(max_state):
            if state not in self.skips:
                Z_renormalized = []
                Z_renormalized_err = []
                for operator in range(self.operators):
                    path = f"{self.irrep}\\Volume_{self.volume}\\t0{self.t0}\\ZvaluesRenormalized\\Z_t0_{self.t0}_reorder_state{state}_op{operator}.jack"
                    a = np.loadtxt(path)
                    val,err = a[0],a[1]
                    Z_renormalized.append(val)
                values = {}
                total = 0
                for i,op in enumerate(new_dict.keys()):
                    if dict_color_ops[op] not in values.keys():
                        values[dict_color_ops[op]] = Z_renormalized[i]
                    else:
                        values[dict_color_ops[op]] += Z_renormalized[i]
                    total += Z_renormalized[i]
                possible_initial_sweep = []
                possible_hits = {}
                totals = 0
                for i,key in enumerate(new_dict.keys()):
                    if abs(Z_renormalized[i] -1) < 0.15:
                        if dict_color_ops[key] not in possible_initial_sweep:
                            possible_initial_sweep.append(dict_color_ops[key])
                            possible_hits[dict_color_ops[key]] = Z_renormalized[i]

                        else:
                            possible_hits[dict_color_ops[key]] += Z_renormalized[i]
                        totals += Z_renormalized[i]
                       
                            
                if len(possible_initial_sweep) == 1:
                    save_color_code_state(color_to_num[possible_initial_sweep[0]],state,self)
                    continue 
                elif len(possible_initial_sweep) > 1:
                    print(state)
                    value = None
                    keys = None
                    for key in possible_hits.keys():
                        
                        if possible_hits[key]/totals > 0.65:
                            value = possible_hits[key]
                            keys = key
                    if value != None:

                        save_color_code_state(color_to_num[keys],state,self)
                        continue
                    else:

                    
                        path = f"{self.irrep}\\Volume_{self.volume}\\t0{self.t0}\\Annotations\\Colorstate{state}.txt"
                        with open(path, "w") as f:
                            f.write("Not possible to automatically decide color, possible options are: \n")
                            f.write(str(possible_initial_sweep))
                    continue


                for key in values.keys():
                    values[key] = values[key]/total
                possible = []
                for key in values.keys():
                    if values[key] > 0.5:
                        possible.append((key,values[key]))
                if len(possible) == 1:
                    save_color_code_state(color_to_num[possible[0][0]],state,self)
                    continue
                elif len(possible) == 0:
                    possible = 0
                    value = 0
                    for key in values.keys():
                        if values[key] > value:
                            possible = key
                            value = values[key]
                    save_color_code_state(color_to_num[possible],state,self) 
                    continue

                   
                else:
                    path = f"{self.irrep}\\Volume_{self.volume}\\t0{self.t0}\\Annotations\\Colorstate{state}.txt"
                    with open(path, "w") as f:
                        f.write("Not possible to automatically decide color, possible options are: \n")
                        f.write(str(possible))
                        
                    
                
            
                
#obtain_mass(54,8,[38],24)
#plot_irrep_mass([8],np.linspace(14,26),[np.array([0,0,0]),np.array([1,0,0]),np.array([1,1,0]),np.array([1,1,1])],[24])

    
    
def plot_sigma_E_levels(particle1,Ls,Ps,irrep_name):
    channel = {'Charm':0,'Strange':0,'Isospin':0,'Charm_Isospin':1,'C_parity':-1}
    xi = 3.44
    energies = []
    multiplcities = []
    Ls = [16,20,24]

    Lss = []
    name_2 = "\sigma"
    ranking = []
    #print(particle1)
    for P in Ps:
        allowed,multiplcity = E_level_in_irrep_sigma_rest(irrep_name,particle1,name_2,P,channel,4,0.74)   
        
        if allowed:
            
            ## check if path exists
            sigma_E = []
            new_L = []
            for L in Ls:
                path = f"sigma/sigma_E_levels_{P[0]}{P[1]}{P[2]}_{L}.txt"
                if  os.path.exists(path):
                    for l in range(len(np.loadtxt(path))):
                        new_L += [L]
                        sigma_E += [np.loadtxt(path)[l]]
                        ranking += [l]
            #print(sigma_E)
            all_particles = p.read_particles('particles_unfl.txt')+p.read_particles('charmonium.txt')+p.read_particles('Ds.txt')
            for particle in all_particles:
                if particle.name == particle1:
                    particle_1 = particle
            mass_1 = float(particle_1.Mass)
            #particle_1_E = np.array([no_int.no_int_e_level_p(mass_1,P,3.444,L) for L in new_L])


            sigma_E = np.array(sigma_E)

            if len(sigma_E) == 0:
                continue
            for m,sigma in enumerate(sigma_E):
                #sigma[0] = sigma[0] - 2*mass_1
                #print(sigma)
                sigma_Es = sigma+ ((2*np.pi)/(new_L[m]*xi))**2*np.linalg.norm(P)**2
                particle_1_E = no_int.no_int_e_level_p(mass_1,P,3.444,new_L[m])
                print(particle_1_E,P)
                E_tot = particle_1_E + sigma_Es
                energies.append([E_tot])
                Lss.append([new_L[m]])
                multiplcities.append(multiplcity)

    return energies,multiplcities,Lss


def E_level_in_irrep_sigma_rest(irrep_name,name_1,name_2,P,channel,Jmax,threshold):
    Flag = False
    
    if 'D' in name_1 and 'D' in name_2:
        if name_1 == p.bar(name_2):
            Flag = True
            general_allowed_irreps_with_symmetry = p.possible_irreps_rest(name_1,name_2,channel,Jmax,threshold)
    group = no_int.identify_momentum_type(P)
    if group == 'O_D^h':
        irreps = gt.identify_ni_at_rest_levels_subduction(name_1,name_2)
    if group == 'Dic_4':
        irreps = gt.identify_ni_Dic4_levels_subduction(name_1,name_2)
    if group == 'Dic_2':
        irreps = gt.identify_ni_Dic2_levels_subduction(name_1,name_2)
    if group == 'Dic_3':
        irreps = gt.identify_ni_Dic3_levels_subduction(name_1,name_2)
    for irrep in irreps:
        name = ''
        count = 0
        num = ''
        for let in irrep:
            if let.isnumeric() and count == 0:
                num += let
                continue
            if let.isalpha():
                count += 1
            name += let
        if name == irrep_name:
            if Flag:
                if  name in general_allowed_irreps_with_symmetry:
                     return True,1 if num == '' else int(num)
                else:
                    return False,0
            else:
                return True,1 if num == '' else int(num)
    return False,0





#obtain_Z_t0(55,54,8,[38])     
#[renormalize(i,54,8,[38]) for i in range(55)]

