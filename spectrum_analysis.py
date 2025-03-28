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
def read_ax_file(filename):
    """
    Function used to read the ax files. These files are produced by semble_vfit and contain the information for the
    Principal correlator plots. The function will read the contents of the file and obtain the relevant pieces.
    
    :param filename: The name of the file to be read
    :return: A dictionary with the data for the plot and the max and min values for the x and y axis
    """
    

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
    """
    Function used to convert data of format "x y" where x,y are floats to the tuple (x,y)
    :param list: A list of strings of format "x y"
    :return: Two lists, one with the x values and one with the y values
    """
    x = []
    y = []
    for i in range(len(list)):
        x += [float(list[i].split('  ')[0])]
        y += [float(list[i].split('  ')[1])]
    return x, y
def list_to_x_y_err(list):
    """
    Function used to convert data of format "x y err" where x,y,err are floats to the tuple (x,y,err)
    :param list: A list of strings of format "x y err"
    :return: Three lists, one with the x values, one with the y values and one with the errors
    """
    x = []
    y = []
    err = []
    for i in range(len(list)):
        x += [float(list[i].split('  ')[0])]
        y += [float(list[i].split('  ')[1])]
        err += [float(list[i].split('  ')[2])]
    return x, y,err
def plot_ax_file(filename,show,save,save_path):
    """
    Function used to plot the data from the ax files. The function will read the data from the file and will combine the 
    data according to a repeatable pattern present in all the files
    (at least the ones corresponding to Principal correlator plots)
    """

    dic,max_x,min_x,max_y,min_y = read_ax_file(filename)

    # Get the points of the 3 fits (1-sigma region) in the region from 0 to the first fitted point
    # In the plot these is shown as a blue line
    blue_1_fit_x = []
    blue_2_fit_x = []
    blue_3_fit_x = []
    blue_1_fit_y = []
    blue_2_fit_y = []
    blue_3_fit_y = []
    blue_1_fit_x, blue_1_fit_y = list_to_x_y(dic[0])
    blue_2_fit_x, blue_2_fit_y = list_to_x_y(dic[1])
    blue_3_fit_x, blue_3_fit_y = list_to_x_y(dic[2])

    # Get the points of the 3 fits (1-sigma region) in the region from the first fitted point to the last fitted point
    # In the plot these is shown as a red line
    red_1_fit_x, red_1_fit_y = list_to_x_y(dic[3])
    red_2_fit_x, red_2_fit_y = list_to_x_y(dic[4])
    red_3_fit_x, red_3_fit_y = list_to_x_y(dic[5])

    # Get the points of the 3 fits (1-sigma region) in the region from the last fitted point to the end of the data
    # In the plot these is shown as a blue line
    blue_4_fit_x, blue_4_fit_y = list_to_x_y(dic[6])
    blue_5_fit_x, blue_5_fit_y = list_to_x_y(dic[7])
    blue_6_fit_x, blue_6_fit_y = list_to_x_y(dic[8])

    # Blue points are the data points that are not included in the fit
    not_in_color = 'teal'

    # Red points are the data points that are included in the fit
    in_color = 'darkred'
    fig,ax = plt.subplots()

    # Set up plot
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.xlim(min_x,max_x)
    plt.ylim(min_y,max_y)

    # Plot the fits

    plt.plot(blue_1_fit_x, blue_1_fit_y, label='blue_1',color = not_in_color,linewidth = 0.5)
    plt.plot(blue_2_fit_x, blue_2_fit_y, label='blue_2',color = not_in_color,linewidth = 0.5)
    plt.plot(blue_3_fit_x, blue_3_fit_y, label='blue_3',color = not_in_color,linewidth = 0.5)
    plt.plot(red_1_fit_x, red_1_fit_y, label='red_1',color = in_color,linewidth = 0.5)
    plt.plot(red_2_fit_x, red_2_fit_y, label='red_2',color = in_color,linewidth = 0.5)
    plt.plot(red_3_fit_x, red_3_fit_y, label='red_3',color = in_color,linewidth = 0.5)
    plt.plot(blue_4_fit_x, blue_4_fit_y, label='blue_4',color = not_in_color,linewidth = 0.5)
    plt.plot(blue_5_fit_x, blue_5_fit_y, label='blue_5',color = not_in_color,linewidth = 0.5)
    plt.plot(blue_6_fit_x, blue_6_fit_y, label='blue_6',color = not_in_color,linewidth = 0.5)

    # Obtain data points included in the fit and data points not included in the fit
    data_in_x,data_in_y,data_in_err = list_to_x_y_err(dic[9])
    data_not_in_x,data_not_in_y,data_not_in_err = list_to_x_y_err(dic[10])

    # Plot the data points
    plt.errorbar(data_in_x,data_in_y,yerr=data_in_err,fmt='o',color=in_color,label='data',markersize=1.5)
    plt.errorbar(data_not_in_x,data_not_in_y,yerr=data_not_in_err,fmt='o',color=not_in_color,label='data_not_in',markersize=1.5)
    
    # Final set up of plot
    title = dic[11][1]

    m = title.split('; ')[1]
    m_val = m.split('\+-')[0].split(')=')[1]
    m_err = m.split('\+-')[1].split('"\n')[0]
    mbefore = m.split(')')[0]
    m_title = f"${mbefore}) = {float(m_val)} \pm {float(m_err)}$"
    before = title.split('; ')[0]
    val = before.split('=')[1]
    chi_title = "$\\chi^2/N_{Dof} = $" + val
    title_def = chi_title + '\n' + m_title
    ax.set_xlabel('t')
    ax.set_ylabel('$\\lambda_n \cdot e^{m_n t}$')
    plt.title(title_def)

    # Show or save the plot
    if show:
        plt.show()
    if save:
        plt.savefig(save_path)
        plt.close()
def plot_ax_file_in_fig(self,filename,axis,state):
    """
    Function used to plot the data from the ax files. The function will read the data from the file and will combine the 
    data according to a repeatable pattern present in all the files
    (at least the ones corresponding to Principal correlator plots)
    """

    dic,max_x,min_x,max_y,min_y = read_ax_file(filename)

    # Get the points of the 3 fits (1-sigma region) in the region from 0 to the first fitted point
    # In the plot these is shown as a blue line
    blue_1_fit_x = []
    blue_2_fit_x = []
    blue_3_fit_x = []
    blue_1_fit_y = []
    blue_2_fit_y = []
    blue_3_fit_y = []
    blue_1_fit_x, blue_1_fit_y = list_to_x_y(dic[0])
    blue_2_fit_x, blue_2_fit_y = list_to_x_y(dic[1])
    blue_3_fit_x, blue_3_fit_y = list_to_x_y(dic[2])

    # Get the points of the 3 fits (1-sigma region) in the region from the first fitted point to the last fitted point
    # In the plot these is shown as a red line
    red_1_fit_x, red_1_fit_y = list_to_x_y(dic[3])
    red_2_fit_x, red_2_fit_y = list_to_x_y(dic[4])
    red_3_fit_x, red_3_fit_y = list_to_x_y(dic[5])

    # Get the points of the 3 fits (1-sigma region) in the region from the last fitted point to the end of the data
    # In the plot these is shown as a blue line
    blue_4_fit_x, blue_4_fit_y = list_to_x_y(dic[6])
    blue_5_fit_x, blue_5_fit_y = list_to_x_y(dic[7])
    blue_6_fit_x, blue_6_fit_y = list_to_x_y(dic[8])

    # Blue points are the data points that are not included in the fit
    not_in_color = 'grey'

    # Red points are the data points that are included in the fit
    in_color = 'darkred'

    # Set up plot
    #axis.spines['right'].set_visible(False)
    #axis.spines['top'].set_visible(False)
    #axis.set_xlim(min_x,max_x)
    #axis.set_ylim(min_y,max_y)

    # Plot the fits

    axis.plot(blue_1_fit_x, blue_1_fit_y, label='blue_1',color = not_in_color,linewidth = 0.5)
    axis.plot(blue_2_fit_x, blue_2_fit_y, label='blue_2',color = not_in_color,linewidth = 0.5)
    axis.plot(blue_3_fit_x, blue_3_fit_y, label='blue_3',color = not_in_color,linewidth = 0.5)
    axis.plot(red_1_fit_x, red_1_fit_y, label='red_1',color = in_color,linewidth = 0.5)
    axis.plot(red_2_fit_x, red_2_fit_y, label='red_2',color = in_color,linewidth = 0.5)
    axis.plot(red_3_fit_x, red_3_fit_y, label='red_3',color = in_color,linewidth = 0.5)
    axis.plot(blue_4_fit_x, blue_4_fit_y, label='blue_4',color = not_in_color,linewidth = 0.5)
    axis.plot(blue_5_fit_x, blue_5_fit_y, label='blue_5',color = not_in_color,linewidth = 0.5)
    axis.plot(blue_6_fit_x, blue_6_fit_y, label='blue_6',color = not_in_color,linewidth = 0.5)

    # Obtain data points included in the fit and data points not included in the fit
    data_in_x,data_in_y,data_in_err = list_to_x_y_err(dic[9])
    data_not_in_x,data_not_in_y,data_not_in_err = list_to_x_y_err(dic[10])

    # Plot the data points
    axis.errorbar(data_in_x,data_in_y,yerr=data_in_err,fmt='o',color=in_color,label='data',markersize=1.5)
    axis.errorbar(data_not_in_x,data_not_in_y,yerr=data_not_in_err,fmt='o',color=not_in_color,label='data_not_in',markersize=1.5)
    
    # Final set up of plot
    title = dic[11][1]

    m = title.split('; ')[1]
    m_val = m.split('\+-')[0].split(')=')[1]
    m_err = m.split('\+-')[1].split('"\n')[0]
    mbefore = m.split(')')[0]
    
    m_title = f"$a_t m_{state} = {float(m_val)} \pm {float(m_err)}$"
    before = title.split('; ')[0]
    val = before.split('=')[1]
    chi_title = "$\\chi^2/N_{Dof} = $" + val
    title_def = f"n = {state}, "+chi_title + '\n' + m_title
    #axis.set_xlabel('t')
    #axis.set_ylabel('$\\lambda_n \cdot e^{m_n t}$')
    hist_axis = axis.inset_axes([0.2,0.55,0.3,0.4])
    self.plot_histogram_state_axis(state,hist_axis)
    axis.text(0.6,0.7,title_def,transform=axis.transAxes,fontsize=8)
    #plt.title(title_def)

    # Show or save the plot
    

    
def calc(file):
    """
    Function used to calculate the average and error of the data in a file
    :param file: The name of the file to be read
    :return: The average and error of the data in the file
    """
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
    """
    Class used to represent the spectrum of a given irrep, volume and t0. The class will contain the information of the
    states, operators and the mass of the states. The class will also contain methods to plot the correlators and the
    histogram of the states.
    """
    def __init__(self,volume,t0,irrep,create_files = False,saveAllPlots = False,saveHistPlots = False, pathData = "",pathOps = ""):
        """
        Constructor for the Spectrum class. The constructor will create the necessary files if create_files is True.
        The constructor will also save all the plots if saveAllPlots is True and will save the histogram plots if
        saveHistPlots is True. For this to work it needs to have first some files from semble_vfit in the directory
        path ={irrep}/Volume_{volume}/t0{t0}. It needs the list of operators (ops.txt) and the folders
        ZJackFiles and MassJackFiles and PrinCorrPlots produced by semble_vfit. Once the necessary files are created after
        a create_files = True, the class will be able to plot the correlators and the histogram plots and successive runs
        can be done with create_files = False, same for the saveAllPlots and saveHistPlots.
        :param volume: The volume of the data
        :param t0: The t0 of the data
        :param irrep: The irrep of the data
        :param create_files: A boolean that determines if the files should be created
        :param saveAllPlots: A boolean that determines if all the plots should be saved
        :param saveHistPlots: A boolean that determines if the histogram plots should be saved
        """
        self.irrep = irrep
        self.volume = volume
        self.t0 = t0
        self.states = []
        self.operators = []
        self.skips = []
        self.create_files = create_files
        if pathOps == "":
            self.pathOps = f"{irrep}/Volume_{volume}/ops.txt"
        else:
            self.pathOps = pathOps
            path_aux = f"{irrep}/Volume_{volume}/ops.txt"
            with open(self.pathOps, "r") as f:
                with open(path_aux, "w") as f2:
                    lines = f.readlines()
                    for operator in lines:
                        f2.write(operator)
        if pathData == "":
            self.pathData = f"{irrep}/Volume_{volume}/t0{t0}"
        else:
            self.pathData = pathData
        self.path_home = f"{irrep}/Volume_{volume}/t0{t0}"
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
        if saveHistPlots:
            print("Saving histogram plots")
            for state in range(self.states):
                if state not in self.skips:
                    self.plot_histogram_state(state,saveHistPlots)
        

    def get_operators(self):
        """
        Function used to get the number of operators used to extract the spectrum. It reads the ops.txt file in the
        directory {irrep}/Volume_{volume}/ops.txt
        :return: The number of operators used to extract the spectrum
        """
        path = self.pathOps
        with open(path, "r") as f:
            line = f.readlines()
            operators = len(line)
        return operators

    def get_states(self):
        """
        Function used to get the states (numbered) from the files in MassJackFiles
        :return: A list of the states present in the files in MassJackFiles
        """

        ## Get the states from the files in MassJackFiles
        if self.create_files:
            path = f"{self.pathData}/MassJackFiles"
        else:
            path =  f"{self.pathData}/Massvalues"
        files = os.listdir(path)
        states = []
        for file in files:
            split1 = file.split('state')[1]
            number = int(split1.split('.')[0])
            states.append(number)
        return states

    def create_names(self):
        """
        Function used to create the names of the files in ZJackFiles
        :return: A list of the names of the files in ZJackFiles
        """
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
        """
        Function used to create the names of the files in ZJackFiles for a given operator across all states.
        Meaning a list of names correspoding to the overlaps of a given operators across all states.
        :param operator: The operator for which the names of the files are to be returned
        :return: A list of the names of the files in ZJackFiles for the given operator
        """
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
        """
        Function used to create the names of the files in ZJackFiles for a given state across all operators.
        Meaning a list of names correspoding to the overlaps of a given state across all operators.
        :param state: The state for which the names of the files are to be saved
        :return: A list of the names of the files in ZJackFiles for the given state
        """
        names = []
        operators = [i for i in range(self.operators)]
        for operator in operators:
            names.append(f"Z_t0_{self.t0}_reorder_state{state}_op{operator}.jack")

        return names
    
    def obtain_Z_t0(self,):
        """
        Function used to calculate the average and error of the data in the files in ZJackFiles. The function will
        calculate the average and error of the data in the files in ZJackFiles and will save the results in the
        directory {irrep}/Volume_{volume}/t0{t0}/Zvalues.
        """

        names = self.create_names()

        for name in names:

            path = f"{self.pathData}/ZJackFiles/{name}"
            path2 = f"{self.irrep}/Volume_{self.volume}/t0{self.t0}/Zvalues/{name}"
            with open(path2, "w") as f:
                f.write(f"{abs(calc(path)[0])} {calc(path)[1]}")

    def renormalize(self,operator):
        """
        Function used to renormalize the data in the files in ZJackFiles. The function will renormalize the data in the
        files in ZJackFiles and will save the results in the directory {irrep}/Volume_{volume}/t0{t0}/ZvaluesRenormalized.
        The normalization is performed such that the overlap for a given state and operator is divided by the maximum of the set
        formed by the overlaps for such operator across all states.
        :param operator: The operator for which the data is to be renormalized
        """
        names = self.create_names_operators(operator)
        values = np.zeros(len(names))
        errors = np.zeros(len(names))
        for name in names:
            path = f"{self.irrep}/Volume_{self.volume}/t0{self.t0}/Zvalues/{name}"
            a = np.loadtxt(path)
            val,err = a[0],a[1]
            values[names.index(name)] = val
            errors[names.index(name)] = err
        maximum = np.max(values)
        for i in range(len(values)):
            values[i] = values[i]/maximum
            errors[i] = errors[i]/maximum
        
        for i,name in enumerate(names):
            path = f"{self.irrep}/Volume_{self.volume}/t0{self.t0}/ZvaluesRenormalized/{name}"
            with open(path, "w") as f:
                f.write(f"{values[i]} {errors[i]}")

    def plot_histogram_state(self,state,save = False):
        """
        Function used to plot the histogram of the overlaps for a given state. The function will read the data from the
        files in ZJackFiles and will plot the histogram of the overlaps of the given state with all operators.
        It will assign colors depending on the operators, if two operators correspond to the same pair of particles, they will
        be assigned the same color. The function will save the plot if save is True.
        :param state: The state for which the histogram is to be plotted
        :param save: A boolean that determines if the histogram plot should be saved
        """
        file = self.pathOps
        names = self.create_names_states(state)
        values = np.zeros(len(names))
        errors = np.zeros(len(names))
        for name in names:
            path = f"{self.irrep}/Volume_{self.volume}/t0{self.t0}/ZvaluesRenormalized/{name}"
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
        v = np.loadtxt(f"{self.irrep}/Volume_{self.volume}/t0{self.t0}/MassValues/{name2}")
        value,error = round(v[0],3),round(v[1],3)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        plt.title(f"State {state} $m= {value} \pm {error}$")

        plt.legend()
        plt.xticks([])
        if not save:
            plt.show()
        else:
            name = f"histogram_state{state}.pdf"
            path = f"{self.irrep}/Volume_{self.volume}/t0{self.t0}/HistogramPlots/{name}"
            plt.savefig(path)
            plt.close()
    def plot_histogram_state_axis(self,state,axis):
        """
        Function used to plot the histogram of the overlaps for a given state. The function will read the data from the
        files in ZJackFiles and will plot the histogram of the overlaps of the given state with all operators.
        It will assign colors depending on the operators, if two operators correspond to the same pair of particles, they will
        be assigned the same color. The function will save the plot if save is True.
        :param state: The state for which the histogram is to be plotted
        :param save: A boolean that determines if the histogram plot should be saved
        """
        file = self.pathOps
        names = self.create_names_states(state)
        values = np.zeros(len(names))
        errors = np.zeros(len(names))
        for name in names:
            path = f"{self.irrep}/Volume_{self.volume}/t0{self.t0}/ZvaluesRenormalized/{name}"
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
        x = np.arange(len(x))
        axis.barh(x,values,yerr=errors,color=colors,label=labels)
        name2 =  f"mass_t0_{self.t0}_reorder_state{state}.jack"
        v = np.loadtxt(f"{self.irrep}/Volume_{self.volume}/t0{self.t0}/MassValues/{name2}")
        value,error = round(v[0],3),round(v[1],3)
        #axis.spines['right'].set_visible(False)
        #axis.spines['top'].set_visible(False)

        #plt.title(f"State {state} $m= {value} \pm {error}$")

        #plt.legend()
        axis.set_xticks([])
        axis.set_yticks([])
        
    
    def plot_correlator(self,state,show,save = False):
        """
        Function used to plot the correlator for a given state. The function will read the data from the files in
        PrinCorrPlots and will plot the correlator for the given state. The function will save the plot if save is True.
        :param state: The state for which the correlator is to be plotted
        :param show: A boolean that determines if the plot should be shown
        :param save: A boolean that determines if the plot should be saved
        """

        if state in self.skips or state >= self.states:
            raise ValueError("State not found")
        name = f"prin_corr_fit_t0{self.t0}_reorder_state{state}.ax"
        path =  f"{self.pathData}/PrinCorrPlots/{name}"
        name = f"prin_corr_fit_t0{self.t0}_reorder_state{state}.pdf"
        save_path = f"{self.irrep}/Volume_{self.volume}/t0{self.t0}/CorrPlots/{name}"
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
            path = f"{self.pathData}/MassJackFiles/{name}"
            path2 = f"{self.irrep}/Volume_{self.volume}/t0{self.t0}/MassValues/{name}"
            with open(path2, "w") as f:
                f.write(f"{abs(calc(path)[0])} {calc(path)[1]}")

    def get_masses(self):
        names = self.create_names_mass()
        
        masses = []
        errors = []
        for name in names:
            path = f"{self.irrep}/Volume_{self.volume}/t0{self.t0}/MassValues/{name}"
            a = np.loadtxt(path)
            val,err = a[0],a[1]
            masses.append(val)
            errors.append(err)
        return masses,errors
    @staticmethod
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

                    ax.errorbar(xs[i],mass[i],yerr=err[i],fmt='o',c=colors[i],markersize=3)
                    ax.text(xs[i],mass[i],nams[i],fontsize = 8)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        plt.xlabel('Volume')
        plt.ylabel('$a_tE_{cm}$',rotation = 0)
        irrep_name = irr 
        name = f"{irrep_name}/Spectrum_{irrep_name}.pdf"
        plt.legend(fontsize = 6,ncol = 2)
        plt.savefig(name)
        plt.show()
    def plot_irrep_mass_no_sigma(spectrums,irrep_name,Ls,Ps,max_state,channel,th):
    
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
        for spectrum in spectrums:
            irreps.append(spectrum.irrep)
            volumes.append(spectrum.volume)
            t0s.append(spectrum.t0)
            skpss.append(spectrum.skips)
            max_states.append(max_state)
        irr = irreps[0]
        for i in range(1,len(irreps)):
            pass
            """if irr != irreps[i]:
                raise ValueError("Different irreps")"""

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

                    ax.errorbar(xs[i],mass[i],yerr=err[i],fmt='o',c=colors[i],markersize=3)
                    ax.text(xs[i],mass[i],nams[i],fontsize = 8)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        plt.xlabel('Volume')
        plt.ylabel('$a_tE_{cm}$',rotation = 0)
        irrep_name = irr 
        name = f"{irrep_name}/Spectrum_{irrep_name}.pdf"
        plt.legend()
        plt.savefig(name)
        plt.show()
    def plot_irrep_mass_no_sigma_grey_out(spectrums,irrep_name,Ls,Ps,max_state,channel,th,included):
    
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
        plt.xlabel('Volume')
        plt.ylabel('$a_tE_{cm}$',rotation = 0)
        irrep_name = irr 
        name = f"{irrep_name}/Spectrum_{irrep_name}.pdf"
        plt.legend()
        #plt.savefig(name)
        plt.show()
    def plot_nine_levels(self,included,irrep,ymin,ymax,xmin,xmax):
        fig,axes = plt.subplots(4,3,figsize=(15,15),sharex=True,sharey=True)
        plt.subplots_adjust(wspace=0.1, hspace=0.1) 
        plt.ylim(ymin,ymax)
        plt.xlim(xmin,xmax)
        x_and_y_axis = axes[0,0]
        title = axes[0,1]
        legend = axes[0,2]
        x_and_y_axis.arrow(15,1,5,0,head_width=0.05, head_length=0.2, fc='k', ec='k')

        x_and_y_axis.arrow(15,1,0,0.2,head_width=0.2, head_length=0.04, fc='k', ec='k')
        #x_and_y_axis.set_xlim(-1,1.5)
        #x_and_y_axis.set_ylim(-1,1.5)
        x_and_y_axis.text(20,0.9,'$t$',fontsize=12)
        x_and_y_axis.text(9,1.2,'$\\lambda_n \cdot e^{m_n t}$',fontsize=12)
        x_and_y_axis.axis('off')

        #title.set_xlim(0,1)
        #title.set_ylim(0,1)
        title.text(10,1.1,f'{irrep} $L/a_s$ = {self.volume}',fontsize=12)
        
        title.axis('off')
        file =self.pathOps
      
     
        x = [f"op{i}" for i in range(self.operators)]
        dict_operators,color_code_dict = color_coding_file(file)
        colors = [color_code_dict[dict_operators[name]] for name in x]
        labels = []
        for op in x:

            if dict_operators[op] not in labels:
                labels.append(dict_operators[op])
            else:
                labels.append('_'+dict_operators[op])
        legend_patches = [mpatches.Patch(color=colors[::-1][i],label=labels[::-1][i]) for i in range(len(labels))]
        legend.legend(handles=legend_patches,loc='center',fontsize=5,  handleheight=1, handlelength=10,frameon=False,ncol = 2)
        legend.axis('off')

        #x_and_y_axis.aspect = 'equal'


        for i in range(9):
            state = included[i]
            name = f"prin_corr_fit_t0{self.t0}_reorder_state{state}.ax"
            path =  f"{self.irrep}/Volume_{self.volume}/t0{self.t0}/PrinCorrPlots/{name}"
            ax = axes[(i+3)//3,i%3]
            plot_ax_file_in_fig(self,path,ax,state)

            #ax.set_title(f"State {state}")
        plt.show()
            

    def plot_irrep_mass_with_get_finite(spectrums,irrep_name,Ls,Ps,max_state,channel,th,path_to_file):
    
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
        """
        
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
        """
            
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
        plt.xlabel('Volume')
        plt.ylabel('$a_tE_{cm}$',rotation = 0)
        irrep_name = irr 
        name = f"{irrep_name}/Spectrum_{irrep_name}.pdf"
        plt.legend()
        #plt.savefig(name)
        plt.show()
    def plot_irrep_mass_superimposed(spectrums,irrep_name,Ls,Ps,max_state,channel,th):
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

        #title = 'Irrep: $ '+irrep_name.split('^')[0] + "^{"+irrep_name.split('^')[1] +C_p+'}$'
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
        #irr = irreps[0]
        for i in range(1,len(irreps)):
            pass
            """if irr != irreps[i]:
                raise ValueError("Different irreps")"""

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

                    ax.errorbar(xs[i],mass[i],yerr=err[i],fmt='o',c=colors[i],markersize=3)
                    ax.text(xs[i],mass[i],nams[i],fontsize = 8)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        plt.xlabel('Volume')
        plt.ylabel('$a_tE_{cm}$',rotation = 0)
        #irrep_name = irr 
        #name = f"{irrep_name}/Spectrum_{irrep_name}.pdf"
        plt.show()    
    
    def obtain_levels_of_color(self,color):
        irr = self.irrep
        volume = self.volume
        t0 = self.t0
        levels = []
        dcol = color_coding_file_simple()
        for i in range(self.states):
            if i not in self.skips:
                name = f"{irr}/Volume_{volume}/t0{t0}/StateColorFiles/state{i}.txt"
                if os.path.exists(name):
                    with open(name) as f:
                        lines = f.readlines()
                        number = int(lines[0])
                        color_state = dcol[number]
                    if color_state == color:
                        levels.append(i)
        return levels
               
    def plot_histogram_of_levels_of_color(self,color,level):
        levels = self.obtain_levels_of_color(color)
        for i,state in enumerate(levels):
            if i != level:
                continue
            file = self.pathOps
            names = self.create_names_states(state)
            values = np.zeros(len(names))
            errors = np.zeros(len(names))
            for name in names:
                path = f"{self.irrep}/Volume_{self.volume}/t0{self.t0}/ZvaluesRenormalized/{name}"
                a = np.loadtxt(path)
                val,err = a[0],a[1]
                values[names.index(name)] = val
                errors[names.index(name)] = err
            x = [f"op{i}" for i in range(self.operators)]
            dict_operators,color_code_dict = color_coding_file(file)
            print(dict_operators)
            
            labels = []
            x_trues = []
            errs_trues = []
            val_trues = []
            for i,ops in enumerate(x):
                if color_code_dict[dict_operators[ops]] == color:
                    x_trues.append(ops)
                    val_trues.append(values[i])
                    errs_trues.append(errors[i])

            
            
            colors = [color_code_dict[dict_operators[name]] for name in x_trues]
            dict_operators2 = operator_identification_plus(file_name=file)
            xnew = [dict_operators2[op] for op in x_trues]

            fig,ax = plt.subplots()
            print(x_trues,xnew,val_trues,errs_trues)

            rgb_color = matplotlib.colors.to_rgb(color)
            ## For each operator in the plot make a new color similar to the original one
            new_colors = []
            alphas = np.linspace(0.1,1,len(x_trues))
            for i in range(len(x_trues)):

                new_colors.append((rgb_color[0],rgb_color[1],rgb_color[2],alphas[i]))

            ax.bar(xnew,val_trues,yerr=errs_trues,color=new_colors,label=xnew)
            name2 =  f"mass_t0_{self.t0}_reorder_state{state}.jack"
            v = np.loadtxt(f"{self.irrep}/Volume_{self.volume}/t0{self.t0}/MassValues/{name2}")
            value,error = round(v[0],3),round(v[1],3)
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)

            plt.title(f"State {state} $m= {value} \pm {error}$")

            plt.legend()
            ## Change the size of the text ticks
            plt.xticks(fontsize = 6)
            plt.xticks([])

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
                name = f"{irr}/Volume_{volume}/t0{t0}/StateColorFiles/state{j}.txt"
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
        new_dict,color_code_dict_new = color_coding_file(self.pathOps)
        dict_color_ops = {}
        for op in new_dict.keys():
            dict_color_ops[op] = color_code_dict_new[new_dict[op]]
        

        for state in range(max_state):
            
            if state not in self.skips:
                Z_renormalized = []
                Z_renormalized_err = []
                for operator in range(self.operators):
                    path = f"{self.irrep}/Volume_{self.volume}/t0{self.t0}/ZvaluesRenormalized/Z_t0_{self.t0}_reorder_state{state}_op{operator}.jack"
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
                    #print(state)
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

                    
                        path = f"{self.irrep}/Volume_{self.volume}/t0{self.t0}/Annotations/Colorstate{state}.txt"
                        with open(path, "w") as f:
                            #print("Not possible to automatically decide color, possible options are: \n")
                            #print(str(possible_initial_sweep),state)
                            f.write("Not possible to automatically decide color, possible options are: \n")
                            f.write(str(possible_initial_sweep))
                            save_color_code_state(0,state,self)
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
                else:
                    path = f"{self.irrep}/Volume_{self.volume}/t0{self.t0}/Annotations/Colorstate{state}.txt"
                    with open(path, "w") as f:
                        f.write("Not possible to automatically decide color, possible options are: \n")
                        f.write(str(possible))
                        #print("Not possible to automatically decide color, possible options are: \n")
                        #print(str(possible),state)
                        save_color_code_state(0,state,self)


                        
                """elif len(possible) == 0:
                    possible = 0
                    value = 0
                    for key in values.keys():
                        if state== 13:
                                print(key,values[key])
                        if values[key] > value:
                           
                            possible = key

                            value = values[key]
                            
                    
                    save_color_code_state(color_to_num[possible],state,self) 
                    continue
                """
                

                   
                
                    
                
            
                
#obtain_mass(54,8,[38],24)
#plot_irrep_mass([8],np.linspace(14,26),[np.array([0,0,0]),np.array([1,0,0]),np.array([1,1,0]),np.array([1,1,1])],[24])

    
    
def plot_sigma_E_levels(particle1,Ls,Ps,irrep_name):
    channel = {'Charm':0,'Strange':0,'Isospin':0,'Charm_Isospin':1,'C_parity':-1}
    xi = 3.44
    energies = []
    multiplcities = []
    errors = []
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
            err_E = []
            for L in Ls:
                path = f"sigma/sigma_E_levels_{P[0]}{P[1]}{P[2]}_{L}.txt"
                if  os.path.exists(path):
                    with open(path) as f:
                        for l,line in enumerate(f):
                            sigma_E += [float(line.split('+/-')[0])]
                            if len(line.split('+/-')) == 2:
                                err_E += [float(line.split('+/-')[1])]
                            else:
                                err_E += [0]
                            new_L += [L]
                            ranking += [l]
                            
                  
            #print(sigma_E)
            all_particles = p.read_particles('Particles/particles_unfl.txt')+p.read_particles('Particles/charmonium.txt')+p.read_particles('Particles/Ds.txt')
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
                err_e = err_E[m]
                errors.append(err_e)

                energies.append([E_tot])
                Lss.append([new_L[m]])
                multiplcities.append(multiplcity)
    return energies,multiplcities,Lss,errors


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
def plot_irrep_mass_in_flight(irrep_name,Ls,Ps,channel,th):
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
                else:
                    ax.plot(Ls,Es+j*padding,color = dict[chan])
        
        E_levels_sigma = no_int.get_E_levels_irrep_in_flight_sigma(c,irrep_name,Ps)
        print(E_levels_sigma)
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
                    else:
                        padding = np.array([0.1 for l in Ess])
                        ax.errorbar(Vss+padding*j,Ess,yerr=Erss,color = color,alpha = 0.3,fmt='o',markersize = 3)

            

        
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
        
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        plt.xlabel('Volume')
        plt.legend()
        plt.ylabel('$a_tE_{cm}$',rotation = 0)
        
        
        path = f"ni irreps/NI Spectrum_{irrep_name}_[{Ps}].pdf"
        plt.savefig(path)
        plt.show()
def bias_analysis_create_structure(irrep_name,V,tmaxs,t0min,t0max):
    t0s = range(t0min,t0max+1)
    for tmax in tmaxs:
        irrep = irrep_name + f"_tmax_{tmax}"

        create_structure(V,t0s,irrep)
def bias_analysis_load_and_create_plots(irrep_name,V,tmaxs,t0min,t0max):
    bias_analysis_create_structure(irrep_name,V,tmaxs,t0min,t0max)
    for tmax in tmaxs:
        irrep = irrep_name + f"_tmax_{tmax}"
        for t0 in range(t0min,t0max+1):


            
            path = f"{irrep}/Volume_{V}/t0{t0}/MassValues"
            spectrum = Spectrum(V,t0,irrep,create_files=False,saveAllPlots=False,saveHistPlots=False)
            mass,err = spectrum.get_masses()
            len_mass = len(mass)
            print(f"Automatic coloring for {irrep} t0 {t0} with {len_mass} states")
            spectrum.automatic_coloring(len_mass)
    
def create_structure(volume,t0s,irrep):
    if not os.path.exists(irrep):
        os.mkdir(irrep)
    if not os.path.exists(f"{irrep}/Volume_{volume}"):
        os.mkdir(f"{irrep}/Volume_{volume}")
    for t0 in t0s:
        path2 = f"{irrep}/Volume_{volume}/t0{t0}"
        if not os.path.exists(path2):
            os.mkdir(path2)
        path3 = f"{irrep}/Volume_{volume}/t0{t0}/StateColorFiles"
        if not os.path.exists(path3):
            os.mkdir(path3)
        path4 = f"{irrep}/Volume_{volume}/t0{t0}/ZValues"
        if not os.path.exists(path4):
            os.mkdir(path4)
        path5 = f"{irrep}/Volume_{volume}/t0{t0}/ZvaluesRenormalized"
        if not os.path.exists(path5):
            os.mkdir(path5)
        path6 = f"{irrep}/Volume_{volume}/t0{t0}/MassValues"
        if not os.path.exists(path6):
            os.mkdir(path6) 
        path8 = f"{irrep}/Volume_{volume}/t0{t0}/CorrPlots"
        if not os.path.exists(path8):
            os.mkdir(path8)
        path9 = f"{irrep}/Volume_{volume}/t0{t0}/Annotations"
        if not os.path.exists(path9):
            os.mkdir(path9)
        path10 = f"{irrep}/Volume_{volume}/t0{t0}/HistogramPlots"
        if not os.path.exists(path10):
            os.mkdir(path10)
#[create_structure(i,[8,9,10,11,12],"T1mM-fewer-djw-tmax30") for i in [24,16,20]]
def bias_analysis_level_matching(irrep_name,V,tmaxs,t0min,t0max,ncolor,tmax_ref,t0_ref,max_state):
    dict_spectrums = {}
    for tmaxs in tmaxs:
        irrep = irrep_name + f"_tmax_{tmaxs}"
        for t0 in range(t0min,t0max+1):
            individual_tag = irrep + f"_t0_{t0}"
            dict_spectrums[individual_tag] = Spectrum(V,t0,irrep,create_files=False,saveAllPlots=False,saveHistPlots=False)
    ref_tag = irrep_name + f"_tmax_{tmax_ref}" + f"_t0_{t0_ref}"
    ref_spectrum = dict_spectrums[ref_tag]
    names_ref = [i for i in range(max_state)]
    skips_ref = ref_spectrum.skips
    tags = list(dict_spectrums.keys())
    tags.remove(ref_tag)
    for skips in skips_ref:
        if skips in names_ref:
            names_ref.remove(skips)
    masses_ref,err_ref= ref_spectrum.get_masses()
    dict_closest_levels = {}
    for state in range(max_state):
        mass = masses_ref[state]
        name = names_ref[state]
        path = f"{irrep}/Volume_{V}/t0{t0_ref}/StateColorFiles/state{name}.txt"
        if os.path.exists(path):
                    with open(path) as f:
                        lines = f.readlines()
                        color_ref= int(lines[0])
        for tag in tags:
            spectrum = dict_spectrums[tag]
            masses,errs = np.array(spectrum.get_masses())
            names = [i for i in range(spectrum.states)]
            skips = spectrum.skips
            for skip in skips:
                if skip in names:
                    names.remove(skip)
            index = np.abs(masses - mass).argmin()  # Find the index of the closest value
            
            path_closest_state = f"{irrep}/Volume_{V}/t0{spectrum.t0}/StateColorFiles/state{names[index]}.txt"
            if os.path.exists(path_closest_state):
                    with open(path_closest_state) as f:
                        lines = f.readlines()
                        color = int(lines[0])
            
            if color != color_ref:
                differences = np.abs(masses - mass)
                tolerable = []
                indices = []
                for j,difference in enumerate(differences):
                    path_closest_state_new = f"{irrep}/Volume_{V}/t0{spectrum.t0}/StateColorFiles/state{names[j]}.txt"
                    if os.path.exists(path_closest_state_new):
                        with open(path_closest_state_new) as f:
                            lines = f.readlines()
                            color_new = int(lines[0])

                    if difference < 0.1 and color_new == color_ref:
                        tolerable.append(mass)
                        indices.append(j)
                if len(tolerable) == 0:
                    print(f"State {state} has no close state in {tag}")
                    continue
                index_s = np.abs(np.array(tolerable)-mass).argmin()
                index = indices[index_s]
                color = color_ref



            t0 = spectrum.t0
            if state not in dict_closest_levels.keys():

                dict_closest_levels[state] = [(tag,index,masses[index],errs[index],color,names[index],t0)]
            else:
                dict_closest_levels[state] += [(tag,index,masses[index],errs[index],color,names[index],t0)]
    return dict_closest_levels
def plot_level_bias_analysis(irrep_name,V,tmaxs,t0min,t0max,ncolor,tmax_ref,t0_ref,max_state,level_of_interest,ax):
    dict_closest_levels =    bias_analysis_level_matching(irrep_name,V,tmaxs,t0min,t0max,ncolor,tmax_ref,t0_ref,max_state)
    Spectrum_ref = Spectrum(V,t0_ref,irrep_name + f"_tmax_{tmax_ref}",create_files=False,saveAllPlots=False,saveHistPlots=False)
    masses_ref,errs_ref = Spectrum_ref.get_masses()
    dcol = color_coding_file_simple()
    t0s = range(t0min,t0max+1)

    padding = [0.02*l for l in range(len(t0s))]

    names = [i for i in range(max_state)]
    skips = Spectrum_ref.skips
    for skip in skips:
        if skip in names:
            names.remove(skip)
    index_of_interest = np.abs(np.array(names) - level_of_interest).argmin()
    mass_state_of_interest = masses_ref[index_of_interest]
    err_state_of_interest = errs_ref[index_of_interest]
    print(dict_closest_levels[index_of_interest])

    for close_state in dict_closest_levels[index_of_interest]:
        tag,index,mass,err,color,name,t0 = close_state
        colors = dcol[color]
        tmax = int(tag.split("_")[-3])
        #print(tmax)
        index_t0s = np.abs(np.array(t0s) - t0).argmin()

        x = tmax + padding[index_t0s]
        ax.errorbar(x,mass,yerr=err,fmt='o',label=name,color = colors)
        #ax.text(x,mass,name,fontsize = 8)
    irrep = irrep_name + f"_tmax_{tmax_ref}"

    path = f"{irrep}/Volume_{V}/t0{t0_ref}/StateColorFiles/state{level_of_interest}.txt"
    if os.path.exists(path):
                with open(path) as f:
                    lines = f.readlines()
                    color_ref= int(lines[0])
    index_t0 = np.abs(np.array(t0s) - t0_ref).argmin()
    ax.errorbar(tmax_ref+padding[index_t0],mass_state_of_interest,yerr=err_state_of_interest,fmt='o',label="State of interest",color = dcol[color_ref])
    ax.set_ylim(mass_state_of_interest+3*err_state_of_interest,mass_state_of_interest-3*err_state_of_interest)
    ticks = [round(s,4) for s in [mass_state_of_interest+2*err_state_of_interest,mass_state_of_interest+err_state_of_interest,mass_state_of_interest,mass_state_of_interest-err_state_of_interest,mass_state_of_interest-2*err_state_of_interest]]
    ax.set_yticks(ticks)
    #ax.text(tmax_ref+padding[index_t0],mass_state_of_interest,str(level_of_interest)+ f" state ref",fontsize = 8)
def plot_nine_levels_bias_analysis(irrep_name,V,tmaxs,t0min,t0max,ncolor,tmax_ref,t0_ref,max_state,levels_of_interest):
    fig, ax = plt.subplots(3,3,sharex=True,figsize=(15,15))
    #plt.subplots_adjust(wspace=0.2, hspace=0.2) 

    for i in range(3):
        for j in range(3):
            plot_level_bias_analysis(irrep_name,V,tmaxs,t0min,t0max,ncolor,tmax_ref,t0_ref,max_state,levels_of_interest[3*i+j],ax[i,j])
            ax[i,j].spines['right'].set_visible(False)
            ax[i,j].spines['top'].set_visible(False)
            

    plt.show()



            
               

def bias_analysis_level_matching_ref_irrep(irrep_ref,irrep_name,V,tmaxs,t0min,t0max,ncolor,t0_ref,max_state):
    dict_spectrums = {}
    for tmaxs in tmaxs:
        irrep = irrep_name + f"_tmax_{tmaxs}"
        for t0 in range(t0min,t0max+1):
            individual_tag = irrep + f"_t0_{t0}"
            dict_spectrums[individual_tag] = Spectrum(V,t0,irrep,create_files=False,saveAllPlots=False,saveHistPlots=False)
    #ref_tag = irrep_name + f"_tmax_{tmax_ref}" + f"_t0_{t0_ref}"

    ref_spectrum = Spectrum(V,t0_ref,irrep_ref,create_files=False,saveAllPlots=False,saveHistPlots=False)
    #ref_spectrum.automatic_coloring(max_state)
    names_ref = [i for i in range(max_state)]
    skips_ref = ref_spectrum.skips
    tags = list(dict_spectrums.keys())
    for skips in skips_ref:
        if skips in names_ref:
            names_ref.remove(skips)
    masses_ref,err_ref= ref_spectrum.get_masses()
    dict_closest_levels = {}
    for state in range(max_state):
        mass = masses_ref[state]
        err = err_ref[state]
        name = names_ref[state]
        path = f"{irrep_ref}/Volume_{V}/t0{t0_ref}/StateColorFiles/state{name}.txt"
        if os.path.exists(path):
                    with open(path) as f:
                        lines = f.readlines()
                        color_ref= int(lines[0])
        else:
            color_ref = 0

        for tag in tags:
            spectrum = dict_spectrums[tag]
            masses,errs = np.array(spectrum.get_masses())
            names = [i for i in range(spectrum.states)]
            skips = spectrum.skips
            for skip in skips:
                if skip in names:
                    names.remove(skip)
            index = np.abs(masses - mass).argmin()  # Find the index of the closest value
            
            path_closest_state = f"{irrep}/Volume_{V}/t0{spectrum.t0}/StateColorFiles/state{names[index]}.txt"
            if os.path.exists(path_closest_state):
                    with open(path_closest_state) as f:
                        lines = f.readlines()
                        color = int(lines[0])
            else:
                color = 0
            
            if color != color_ref:
                if state == 0:
                    print(color,color_ref)
                differences = np.abs(masses - mass)
                tolerable = []
                indices = []
                colors = []
                for j,difference in enumerate(differences):
                    path_closest_state_new = f"{irrep}/Volume_{V}/t0{spectrum.t0}/StateColorFiles/state{names[j]}.txt"
                    if os.path.exists(path_closest_state_new):
                        with open(path_closest_state_new) as f:
                            lines = f.readlines()
                            color_new = int(lines[0])
                        

                    else:
                        color_new = 0
                    colors += [color_new]
                    
                    if difference < 0.01 and color_new == color_ref:
                        
                        tolerable.append(masses[j])
                        indices.append(j)
                if len(tolerable) == 0:
                    #print(f"State {state} has no close state in {tag}")
                    for j,difference in enumerate(differences):
                        
                        path_closest_state_new = f"{irrep}/Volume_{V}/t0{spectrum.t0}/StateColorFiles/state{names[j]}.txt"
                        if os.path.exists(path_closest_state_new):
                            with open(path_closest_state_new) as f:
                                lines = f.readlines()
                                color_new = int(lines[0])
                        else:
                            color_new = 0
                        
                        if difference < 0.03 and color_new == 0:
                            tolerable.append(masses[j])
                            indices.append(j)
                    if len(tolerable) == 0:
                        continue


                    
                
                index_s = np.abs(np.array(tolerable)-mass).argmin()
                index = indices[index_s]
                color = colors[index]
               


            t0 = spectrum.t0
            if state not in dict_closest_levels.keys():

                dict_closest_levels[state] = [(tag,index,masses[index],errs[index],color,names[index],t0)]
            else:
                dict_closest_levels[state] += [(tag,index,masses[index],errs[index],color,names[index],t0)]
    return dict_closest_levels
def plot_level_bias_analysis_ref_irrep(irrep_ref,irrep_name,V,tmaxs,t0min,t0max,ncolor,t0_ref,max_state,level_of_interest,ax,dict_closest_levels,marker):
    
    Spectrum_ref = Spectrum(V,t0_ref,irrep_ref,create_files=False,saveAllPlots=False,saveHistPlots=False)
    masses_ref,errs_ref = Spectrum_ref.get_masses()
    dcol = color_coding_file_simple()
    t0s = range(t0min,t0max+1)

    padding = [0.1*l for l in range(len(t0s))]

    names = [i for i in range(max_state)]
    skips = Spectrum_ref.skips
    for skip in skips:
        if skip in names:
            names.remove(skip)
    index_of_interest = np.abs(np.array(names) - level_of_interest).argmin()
    mass_state_of_interest = masses_ref[index_of_interest]
    err_state_of_interest = errs_ref[index_of_interest]
    #print(dict_closest_levels[index_of_interest])
    states =[]
    err_states = []
   
    if index_of_interest not in dict_closest_levels.keys():
        dict_closest_levels[index_of_interest] = []
    for close_state in dict_closest_levels[index_of_interest]:

        tag,index,mass,err,color,name,t0 = close_state
        colors = dcol[color]
        tmax = int(tag.split("_")[-3])
        #print(tmax)
        index_t0s = np.abs(np.array(t0s) - t0).argmin()

        x = tmax + padding[index_t0s]
        ax.errorbar(x,mass,yerr=err,label=name,color = colors,marker = marker,markersize = 3)
        states.append(mass)
        err_states.append(err)
       
        #ax.text(x,mass,name,fontsize = 8)
    irrep = irrep_ref

    path = f"{irrep}/Volume_{V}/t0{t0_ref}/StateColorFiles/state{level_of_interest}.txt"
    if os.path.exists(path):
                with open(path) as f:
                    lines = f.readlines()
                    color_ref= int(lines[0])
    else:
        color_ref = 0
    index_t0 = np.abs(np.array(t0s) - t0_ref).argmin()
    xs = np.linspace(min(tmaxs),max(tmaxs)+1,100)
    ys =np.array([mass_state_of_interest for l in xs])
    yerrs =np.array( [err_state_of_interest for l in xs])
    ax.plot(xs,ys,color = dcol[color_ref],linestyle = '--',linewidth = 0.5)
    ax.fill_between(xs,ys-yerrs,ys+yerrs,color = dcol[color_ref],alpha = 0.15)
    #ax.errorbar(tmax_ref+padding[index_t0],mass_state_of_interest,yerr=err_state_of_interest,fmt='o',label="State of interest",color = dcol[color_ref])
    ax.set_ylim(round(mass_state_of_interest-3*err_state_of_interest,4),round(mass_state_of_interest+3*err_state_of_interest,4))
    ticks = [round(s,4) for s in [mass_state_of_interest+2*err_state_of_interest,mass_state_of_interest+err_state_of_interest,mass_state_of_interest,mass_state_of_interest-err_state_of_interest,mass_state_of_interest-2*err_state_of_interest]]
    ax.set_yticks(ticks)
    return states,err_states
    #ax.text(tmax_ref+padding[index_t0],mass_state_of_interest,str(level_of_interest)+ f" state ref",fontsize = 8)
def plot_nine_levels_bias_analysis_ref_irrep(ref_irrep,irrep_names,V,tmaxs,t0min,t0max,ncolor,t0_ref,max_state,levels_of_interest,markers):
    fig, ax = plt.subplots(3,3,sharex=True,figsize=(15,15))
    dcit_closest_levelss = []
    for l in range(len(irrep_names)):
        irrep_name = irrep_names[l]
        dcit_closest_levelss.append(bias_analysis_level_matching_ref_irrep(ref_irrep,irrep_name,V,tmaxs,t0min,t0max,ncolor,t0_ref,max_state))
    #plt.subplots_adjust(wspace=0.2, hspace=0.2) 

    for i in range(3):
        for j in range(3):
            for k in range(len(irrep_names)):
                irrep_name = irrep_names[k]
                dcit_closest_levels = dcit_closest_levelss[k]
                marker = markers[k]
                
                plot_level_bias_analysis_ref_irrep(ref_irrep,irrep_name,V,tmaxs,t0min,t0max,ncolor,t0_ref,max_state,levels_of_interest[3*i+j],ax[i,j],dict_closest_levels=dcit_closest_levels,marker = marker)
            ax[i,j].spines['right'].set_visible(False)
            ax[i,j].spines['top'].set_visible(False)
            

    plt.show()
        ## Find same state in other spectrums



            

    
            
     
    