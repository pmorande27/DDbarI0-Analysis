import os
import particle as P
def read_model_params(path_read):
    model_params = {}
    with open(path_read, 'r') as r:
        for i,line in enumerate(r.readlines()):
            if i == 0:
                continue
            if line[0] == '#':
                elems = line.split()
                JsP_str = elems[1]
                J = ""
                P = ""
                for let in JsP_str:
                    if let.isdigit() or let == ".":
                        J += let
                    else:
                        P += let
                twoJ = int(2*float(J))
                P = 1 if P == "+" else -1
                model_params[(twoJ,P)] = {}
            else:
                elements = line.split()
                key = elements[0]
                value = elements[1]
                model_params[(twoJ,P)][key] = value
    return model_params
def read_PartwialWaves(path_read):
    partial_waves = {}
    with open(path_read, 'r') as r:
        for i,line in enumerate(r.readlines()):
            if line[0] == '#':
                string = line.split()[1]
                Jstring = ''
                Pstring = ''
                for let in string:
                    if let.isdigit() or let == ".":
                        Jstring += let
                    else:
                        Pstring += let
                twoJ = int(2*float(Jstring))
                P = 1 if Pstring == '+' else -1
                partial_waves[(twoJ,P)] = {}
            else:
                string = line.split()
                name_1 = string[0]
                name_2 = string[1]
                waves = string[2].replace("[","").replace("]","").replace("'","").replace("{","").replace("}","").split(",")
                
                for i,wave in enumerate(waves):
                    upper = False
                    spin = ""
                    other = ""
                    for let in wave:
                        if let == "^":
                            continue
                        if let.isalpha():
                            upper = True
                        if not upper:
                            spin += let
                        if upper:
                            other += let
                    waves[i] = spin+"^"+other
                partial_waves[(twoJ,P)][(name_1,name_2)] = waves
    return partial_waves
def create_combinations_partial_wave(partial_waves_given_JP):
    possible_combinations = []
    total_single_channels = []
    for channel in partial_waves_given_JP.keys():
        for partial_wave in partial_waves_given_JP[channel]:
            name_1 = channel[0]
            name_2 = channel[1]
            chan = name_1 + ":" + name_2 + "/" + partial_wave
            total_single_channels.append(chan)
    for i in range(len(total_single_channels)):
        for j in range(i,len(total_single_channels)):
            

            possible_combinations.append((total_single_channels[j]+"|"+total_single_channels[i]))



        
    return total_single_channels,possible_combinations
def set_up_initial_values(path_model_params):
    model_params = read_model_params(path_model_params)
    partial_waves = read_PartwialWaves("./Data/PartialWaves.dat")
    print(model_params)
    print(partial_waves)
    path_write = "./Data/InitialValues.dat"
    with open(path_write,"w") as f:
        for JP in model_params.keys():
            total_single_channels,possible_combinations = create_combinations_partial_wave(partial_waves[JP])
            n_poles = int(model_params[JP]["n_poles"])
            poly_order = int(model_params[JP]["poly_order"])
            J = float(JP[0]/2)
            P = JP[1]
            P = "+" if P == 1 else "-"
            for n in range(n_poles):
                key = "m_pole"+str(n)+f"_{J}_{P}"
                m = input("Enter the value for {0} in ${1}^{2}$: ".format(key,J,P))
                value = (key,m)
                f.write(str(value)+"\n")
                
                key_err = "m_pole"+str(n)+"_err"+f"_{J}_{P}"
                m_err = input("Enter the value for {0} in ${1}^{2}$: ".format(key_err,J,P))
                value = (key_err,m_err)
                f.write(str(value)+"\n")
                key_fix = "m_pole"+str(n)+"_fix"+f"_{J}_{P}"
                m_fix = input("Enter the value for {0} in ${1}^{2}$: ".format(key_fix,J,P))
                value = (key_fix,m_fix)
                f.write(str(value)+"\n")
                key_limits = "m_pole"+str(n)+"_limits"+f"_{J}_{P}"
                m_limits = input("Enter the value for {0} in ${1}^{2}$: ".format(key_limits,J,P))
                value = (key_limits,m_limits)
                f.write(str(value)+"\n")

                for i,partial_wave in enumerate(total_single_channels):
                    key = f"g_{partial_wave}_pole{n}"
                    g = input("Enter the value for {0} in ${1}^{2}$: ".format(key,J,P))
                    #f.write("{0}:{1}\n".format(key,g))
                    value = (key,g)
                    
                    f.write(str(value)+"\n")
                    key_err = f"g_{partial_wave}_pole{n}_err"
                    g_err = input("Enter the value for {0} in ${1}^{2}$: ".format(key_err,J,P))
                    value = (key_err,g_err)
                    f.write(str(value)+"\n")
                    key_fix = f"g_{partial_wave}_pole{n}_fix"
                    g_fix = input("Enter the value for {0} in ${1}^{2}$: ".format(key_fix,J,P))
                    value = (key_fix,g_fix)
                    f.write(str(value)+"\n")
            for i in range(poly_order+1):
                for wave_combination in possible_combinations:
                    key = f"gamma_{wave_combination}_order{i}"
                    b = input("Enter the value for {0} in ${1}^{2}$: ".format(key,J,P))
                    value = (key,b)
                    f.write(str(value)+"\n")
                    key_err = f"gamma_{wave_combination}_order{i}_err"
                    b_err = input("Enter the value for {0} in ${1}^{2}$: ".format(key_err,J,P))
                    value = (key_err,b_err)
                    f.write(str(value)+"\n")
                    key_fix = f"gamma_{wave_combination}_order{i}_fix"
                
                    b_fix = input("Enter the value for {0} in ${1}^{2}$: ".format(key_fix,J,P))
                    value = (key_fix,b_fix)
                    f.write(str(value)+"\n")
                    

            #print(key,model_params[JP][key])
def fit_parameters(fit_params):
    path = "./Data"
    if os.path.exists(path) == False:
        os.mkdir(path)
    filename = path+"/FitParameters.dat"
    with open(filename,"w") as f:
        f.write("# Fit Parameters\n")
        for keys in fit_params.keys():
            f.write("# {0}:{1}\n".format(keys,fit_params[keys]))
def LatticeParameters(Vs,xi,xi_err):
    path = "./Data"
    if os.path.exists(path) == False:
        os.mkdir(path)
    filename = path+"/LatticeParameters.dat"

    with open(filename,"w") as f:
        f.write("# V xi xi_err\n")

        for V in Vs:
        
            f.write("{0} {1} {2}\n".format(V,xi,xi_err))
#LatticeParameters([16,20,24],3.444,0.006)

def model_parameters(model_params):
    path = "./Data"
    if os.path.exists(path) == False:
        os.mkdir(path)
    filename = path+"/ModelParameters.dat"
    with open(filename,"w") as f:
        f.write("# Model Parameters\n")
        for Js in model_params.keys():
            f.write("# {0}\n".format(Js))
            for key in model_params[Js].keys():
                f.write("{0} {1}\n".format(key,model_params[Js][key]))
def Hadrons():
    path = "./Data"
    if os.path.exists(path) == False:
        os.mkdir(path)
    filename = path+"/Hadrons.dat"
    all_particles = []
    path_particles = "./Particles"
    for particle_file in os.listdir(path_particles):
        all_particles += P.read_particles(path_particles+"/"+particle_file)
    with open(filename,"w") as f:
        for p in all_particles:
        
            f.write("{0} {1} {2} {3} {4}\n".format(p.name,p.Mass,p.J,p.Parity,p.mass_err))

def SpectrumInfo(SpectrumbyIrrepsbyVolume,Volumes,Irreps,excluded_levels,threshold,paths):
    path = "./Data"
    if os.path.exists(path) == False:
        os.mkdir(path)
    filename = path+"/Spectrum_formatted.dat"
    with open(filename,"w") as f:
        #   f.write("# Volume Irrep Level Energy\n")
        for irrep in Irreps:
            for V in Volumes:
                spectrum = SpectrumbyIrrepsbyVolume[irrep][V]
                exclude_levels = excluded_levels[irrep][V]
                states = spectrum.states
                skips = spectrum.skips
                real_states = [i for i in range(states) if i not in skips]
                masses,errs = spectrum.get_masses()
                for i in range(len(real_states)):
                    if i not in exclude_levels:
                        if masses[i] < threshold:
                            irrep_name = irrep[:-1]
                            momentum = "000" # need to fix this for general case
                            num = i
                            t0 = spectrum.t0
                            path_data = paths[irrep][V] + f"/t0{t0}/MassJackFiles/mass_t0_{t0}_reorder_state{num}.jack"

                            f.write("{0} {1} {2}  {3} {4}\n".format(V,momentum,irrep,num,path_data))
    filename = path+"/paths.dat"

    with open(filename,"w") as f:
        #f.write("# Volume Irrep Level Energy\n")
        for irrep in Irreps:
            for V in Volumes:
                spectrum = SpectrumbyIrrepsbyVolume[irrep][V]
                exclude_levels = excluded_levels[irrep][V]
                states = spectrum.states
                skips = spectrum.skips
                real_states = [i for i in range(states) if i not in skips]
                masses,errs = spectrum.get_masses()
                for i in range(len(real_states)):
                    if i not in exclude_levels:
                        if masses[i] < threshold:
                            irrep_name = irrep[:-1]
                            momentum = "000" # need to fix this for general case
                            num = i
                            t0 = spectrum.t0
                            path_data = paths[irrep][V] + f"/t0{t0}/MassJackFiles/mass_t0_{t0}_reorder_state{num}.jack"

                            f.write("{0} \n".format(path_data))

def pairsofHadrons(channel,included_indices,threshold):
    path = "./Data"
    if os.path.exists(path) == False:
        os.mkdir(path)
    filename = path+"/PairsOfHadrons.dat"
    chan,_ = P.channels(channel,threshold)
    Ds = P.Ds(channel,threshold)
    for key in Ds.keys():
        if key in chan.keys():
            chan[key] = chan[key] + Ds[key]
        else:
            chan[key] = Ds[key]
    dic_channels = dict(sorted(chan.items(), key=lambda item: item[1]))
    print(dic_channels)
    with open(filename,"w") as f:
        for i,key in enumerate(dic_channels.keys()):
                if i in included_indices:
                    name_1 = key[0]
                    name_2 = key[1]
                    """formatted_name_1 = ''
                    for letter in name_1:
                        if name_1 !="'":
                            formatted_name_1 += letter
                        else:
                            formatted_name_1 += ""
                    formatted_name_2 = ''
                    for letter in name_2:
                        if name_2 !="'":
                            formatted_name_2 += letter
                        else:
                            formatted_name_2 += """""


                    f.write("{0} {1}\n".format(name_1,name_2))



def set_up_partial_waves(channel,Jmax,threshold,target_Js,target_pairs):
    path = "./Data"
    if os.path.exists(path) == False:
        os.mkdir(path)
    filename = path+"/PartialWaves.dat"
    
    partial_waves_dic =P.partial_waves(channel,Jmax,threshold)
    t_pairs = []
    with open(target_pairs,"r") as f:
        for line in f:
            name_1 = line.split()[0] 
            name_2 = line.split()[1]
            

            t_pairs.append("$"+name_1+' '+name_2+"$")
    # Repeat for particle pairs with additional symmetry constraints
    Ds_partial_waves_dic = P.Ds_partial_waves(channel,Jmax,threshold)
    for i in Ds_partial_waves_dic.keys():
        if i in partial_waves_dic.keys():
            partial_waves_dic[i].update(Ds_partial_waves_dic[i])
        else:
            partial_waves_dic[i] = Ds_partial_waves_dic[i]
    new_dict = {}
    for key in partial_waves_dic.keys():
        if key in target_Js:
            for pairs in partial_waves_dic[key].keys():
                if pairs in t_pairs:
                    if key not in new_dict.keys():
                        new_dict[key] = {}
                    pair_format = pairs.replace("$","")
                    elements = pair_format.split()
                    tuple = (elements[0],elements[1])
                    new_str = elements[0] + " " + elements[1]
                    new_dict[key][new_str] = partial_waves_dic[key][pairs]
    
    with open(filename,"w") as f:
        for key in new_dict.keys():
            f.write("# {0}\n".format(key))
            for pair in new_dict[key].keys():
                waves = new_dict[key][pair]
                for i,wave in enumerate(waves):
                    waves[i] = wave.replace("$","")
                    
                f.write("{0} {1}\n".format(pair,waves))

            

def main():
    channel = {'Charm':0,'Strange':0,'Isospin':0,'Charm_Isospin':1,'C_parity':-1}
    target_pairs = "./Data/PairsOfHadrons.dat"
    target_JS = ["1.0-"]
    Hadrons()
    pairsofHadrons(channel,[0,1],0.74)
    set_up_partial_waves(channel,1,0.74,target_JS,target_pairs)
    model_parameters({"1.0-":{'n_poles':2,'poly_order':-1,"chew_man":"true"
    ,"chew_man_sub_point": "pole0"}})
   
#set_up_initial_values("./Data/ModelParameters.dat")

#set_up_initial_values("./Data/ModelParameters.dat")
