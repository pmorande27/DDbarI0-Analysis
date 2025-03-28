
import sys
def read_fit_params(path_read):
    with open(path_read, 'r') as r:
        fit_params = {}
        for i,line in enumerate(r.readlines()):
            if i == 0:
                continue
            if line[0] == '#':
                elems = line.split("#")[1].split(':')

                key = elems[0].strip()
                value = elems[1]
                fit_params[key] = float(value)
    return fit_params
def create_splines_xml(path_write, path_read,path_read_lattice_params,ell_max):
    Vs = []
    try:
        fit_params = read_fit_params("./Data/FitParameters.dat")
        EcmMax = fit_params["EcmMax"]
    except:
        EcmMax = 0.25

    
    with open(path_read_lattice_params, 'r') as r_params:
        for i,line in enumerate(r_params.readlines()):
            if i != 0:
                elements = line.split()
                V = float(elements[0])
                xi = float(elements[1])
                xi_err = float(elements[2])
                Vs.append(V)
    

    
    with open(path_read, 'r') as r:
        with open(path_write, 'w') as f:
            f.write('<?xml version="1.0"?>\n')

            f.write('<MakeSplines>\n')
            f.write('   <control>\n')
            f.write(f'       <Ecm_max>{EcmMax}</Ecm_max>\n')
            f.write('       <n_points>200</n_points>\n')
            f.write('       <nsq_exp>50</nsq_exp>\n')
            f.write('       <nsq_int>50</nsq_int>\n')
            f.write('   </control>\n')
            f.write('   <lattices>\n')
            f.write('       <xi>{0}</xi>\n'.format(xi))
            Vstr = ''
            for V in Vs:
                
                Vstr += '{0} '.format(V)
            Vstr.strip()
            f.write('       <Ls>{0}</Ls>\n'.format(Vstr))
            f.write('   </lattices>\n')
            f.write('   <hadron_pairs>\n')
            for line in r.readlines():
                f.write('       <elem>\n')
                elements = line.split()
                name_one = elements[0]
                name_two = elements[1]
                f.write('           <hadron1>{0}</hadron1>\n'.format(name_one))
                f.write('           <hadron2>{0}</hadron2>\n'.format(name_two))
                f.write('           <ell_max>{0}</ell_max>\n'.format(ell_max))
                f.write('           <ds>000</ds>\n') # need to fix this for general case not at rest
                f.write('       </elem>\n')
            f.write('   </hadron_pairs>\n')
            f.write('</MakeSplines>\n')
def create_hadron_xml(path_write, path_read):
    with open(path_read, 'r') as r:
        with open(path_write, 'w') as f:
            f.write('<?xml version="1.0"?>\n')
            f.write('<Hadrons>\n')
            for line in r.readlines():
                f.write("<elem>\n")
                elements = line.split()
                name = elements[0]
                mass = float(elements[1])
                mass_err = float(elements[4])
                twoS = int(2*float(elements[2]))
                parity = int(elements[3])
                f.write(f"  <name>{name}</name>\n")
                f.write(f"  <twoS>{twoS}</twoS>\n")
                f.write(f"  <parity>{parity}</parity>\n")
                f.write(f"  <mass>{mass}</mass>\n")
                f.write(f"  <mass_err>{mass_err}</mass_err>\n")
                f.write("</elem>\n")

        
            f.write('</Hadrons>\n')
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
def create_amplitude_xml(path_read,path_write,model,model_params_path):
    names = {0:'S',1:'P',2:'D',3:'F',4:'G',5:'H',6:'I',7:'J',8:'K',9:'L',10:'M',11:'N',12:'O'}
    max_l =4
    reversed_names = {v: k for k, v in names.items()}
    model_params = read_model_params(model_params_path)
    with open(path_read, 'r') as r:
        with open(path_write, 'w') as f:

            
            f.write('<?xml version="1.0"?>\n')
            f.write('<Amplitudes>\n')

            f.write('   <PartialWaves>\n'.format(model))
            count = 0
            lines = r.readlines()
            for i,line in enumerate(lines):
                if line[0] == '#':
                    if count != 0:
                        f.write("               </Channels>\n")
                        key = (twoJ,P)
                        n_poles = model_params[key]["n_poles"]
                        poly_order = model_params[key]["poly_order"]
                        chew_man = model_params[key]["chew_man"]
                        chew_man_sub_point = model_params[key]["chew_man_sub_point"]
                        f.write(f"               <n_poles>{n_poles}</n_poles>\n")
                        f.write(f"               <poly_order>{poly_order}</poly_order>\n")
                        f.write(f"               <chew_man>{chew_man}</chew_man>\n")
                        f.write(f"               <chew_man_sub_point>{chew_man_sub_point}</chew_man_sub_point>\n")
                        f.write(f"          </fixed_params>\n")
                        f.write('       </elem>\n')
                    f.write('       <elem>\n')
                    count += 1
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
                    f.write(f"          <twoJ>{twoJ}</twoJ>\n")
                    f.write(f"          <P>{P}</P>\n")
                    f.write(f"          <model>{model}</model>\n")
                    f.write(f"          <fixed_params>\n")
                    f.write("               <Channels>\n")
                else:
                    elements = line.split()
                    name_one = elements[0]
                    name_two = elements[1]
                    partial_waves = elements[2]
                    partial_waves = partial_waves.replace("[","")
                    partial_waves = partial_waves.replace("]","")
                    partial_waves_list = partial_waves.split(",")
                    for wave in partial_waves_list:
                        wave = wave.replace("'","")
                        twoS = int(wave[1])-1
                        start = False
                        Finish = False
                        J = ""
                        for let in wave:

                            if let == "{":
                                start = True
                                continue
                            elif let == "}":
                                Finish = True
                                continue
                            elif start and not Finish:
                                J+=let      
                        two_J = int(2*float(J))
                        lstring = wave[2]
                        l = reversed_names[lstring]


                        f.write("                   <elem>\n")
                        
                        f.write(f"                       <hadron1>{name_one}</hadron1>\n")
                        f.write(f"                       <hadron2>{name_two}</hadron2>\n")
                        f.write(f"                       <twoS>{twoS}</twoS>\n")
                        f.write(f"                       <ell>{l}</ell>\n")
                        f.write(f"                       <twoJ>{two_J}</twoJ>\n")
                        f.write("                   </elem>\n")  


                if i == len(lines) - 1:
                    key = (twoJ,P)
                    f.write("               </Channels>\n")
                    n_poles = model_params[key]["n_poles"]
                    poly_order = model_params[key]["poly_order"]
                    chew_man = model_params[key]["chew_man"]
                    chew_man_sub_point = model_params[key]["chew_man_sub_point"]
                    f.write(f"               <n_poles>{n_poles}</n_poles>\n")
                    f.write(f"               <poly_order>{poly_order}</poly_order>\n")
                    f.write(f"               <chew_man>{chew_man}</chew_man>\n")
                    f.write(f"               <chew_man_sub_point>{chew_man_sub_point}</chew_man_sub_point>\n")
                    f.write(f"          </fixed_params>\n")
                    f.write('       </elem>\n')
                

                        
                
               
            f.write('   </PartialWaves>\n')
            f.write('</Amplitudes>\n')
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
def read_parameter_initial_values(path_read):
    with open(path_read, 'r') as r:
        params = {}
        for i,line in enumerate(r.readlines()):
           
            
            elems = line.replace("\n","").replace("(","").replace(")","").replace("'","").replace("\\\\","\\").split(",")
            key = elems[0]
            value = elems[1]
            params[key] = value
    print(params)
    return params
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


    print(possible_combinations)

        
    return total_single_channels,possible_combinations
            
def create_amplitude_params_xml(path_read_partial_waves,path_write,path_read_params):
    
    model_params = read_model_params(path_read_params)
    partial_waves = read_PartwialWaves(path_read_partial_waves)
    intial_values = read_parameter_initial_values("./Data/InitialValues.dat")

    with open(path_read_partial_waves, 'r') as r:
        with open(path_write, 'w') as f:
            f.write('<?xml version="1.0"?>\n')
            f.write('<Params>\n')
            f.write('   <PartialWaveParams>\n')
            count  = 0
            lines = r.readlines()
            for keys in partial_waves.keys():
                f.write('       <elem>\n')
                f.write('          <twoJ>{0}</twoJ>\n'.format(keys[0]))
                f.write('          <P>{0}</P>\n'.format(keys[1]))
                f.write('          <params>\n')
                total_partial_waves,combintation_partial_waves = create_combinations_partial_wave(partial_waves[keys])
                n_poles = model_params[keys]["n_poles"]
                poly_order = model_params[keys]["poly_order"]
                J = keys[0]/2
                P = "+" if keys[1]== 1 else "-"
                for n_pole in range(int(n_poles)):
                    key_m = f"m_pole{n_pole}_{J}_{P}"
                    key_m_err = f"m_pole{n_pole}_err_{J}_{P}"
                    key_m_fixed = f"m_pole{n_pole}_fix_{J}_{P}"
                    key_m_limits = f"m_pole{n_pole}_limits_{J}_{P}"
                    f.write('               <elem>\n')
                    f.write(f"                  <name>m_pole{n_pole}</name>\n")
                    f.write(f"                  <value>{intial_values[key_m]}</value>\n")
                    f.write(f"                  <error>{intial_values[key_m_err]}</error>\n")
                    f.write(f"                  <limits>{intial_values[key_m_limits]}</limits>\n")
                    f.write(f'                  <fixed>{intial_values[key_m_fixed]}</fixed>\n')

                    f.write('               </elem>\n')
                    for i,partial_wave in enumerate(total_partial_waves):
                        key_g = f"g_{partial_wave}_pole{n_pole}"
                        key_g_err = f"g_{partial_wave}_pole{n_pole}_err"
                        key_g_fixed = f"g_{partial_wave}_pole{n_pole}_fix"
                        f.write('               <elem>\n')
                        f.write(f"                  <name>g_{partial_wave}_pole{n_pole}</name>\n")
                        f.write(f"                  <value>{intial_values[key_g]}</value>\n")
                        f.write(f"                  <error>{intial_values[key_g_err]}</error>\n")

                        f.write(f'                  <fixed>{intial_values[key_g_fixed]}</fixed>\n')
                        f.write('               </elem>\n')
                for order in range(int(poly_order)+1):
                    for i,partial_wave_combo in enumerate(combintation_partial_waves):
                        key_gamma = f"gamma_{partial_wave_combo}_order{order}"
                        key_gamma_err = f"gamma_{partial_wave_combo}_order{order}_err"
                        key_gamma_fixed = f"gamma_{partial_wave_combo}_order{order}_fix"
                        f.write('               <elem>\n')
                        f.write(f"                  <name>gamma_{partial_wave_combo}_order{order}</name>\n")
                        f.write(f"                  <value>{intial_values[key_gamma]}</value>\n")
                        f.write(f"                  <error>{intial_values[key_gamma_err]}</error>\n")
                        f.write(f'                  <fixed>{intial_values[key_gamma_fixed]}</fixed>\n')
                        f.write('               </elem>\n')
                    
          
         


                f.write('          </params>\n')
                f.write('       </elem>\n')

                

            f.write('</PartialWaveParams>\n')
           
            f.write('</Params>\n')
    

def main(args):
    if args[1] == 'create_splines_xml':
        path_read = "./Data/PairsOfHadrons.dat"
        path_write = "./splines/Splines.xml"
        path_lattice_params = "./Data/LatticeParams.dat"
        create_splines_xml(path_write,path_read,path_lattice_params,4)
    elif args[1] == 'create_hadron_xml':
        path_read = "./Data/Hadrons.dat"
        path_write = "./hadrons/Hadrons.xml"

        create_hadron_xml(path_write,path_read)
    elif args[1] == 'create_amplitude_xml':
        path_read = "./Data/PartialWaves.dat"
        path_write = "./amplitudes/Amplitude.xml"
        model = "kmatrix_poles_poly_fast"
        model_params = "./Data/ModelParameters.dat"
        create_amplitude_xml(path_read,path_write,model,model_params)
    elif args[1] == 'create_amplitude_params_xml':
        path_read = "./Data/PartialWaves.dat"
        path_write = "./amplitudes/Params.xml"
        model_params = "./Data/ModelParameters.dat"
        create_amplitude_params_xml(path_read,path_write,model_params)
    else:
        print("Invalid command")
if __name__ == '__main__':
    print(read_fit_params("./Data/FitParameters.dat"))
    #main(sys.argv)
path_read = "./Data/PartialWaves.dat"
path_write = "./Data/Params.xml"
model_params = "./Data/ModelParameters.dat"
#create_amplitude_params_xml(path_read,path_write,model_params)
