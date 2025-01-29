import numpy as np
import matplotlib.pyplot as plt
import particle as p
import grouptheory as gt
import os
import colors
def norm(P):
    return np.sqrt(P[0]**2+P[1]**2+P[2]**2)
def no_int_e_level_p(m_1,P,xi,L):
    return np.sqrt(m_1**2+ (2*np.pi/(L*xi))**2*norm(P)**2)
def non_int_elevel_at_rest(name_1,name_2,P,xi,Ls):
    all_particles = p.read_particles('Particles/particles_unfl.txt')+p.read_particles('Particles/charmonium.txt')+p.read_particles('Particles/Ds.txt')
    count = 0
    for particle in all_particles:
        if particle.name == name_1:
            particle_1 = particle
            count += 1
        if particle.name == name_2:
            particle_2 = particle
            count += 1
    mass_1 = float(particle_1.Mass)
    mass_2 = float(particle_2.Mass)
    e_1 = np.array([no_int_e_level_p(mass_1,P,xi,L) for L in Ls])
    e_2 = np.array([no_int_e_level_p(mass_2,-P,xi,L) for L in Ls])
    E_cm_sq = (e_1+e_2)**2 -  np.array([(2*np.pi/(L*xi))**2 *norm(P-P)**2 for L in Ls])

    return np.sqrt(E_cm_sq)

#print(non_int_elevel_at_rest('\eta_c','\eta',np.array([1,0,0]),3.444,np.linspace(16,24)))
def identify_momentum_type(P):
    if P[0]==0 and P[1] == 0 and P[2] == 0:
        return 'O_D^h'
    elif P[1] == 0 and P[2] == 0:
        return 'Dic_4'
    elif P[2] == 0 and P[1] == P[0]:
        return 'Dic_2'
    elif P[2] == P[1] and P[1] == P[0]:
        return 'Dic_3'
    else:
        raise ValueError('Momentum type not recognized')
def E_level_in_irrep(irrep_name,name_1,name_2,Ls,P,channel,Jmax,threshold):
    xi = 3.444
    Flag = False
    if 'D' in name_1 and 'D' in name_2:
        if name_1 == p.bar(name_2):
            Flag = True
            general_allowed_irreps_with_symmetry = p.possible_irreps_rest(name_1,name_2,channel,Jmax,threshold)

    group = identify_momentum_type(P)
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
def E_levels_in_irrep_P(irrep_name,channels,Ls,P,channel,Jmax,threshold):
    Es = []
    xi = 3.444
    indeces = []
    multiplicities = []
    for i,(name_1,name_2) in enumerate(channels):
        statement,multiplicity = E_level_in_irrep(irrep_name,name_1,name_2,Ls,P,channel,Jmax,threshold)
        if statement:
            multiplicities.append(multiplicity)
            indeces.append(i)
            Es.append(non_int_elevel_at_rest(name_1,name_2,P,xi,Ls))
    return Es,indeces,multiplicities
def plot_irrep(irrep_name,Ls,Ps):
    channel = {'Charm':0,'Strange':0,'Isospin':0,'Charm_Isospin':1,'C_parity':-1}
    colors = ['red','blue','green','orange','purple','pink','black','brown','teal','cyan','magenta','grey','lime','olive','yellow','navy','maroon','aqua','fuchsia','silver','red','blue','green','orange','purple','pink','black','brown','yellow','cyan','magenta','grey','lime','olive','teal','navy','maroon','aqua','fuchsia','silver']
    C_p = '+' if channel['C_parity'] == 1 else '-'
    threshold = 0.72
    channelsDs = p.Ds(channel,threshold)
    chh = {}#p.channels(channel,threshold)[0]
    for key in channelsDs.keys():
        chh[key] = channelsDs[key]
    dict = {}
    channels = chh.keys()
    for i,chan in enumerate(channels):
        dict[chan] = colors[i]
    dict[('\\eta_c','\\eta')] = 'red'
    dict[('\\eta_c','\\eta\'')] = 'pink'
    dict['D','\\bar{D}'] = 'green'
    dict['D_s','\\bar{D}_s'] = 'lime'
    dict['D^*','\\bar{D}^*'] = 'blue'
    #dict['D_s^*','\\bar{D}_s^*'] = 'cyan'
    dict[('\\psi','w')] = 'magenta'
    dict[('\\psi','\phi')] = 'violet'
    dict[('D','\\bar{D}^*')] = 'grey'
    dict[('\\eta_c','\\sigma')] = 'maroon'
    dict[('\\psi','\\sigma')] = 'grey'

    

    


    dict2 = {}
    for i,chan in enumerate(channels):
        dict2[i] = chan
    ch = list(channels)
    counts = {}
    for i,chan in enumerate(channels):
        counts[chan] = 0
    labels = []
    for P in Ps:

        Es,indeces,multiplcities = E_levels_in_irrep_P(irrep_name,channels,Ls,P,channel,4,threshold)
        # Assign a color to each channel
        #colors = plt.cm.get_cmap('hsv', len(ch))

        for i,E in enumerate(Es):
            counts[ch[indeces[i]]] += 1

            for j in range(multiplcities[i]):
                
                padding = np.array([0.001 for l in E])
                label = '$'+ ch[indeces[i]][0] + ' ' +ch[indeces[i]][1] +'$'
                if label not in labels:
                    plt.plot(Ls,E+j*padding,color =dict[dict2[indeces[i]]],linewidth = 0.7,label = label)

                    labels.append(label)
                else:
                    plt.plot(Ls,E+j*padding,color = dict[dict2[indeces[i]]],linewidth = 0.7)

    plt.xlim(14,26)
    plt.ylim(0.6,0.74)
    plt.xticks([16,20,24])
    plt.yticks([0.62,0.64,0.66,0.68,0.70,0.72,0.74])

    for i,ch in enumerate(channels):
        if counts[ch] != 0:
            E = non_int_elevel_at_rest(ch[0],ch[1],np.array([0,0,0]),3.444,Ls)

            plt.plot(Ls,E,color = dict[dict2[i]],linewidth = 0.7,linestyle = '--')

    plt.legend()
    title = 'Irrep: $ '+irrep_name.split('^')[0] + "^{"+irrep_name.split('^')[1] +C_p+'}$'
    plt.title(title)
    plt.savefig('noint_'+irrep_name+'.svg')
    plt.show()

#[plot_irrep(i,np.linspace(14,26),[np.array([0,0,0]),np.array([1,0,0]),np.array([1,1,0]),np.array([1,1,1])]) for i in ['T_1^-']]
        

def get_levels_irrep_in_flight(channel,irrep_name,file_name,momentum):
    channel_name = ""
    irrep_momentum = momentum+"_"+irrep_name
    if momentum == '000':
        for i, key in enumerate(channel.keys()):
            

            if i == len(channel.keys())-1:
                channel_name += str(key) +" " + str(channel[key])
            else:
                channel_name += str(key) +" " + str(channel[key]) + ', '
        path = "moving_frames/"+channel_name + '/' + file_name
        dic = {}
        with open(path,'r') as file:
            lines = file.readlines()
            for l in lines:
                target_irrep = l.split("  ")[1].strip()
                if target_irrep == irrep_momentum:
                    try: 
                        ops = l.split("   ")[3].strip()
                        op1 = ops.split("xx")[0]
                        
                        op2 = ops.split("xx")[1]
                    except:
                        ops = l.split("   ")[4].strip()
                        op1 = ops.split("xx")[0]
                        op2 = ops.split("xx")[1]
                    p1 = colors.identify_particle_plus_plus(op1)
                    p2 = colors.identify_particle_plus_plus(op2)
                    energy = float(l.split(target_irrep)[1].strip().split("  ")[1].strip())
                    mom1 = op1.split('__')[1]
                    mom2 = op2.split('__')[1]
                    
                    if "f_proj0" not in op1 and  "f_proj0" not in op2:
                        if (p1,p2,mom1,mom2) not in dic.keys():
                            dic[(p1,p2,mom1,mom2)] = 1
                        else:
                            dic[(p1,p2,mom1,mom2)] += 1
        
    if momentum != '000':
        for i, key in enumerate(channel.keys()):
            

            if i == len(channel.keys())-1:
                channel_name += str(key) +" " + str(channel[key])
            else:
                channel_name += str(key) +" " + str(channel[key]) + ', '
        path = "moving_frames/"+channel_name + '/' + file_name
        dic = {}
        with open(path,'r') as file:
            lines = file.readlines()
            for l in lines:
                target_irrep = l.split("  ")[1].strip()
                if target_irrep == irrep_momentum:
                    ops = l.split("   ")[4].strip()
                    op1 = ops.split("xx")[0]
                    op2 = ops.split("xx")[1]
                    p1 = colors.identify_particle_plus_plus(op1)
                    p2 = colors.identify_particle_plus_plus(op2)
                    energy = float(l.split(target_irrep)[1].strip().split("  ")[1].strip())
                    mom1 = op1.split('__')[1]
                    mom2 = op2.split('__')[1]
                    if "f_proj0" not in op1 and "f_proj0" not in op2:
                        if (p1,p2,mom1,mom2) not in dic.keys():
                            dic[(p1,p2,mom1,mom2)] = 1
                        else:
                            dic[(p1,p2,mom1,mom2)] += 1
    return dic
def get_levels_irrep_in_flight_sigma(channel,irrep_name,file_name,momentum):
    channel_name = ""
    irrep_momentum = momentum+"_"+irrep_name
    if momentum == '000':
        for i, key in enumerate(channel.keys()):
            

            if i == len(channel.keys())-1:
                channel_name += str(key) +" " + str(channel[key])
            else:
                channel_name += str(key) +" " + str(channel[key]) + ', '
        path = "moving_frames/"+channel_name + '/' + file_name
        dic = {}
        with open(path,'r') as file:
            lines = file.readlines()
            for l in lines:
                target_irrep = l.split("  ")[1].strip()
                if target_irrep == irrep_momentum:
                    try: 
                        ops = l.split("   ")[3].strip()
                        op1 = ops.split("xx")[0]
                        
                        op2 = ops.split("xx")[1]
                    except:
                        ops = l.split("   ")[4].strip()
                        op1 = ops.split("xx")[0]
                        op2 = ops.split("xx")[1]
                    p1 = colors.identify_particle_plus_plus(op1)
                    p2 = colors.identify_particle_plus_plus(op2)
                    energy = float(l.split(target_irrep)[1].strip().split("  ")[1].strip())
                    mom1 = op1.split('__')[1]
                    mom2 = op2.split('__')[1]
                    
                    if "f_proj0" in op1 or  "f_proj0" in op2:
                        if (p1,p2,mom1,mom2) not in dic.keys():
                            dic[(p1,p2,mom1,mom2)] = 1
                        else:
                            dic[(p1,p2,mom1,mom2)] += 1
        
    if momentum != '000':
        for i, key in enumerate(channel.keys()):
            

            if i == len(channel.keys())-1:
                channel_name += str(key) +" " + str(channel[key])
            else:
                channel_name += str(key) +" " + str(channel[key]) + ', '
        path = "moving_frames/"+channel_name + '/' + file_name
        dic = {}
        with open(path,'r') as file:
            lines = file.readlines()
            for l in lines:
                target_irrep = l.split("  ")[1].strip()
                if target_irrep == irrep_momentum:
                    ops = l.split("   ")[4].strip()
                    op1 = ops.split("xx")[0]
                    op2 = ops.split("xx")[1]
                    p1 = colors.identify_particle_plus_plus(op1)
                    p2 = colors.identify_particle_plus_plus(op2)
                    energy = float(l.split(target_irrep)[1].strip().split("  ")[1].strip())
                    mom1 = op1.split('__')[1]
                    mom2 = op2.split('__')[1]
                    #print(op1,op2)
                    if "f_proj0" in op1 or "f_proj0" in op2:

                        if (p1,p2,mom1,mom2) not in dic.keys():
                            dic[(p1,p2,mom1,mom2)] = 1
                        else:
    
                            dic[(p1,p2,mom1,mom2)] += 1
    return dic
def get_E_levels_irrep_in_flight_sigma(channel,irrep_name,momentum):
    file_name = "energies_no_sym.txt"
    levels = {}
    xi = 3.444

    dic = get_levels_irrep_in_flight_sigma(channel,irrep_name,file_name,momentum)
    additional_levels_file = "energies_sym.txt"
    additional_dic = get_levels_irrep_in_flight_sigma(channel,irrep_name,additional_levels_file,momentum)
    all_particles = p.read_particles('Particles/particles_unfl.txt')+p.read_particles('Particles/charmonium.txt')+p.read_particles('Particles/Ds.txt')

    for key in additional_dic.keys():
        if key in dic.keys():
            dic[key] += additional_dic[key]
        else:
            dic[key] = additional_dic[key]
    for key in dic.keys():
        p1_name = key[0]
        p2_name = key[1]
        mom1 = key[2]
        mom2 = key[3]
        multiplicity = dic[key]
        if p2_name == "f_0":
            sigma = "\sigma"
            mom_sigma = mom2
            no_sigma = p1_name
            mom_no_sigma = mom1
        else:
            sigma = "\sigma"
            mom_sigma = mom1
            no_sigma = p2_name
            mom_no_sigma = mom2
        for particle in all_particles:
            if particle.name == sigma:
                sigma_particle = particle
            if particle.name == no_sigma:

                no_sigma_particle = particle
        possible_volumes = [16,20,24]
        
        realised_volumes = []
        realised_Energies = []
        realised_Errors = []
        for volume in possible_volumes:
            path = f"sigma/sigma_E_levels_{mom_sigma}_{volume}.txt"
            if os.path.exists(path):
                realised_volumes.append(volume)
                with open(path) as f:
                        sigma_Es = []
                        err_Es = []
                        for l,line in enumerate(f):
                            P_sigma = get_mom_from_text(mom_sigma)
                            E_sigma_lat =np.sqrt( float(line.split('+/-')[0])**2+((2*np.pi)/(volume*xi))**2*np.linalg.norm(P_sigma)**2)
                            
                            if len(line.split('+/-')) == 2:
                                err_E = float(line.split('+/-')[1])
                            else:
                                err_E = 0
                            err_p = float(line.split('+/-')[0])/np.sqrt( float(line.split('+/-')[0])**2+((2*np.pi)/(volume*xi))**2*np.linalg.norm(P_sigma)**2)*err_E

                            particle_no_sigma_E = no_int_e_level_p(no_sigma_particle.Mass,get_mom_from_text(mom_no_sigma),3.444,volume)
                            E_tot_lat = E_sigma_lat + particle_no_sigma_E
                            total_P = get_mom_from_text(momentum)
                            E_tot_cm = np.sqrt(E_tot_lat**2 - ((2*np.pi)/(volume*xi))**2*np.linalg.norm(total_P)**2)
                            err_tot_propagated =(E_tot_lat /(np.sqrt(E_tot_lat**2 - ((2*np.pi)/(volume*xi))**2*np.linalg.norm(total_P)**2)))*err_p

                            err_Es += [err_tot_propagated]
                            sigma_Es += [E_tot_cm]
                            
                            P_no_sigma = get_mom_from_text(mom_no_sigma)

                            
                        realised_Energies.append((sigma_Es))
                        realised_Errors.append((err_Es))
        levels[key] = (realised_volumes,realised_Energies,realised_Errors,multiplicity)
    
    return levels

        


        




def particles_to_operator_name_for_xml_file(particle):
    if particle == '\pi':
        return 'pion_proj0_J0'
    if particle == 'K':
        return 'kaon_proj0_J0'
    if particle == 'D':
        return 'Dneg_proj0_J0'
    if particle == 'D_s':
        return 'Dsneg_proj0_J0'
    if particle == 'D^*':
        return 'Dneg_proj1_J1'
    if particle == 'D^*_s':
        return 'Dsneg_proj1_J1'
    if particle == '\\bar{D}':
        return 'Dbarneg_proj0_J0'
    if particle == '\\bar{D}_s':
        return 'Dsbarneg_proj0_J0'
    if particle == '\\bar{D}^*':
        return 'Dbarneg_proj1_J1'
    if particle == '\\bar{D}^*_s':
        return 'Dsbarneg_proj1_J1'
    if particle == '\eta_c':
        return 'etace_proj0_J0'
    if particle == "\eta_{c}'":
        return 'etace_proj1_J0'
    if particle == '\psi':
        return 'psice_proj0_J1'
    if particle == "\psi'":
        return 'psice_proj1_J1'
    if particle == '\eta':
        return 'eta_proj0_J0'
    if particle == "\eta'":
        return 'eta_proj1_J0'
    if particle == 'w':
        return 'omega_proj0_J1'
    if particle == "\phi":
        return 'omega_proj1_J1'
    if particle == 'f_0':
        return 'f_'
    if particle == 'h_c':
        return 'hce_proj0_J1'
    if particle == '\chi_{c0}':
        return 'chice_proj0_J0'
    if particle == '\chi_{c1}':
        return 'chice_proj1_J1'
    if particle == '\chi_{c2}':
        return 'chice_proj2_J2'
    if particle == '\sigma':
        return 'f_proj0_J0'
    else:
        print(particle)
    
        
        

    
def generate_xml_pairs_no_symmetry(channel,threshold):
    channelsDs = p.Ds(channel,threshold)
    chh = p.channels(channel,threshold)[0]
    particles = p.read_particles('Particles/Ds.txt') + p.read_particles('Particles/charmonium.txt') + p.read_particles('Particles/particles_unfl.txt')
    print("<Masses>")
    for part in particles:
        """if part.name == '\sigma':
            continue"""
        line = "<elem><Key>"+particles_to_operator_name_for_xml_file(part.name)+"</Key><Val>"+str(part.Mass)+"</Val></elem>"
        print(line)
    print("</Masses>")
    for key in channelsDs.keys():
        chh[key] = channelsDs[key]
    xml_pairs = []
    for key in chh.keys():
        if key[0] == p.bar(key[1]):
            continue
        particle1 = particles_to_operator_name_for_xml_file(key[0])
        particle2 = particles_to_operator_name_for_xml_file(key[1])
        
            

        xml_pairs.append((particle1,particle2))
    
    #print(xml_pairs)
    print("<Operators>")
    for pair in xml_pairs:
        line = "<elem><First>" + pair[0] + "</First>" + "<Second>" + pair[1] + "</Second>"+ "</elem>"
        print(line)
    print("</Operators>")

def generate_xml_pairs_with_symmetry(channel,threshold):
    channelsDs = p.Ds(channel,threshold)
    chh = p.channels(channel,threshold)[0]
    particles = p.read_particles('Particles/Ds.txt')
    print("<Masses>")
    for part in particles:
        line = "<elem><Key>"+particles_to_operator_name_for_xml_file(part.name)+"</Key><Val>"+str(part.Mass)+"</Val></elem>"
        print(line)
    print("</Masses>")
    for key in channelsDs.keys():
        chh[key] = channelsDs[key]
    xml_pairs = []
    for key in chh.keys():
        if key[0] == "\sigma" or key[1] == '\sigma':
            continue
        if key[0] != p.bar(key[1]):
            continue
        particle1 = particles_to_operator_name_for_xml_file(key[0])
        particle2 = particles_to_operator_name_for_xml_file(key[1])
        
            

        xml_pairs.append((particle1,particle2))
    
    #print(xml_pairs)
    print("<Operators>")
    for pair in xml_pairs:
        line = "<elem><First>" + pair[0] + "</First>" + "<Second>" + pair[1] + "</Second>"+ "</elem>"
        print(line)
    print("</Operators>")
    "<target_irrepmom_sym>1</target_irrepmom_sym>"
def non_int_elevel(name_1,name_2,P_1,P_2,total_P,xi,Ls):
    all_particles = p.read_particles('Particles/particles_unfl.txt')+p.read_particles('Particles/charmonium.txt')+p.read_particles('Particles/Ds.txt')
    count = 0
    for particle in all_particles:
        if particle.name == name_1:
            particle_1 = particle
            count += 1
        if particle.name == name_2:
            particle_2 = particle
            count += 1
    mass_1 = float(particle_1.Mass)
    mass_2 = float(particle_2.Mass)
    e_1 = np.array([no_int_e_level_p(mass_1,P_1,xi,L) for L in Ls])
    e_2 = np.array([no_int_e_level_p(mass_2,P_2,xi,L) for L in Ls])
    E_cm_sq = (e_1+e_2)**2 -  np.array([(2*np.pi/(L*xi))**2 *norm(total_P)**2 for L in Ls])

    return np.sqrt(E_cm_sq)
def get_mom_from_text(text):
    mom  =[]
    for i in range(len(text)):
        mom.append(int(text[i]))
    return mom

def get_E_levels_in_flight(channel,target_irrep,total_momentum,L):
    dic_levels = get_levels_irrep_in_flight(channel,target_irrep,"energies_sym.txt",total_momentum)
    additional_levels = get_levels_irrep_in_flight(channel,target_irrep,"energies_no_sym.txt",total_momentum)
    levels = {}
    for key in additional_levels.keys():
        if key in dic_levels.keys():
            dic_levels[key] += additional_levels[key]
        else:
            dic_levels[key] = additional_levels[key]
    for key in dic_levels.keys():
        Es = non_int_elevel(key[0],key[1],get_mom_from_text(key[2]),get_mom_from_text(key[3]),get_mom_from_text(total_momentum),3.444,L)
        mutliplicity = dic_levels[key]
        levels[key] = (Es,mutliplicity)
    return levels
def create_operator_list(channel,irrep_name,momentum,threshold):
    file_name = "energies_no_sym.txt"
    dic = get_operator_list_irrep_threshold(channel,irrep_name,file_name,momentum,threshold)
    additional_levels_file = "energies_sym.txt"
    additional_dic = get_operator_list_irrep_threshold(channel,irrep_name,additional_levels_file,momentum,threshold)
    for key in additional_dic.keys():
        if key in dic.keys():
            dic[key] += additional_dic[key]
        else:
            dic[key] = additional_dic[key]
    new_file_path = "ops_list/"+irrep_name+"_"+momentum+".txt"
    with open(new_file_path,'w') as file:
        for key in dic.keys():
            file.write(dic[key]+'\n')
    

def get_operator_list_irrep_threshold(channel,irrep_name,file_name,momentum,threshold):
    channel_name = ""
    irrep_momentum = momentum+"_"+irrep_name
    if momentum == '000':
        for i, key in enumerate(channel.keys()):
            

            if i == len(channel.keys())-1:
                channel_name += str(key) +" " + str(channel[key])
            else:
                channel_name += str(key) +" " + str(channel[key]) + ', '
        path = "moving_frames/"+channel_name + '/' + file_name
        dic = {}
        with open(path,'r') as file:
            lines = file.readlines()
            for l in lines:
                target_irrep = l.split("  ")[1].strip()
                if target_irrep == irrep_momentum:
                    try: 
                        ops = l.split("   ")[3].strip()
                        op1 = ops.split("xx")[0]
                        
                        op2 = ops.split("xx")[1]
                    except:
                        ops = l.split("   ")[4].strip()
                        op1 = ops.split("xx")[0]
                        op2 = ops.split("xx")[1]
                    p1 = colors.identify_particle_plus_plus(op1)
                    p2 = colors.identify_particle_plus_plus(op2)
                    irrep1 = op1.split('_')[3]
                    irrep2 = op2.split('_')[3]
                    #print(irrep1,irrep2)
                    energy = float(l.split(target_irrep)[1].strip().split("  ")[1].strip())
                    mom1 = op1.split('__')[1]
                    mom2 = op2.split('__')[1]
                    new_l = target_irrep + " " + ops
                    
                    if energy < threshold:
                        dic[(p1,p2,mom1,mom2,irrep1,irrep2)] = new_l
                   
        
    if momentum != '000':
        for i, key in enumerate(channel.keys()):
            

            if i == len(channel.keys())-1:
                channel_name += str(key) +" " + str(channel[key])
            else:
                channel_name += str(key) +" " + str(channel[key]) + ', '
        path = "moving_frames/"+channel_name + '/' + file_name
        dic = {}
        with open(path,'r') as file:
            lines = file.readlines()
            for l in lines:
                target_irrep = l.split("  ")[1].strip()
                if target_irrep == irrep_momentum:
                    ops = l.split("   ")[4].strip()
                    op1 = ops.split("xx")[0]
                    op2 = ops.split("xx")[1]
                    p1 = colors.identify_particle_plus_plus(op1)
                    p2 = colors.identify_particle_plus_plus(op2)
                    #print(op1,op2)
                    irrep1 = op1.split('_')[3]
                    irrep2 = op2.split('_')[3]
                    #print(irrep1,irrep2)
                    energy = float(l.split(target_irrep)[1].strip().split("  ")[1].strip())
                    mom1 = op1.split('__')[1]
                    mom2 = op2.split('__')[1]
                    new_l = target_irrep + " " + ops
                    if energy < threshold:
                        dic[(p1,p2,mom1,mom2,irrep1,irrep2)] = new_l
    return dic
def purge_sym_no_int_levels(channel):
    path = "moving_frames/"
    for i, key in enumerate(channel.keys()):
        if i == len(channel.keys())-1:
            path += str(key) +" " + str(channel[key])
        else:
            path += str(key) +" " + str(channel[key]) + ", "
    dic = {}
    deleted = 0
    deleted_list = []
    with open(path + "/energies_sym.txt",'r') as file:
        lines = file.readlines()
        for l in lines:
            try: 
                ops = l.split("   ")[3].strip()
                op1 = ops.split("xx")[0]
                        
                op2 = ops.split("xx")[1]
            except:
                ops = l.split("   ")[4].strip()
                op1 = ops.split("xx")[0]
                op2 = ops.split("xx")[1]
            particle1 = op1.split('_p')[0]+"_p"+op1.split('_p')[1]
            particle2 = op2.split('_p')[0]+"_p"+op1.split('_p')[1]
            irrep1 = op1.split('_p')[2].split('_')[1]
            irrep2 = op2.split('_p')[2].split('_')[1]

            condition = True
            mom1 = op1.split('_p')[2].split('_')[0]
            mom2 = op2.split('_p')[2].split('_')[0]
           
            ittep_t = l.split("   ")[1]
            tuple1 = (ittep_t,particle1,particle2,irrep1,irrep2,mom1,mom2)
            tuple2 = (ittep_t,particle1,particle2,irrep2,irrep1,mom2,mom1)
            if tuple2 not in dic.keys():
                ## Need to test this further
                if  tuple1 not in dic.keys():
                    dic[(ittep_t,particle1,particle2,irrep1,irrep2,mom1,mom2)] = [l]
                else:
                    dic[(ittep_t,particle1,particle2,irrep1,irrep2,mom1,mom2)].append(l)

            else:
                deleted += 1
                deleted_list.append(l)
    print(deleted)
    for dele in deleted_list:
        print(dele)
    with open(path + "/energies_sym.txt",'w') as file:
        for key in dic.keys():
            for l in dic[key]:
                file.write(l)

            

