import numpy as np
import matplotlib.pyplot as plt
import particle as p
import grouptheory as gt
def norm(P):
    return np.sqrt(P[0]**2+P[1]**2+P[2]**2)
def no_int_e_level_p(m_1,P,xi,L):
    return np.sqrt(m_1**2+ (2*np.pi/(L*xi))**2*norm(P)**2)
def non_int_elevel_at_rest(name_1,name_2,P,xi,Ls):
    all_particles = p.read_particles('particles_unfl.txt')+p.read_particles('charmonium.txt')+p.read_particles('Ds.txt')
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
        


    
