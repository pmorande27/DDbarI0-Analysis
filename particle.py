import quantum_numbers
import angular_mom
import matplotlib.pyplot as plt
import grouptheory as gt
class Particle(object):
    """
    Object that contains the information of a particle, useful to package all the information
    """
    def __init__(self,name,Isospin,Isospin_charm, C_parity,Charm,Strange,Mass,Parity,J):
        """
        Constructor of the Particle object
        :param name: Name of the particle
        :param Isospin: Isospin of the particle
        :param Isospin_charm: In HadSpec we work with a degenerate pair of charm quarks, charm and e-quark, having therefore an additional SU(2)
                              symmetry and an associated conserved number called Isospin_charm. This number is 0 for the charm quark and 1 for 
                              the e-quark. This is performed to work with in a theory which effecitvely forbbids c cbar pairs to annhielate by working
                              with Isospin_charm = 0 channels.
        :param C_parity: C parity of the particle
        :param Charm: Charmness of the particle
        :param Strange: Strangeness of the particle
        :param Mass: Mass of the particle
        :param Parity: Parity of the particle
        :param J: Intrinsic Angular momentum of the particle
        """
        self.name = name
        self.Parity = Parity
        self.Isospin_charm = Isospin_charm
        self.C_parity = C_parity
        self.Charm = Charm
        self.Strange = Strange
        self.Mass = float(Mass)
        self.Isospin = Isospin
        self.J = J
    def __str__(self):
        """
        String representation of the Particle object
        :return: String representation of the Particle object
        """
        Parity = '+' if self.Parity == 1 else '-'
        if self.C_parity == 1:
            C_parity = '+'
        elif self.C_parity == -1:
            C_parity = '-'
        else:
            C_parity = ''
        return self.name + ' ' + str(self.Isospin) + ' ' +  str(self.J) + '^{'+str(Parity)+str(C_parity) + '} '
def two_particles_inchannel(particel_one,particle_two,channel_charm,channel_strange,channel_isospin,channel_charm_isospin,channel_C_parity):
    """
    Function that checks if a state formed by a pair of particles are in a given channel, 
    meaning that they can produce the quantum numbers of the channel.
    :param particel_one: Particle object corresponding to the first particle
    :param particle_two: Particle object corresponding to the second particle
    :param channel_charm: Charmness of the channel
    :param channel_strange: Strangeness of the channel
    :param channel_isospin: Isospin of the channel
    :param channel_charm_isospin: Isospin_charm of the channel
    :param channel_C_parity: C parity of the channel
    :return: True if the pair of particles are in the channel, False otherwise
    """
    # Calculate possible Isospin combinations
    Isospins = quantum_numbers.Isospin(particel_one.Isospin,particle_two.Isospin)

    # Calculate possible Isospin_charm combinations
    charm_Isospins = quantum_numbers.Isospin(particel_one.Isospin_charm,particle_two.Isospin_charm)

    # Calculate Charmness of the pair
    Charm = quantum_numbers.Charm(particel_one.Charm,particle_two.Charm)

    # Calculate Strangeness of the pair
    Strange = quantum_numbers.Strangeness(particel_one.Strange,particle_two.Strange)

    # If C_parity is conserved check for C_parity of the pair
    if particel_one.Charm == 0 and particle_two.Charm == 0 and particel_one.Strange == 0 and particle_two.Strange == 0:
        C_parity =quantum_numbers.Charge_conjugation(particel_one.C_parity,particle_two.C_parity)
        # Check if the pair of particles are in the channel
        if channel_charm== Charm and channel_strange == Strange and channel_isospin in Isospins and channel_charm_isospin in charm_Isospins and channel_C_parity == C_parity:
            return True
        else:
            return False
    # If C_parity not conserved just ignore it and check for the rest of the properties.
    else:
        if channel_charm== Charm and channel_strange == Strange and channel_isospin in Isospins and channel_charm_isospin in charm_Isospins:
            return True
        else:
            return False


def three_particles_inchannel(particel_one,particle_two,particle_three,channel_charm,channel_strange,channel_isospin,channel_charm_isospin,channel_C_parity):
    """
    Function that checks if a state formed by a pair of particles are in a given channel,
    meaning that they can produce the quantum numbers of the channel.
    :param particel_one: Particle object corresponding to the first particle
    :param particle_two: Particle object corresponding to the second particle
    :param particle_three: Particle object corresponding to the third particle
    :param channel_charm: Charmness of the channel
    :param channel_strange: Strangeness of the channel
    :param channel_isospin: Isospin of the channel
    :param channel_charm_isospin: Isospin_charm of the channel
    :param channel_C_parity: C parity of the channel
    :return: True if the pair of particles are in the channel, False otherwise
    """
    
    # Calculate possible Isospin combinations
    Isospins = quantum_numbers.Isospin_three_particles(particel_one.Isospin,particle_two.Isospin,particle_three.Isospin)
    
    # Calculate possible Isospin_charm combinations
    charm_Isospins = quantum_numbers.Isospin_three_particles(particel_one.Isospin_charm,particle_two.Isospin_charm,particle_three.Isospin_charm)
    
    # Calculate Charmness of the trio
    Charm = quantum_numbers.Charm_three_particles(particel_one.Charm,particle_two.Charm,particle_three.Charm)
    
    # Calculate Strangeness of the trio
    Strange = quantum_numbers.Strangeness_three_particles(particel_one.Strange,particle_two.Strange,particle_three.Strange)
    
    # If C_parity is conserved check for C_parity of the trio
    C_parity = quantum_numbers.Charge_conjugation_three_particles(particel_one.C_parity,particle_two.C_parity,particle_three.C_parity)
    if particel_one.Charm == 0 and particle_two.Charm == 0 and particle_three.Charm == 0 and particel_one.Strange == 0 and particle_two.Strange == 0 and particle_three.Strange == 0:
        if channel_charm== Charm and channel_strange == Strange and channel_isospin in Isospins and channel_charm_isospin in charm_Isospins and channel_C_parity == C_parity:
            return True
        else:
            return False
    
    # If C_parity not conserved just ignore it and check for the rest of the properties.
    else:
        if channel_charm== Charm and channel_strange == Strange and channel_isospin in Isospins and channel_charm_isospin in charm_Isospins:
            return True
        else:
            return False

def read_particles(filename):
    """
    Function that reads the particles from a file and creates a list of Particle objects
    :param filename: Name of the file where the particles are stored
    :return: List of Particle objects
    """
    particles = []
    with open(filename) as f:
        for i,line in enumerate(f):
            if i == 0:
                continue
            data = line.split()
            particles.append(Particle(data[0],float(data[1]),float(data[2]),int(data[3]),int(data[4]),float(data[5]),data[6],int(data[7]),data[8]))
    return particles
def channels(channel,threshold):
    """
    Function that calculates the pairs (and trios) of particles that can form a given channel (set of quantum numbers) and a 
    given threshold energy
    :param channel: Dictionary containing the quantum numbers of the channel
    :param threshold: Threshold energy
    :return: Dictionary containing the pairs of particles that can form the channel and their threshold energies and a dictionary
    containing the trios of particles that can form the channel and their threshold energies
    """
    
    particles = read_particles('Particles/particles_unfl.txt')
    charmonmium = read_particles('Particles/charmonium.txt')
    names = []
    masses = []
    Channel_charm = channel['Charm']
    Channel_strange = channel['Strange']
    Channel_isospin = channel['Isospin']
    Channel_charm_isospin = channel['Charm_Isospin']
    Channel_C_parity = channel['C_parity']
    all_particles = charmonmium + particles
    for p1 in all_particles:
        for p2 in all_particles:
            if two_particles_inchannel(p1,p2,Channel_charm,Channel_strange,Channel_isospin,Channel_charm_isospin,Channel_C_parity):
                names.append((p1.name,p2.name))
                masses.append((p1.Mass+p2.Mass))
    seen = []
    real_m = []
    def modified_in(obj,list):
        for index,i in enumerate(list):
            if (obj[0] == i[0] and obj[1] == i[1]) or (obj[0] == i[1] and obj[1] == i[0]):
                return True,index
                
        return False,-1
    for index,i in enumerate(names):
        val, a = modified_in(i,seen)
        if not val:
            seen.append(i)
            real_m.append(masses[index])
    dic = {}
    for index,i in enumerate(seen):
        if real_m[index] < threshold:
            dic[i] = real_m[index]
    sorted_dic = dict(sorted(dic.items(), key=lambda item: item[1]))


    trios = []
    trios_mass = []
    all_particles = charmonmium + particles
    for cham in all_particles:
        for l in  all_particles:
            for m in all_particles:
                if three_particles_inchannel(cham,l,m,Channel_charm,Channel_strange,Channel_isospin,Channel_charm_isospin,Channel_C_parity):
                    trios.append((cham.name,l.name,m.name))
                    trios_mass.append(cham.Mass+l.Mass+m.Mass)
    seen2 = []
    real_m2 = []
    def modified_in2(obj,list):
        for index,i in enumerate(list):
            if (obj[0] == i[0] and obj[1] == i[1] and obj[2] == i[2]) or (obj[0] == i[1] and obj[1] == i[0] and obj[2] == i[2]) or (obj[0] == i[2] and obj[1] == i[1] and obj[2] == i[0]) or (obj[0] == i[1] and obj[1] == i[2] and obj[2] == i[0]) or (obj[0] == i[2] and obj[1] == i[0] and obj[2] == i[1]) or (obj[0] == i[0] and obj[1] == i[2] and obj[2] == i[1]):
                return True,index
                
        return False,-1
    for index,i in enumerate(trios):
        val, a = modified_in(i,seen2)
        if not val:
            seen2.append(i)
            real_m2.append(trios_mass[index])
    dic2 = {}
    for index,i in enumerate(seen2):
        if real_m2[index] < threshold:
            dic2[i] = real_m2[index]
    sorted_dic2 = dict(sorted(dic2.items(), key=lambda item: item[1]))
    return sorted_dic,sorted_dic2
"""
def all_channels():
   
    
    particles = read_particles('Particles/particles_unfl.txt')
    charmonmium = read_particles('Particles/charmonium.txt')
    Ds = read_particles('Particles/Ds.txt')
    names = []
    masses = []
    Channel_charm = 0
    Channel_strange = 0
    Channel_isospin = 0
    Channel_charm_isospin = 1
    Channel_C_parity = 1
    all_particles = charmonmium + particles + Ds
    for p1 in all_particles:
        for p2 in all_particles:
            if two_particles_inchannel(p1,p2,Channel_charm,Channel_strange,Channel_isospin,Channel_charm_isospin,Channel_C_parity):
                names.append((p1.name,p2.name))
                masses.append((p1.Mass+p2.Mass))
    seen = []
    real_m = []
    def modified_in(obj,list):
        for index,i in enumerate(list):
            if (obj[0] == i[0] and obj[1] == i[1]) or (obj[0] == i[1] and obj[1] == i[0]):
                return True,index
                
        return False,-1
    for index,i in enumerate(names):
        val, a = modified_in(i,seen)
        if not val:
            seen.append(i)
            real_m.append(masses[index])
    dic = {}
    for index,i in enumerate(seen):
        if real_m[index] < 0.73:
            dic[i] = real_m[index]
    sorted_dic = dict(sorted(dic.items(), key=lambda item: item[1]))


    trios = []
    trios_mass = []
    all_particles = charmonmium + particles + Ds
    for cham in all_particles:
        for l in  all_particles:
            for m in all_particles:
                if three_particles_inchannel(cham,l,m,Channel_charm,Channel_strange,Channel_isospin,Channel_charm_isospin,Channel_C_parity):
                    trios.append((cham.name,l.name,m.name))
                    trios_mass.append(cham.Mass+l.Mass+m.Mass)
    seen2 = []
    threshold = 0.73

    real_m2 = []
    def modified_in2(obj,list):
        for index,i in enumerate(list):
            if (obj[0] == i[0] and obj[1] == i[1] and obj[2] == i[2]) or (obj[0] == i[1] and obj[1] == i[0] and obj[2] == i[2]) or (obj[0] == i[2] and obj[1] == i[1] and obj[2] == i[0]) or (obj[0] == i[1] and obj[1] == i[2] and obj[2] == i[0]) or (obj[0] == i[2] and obj[1] == i[0] and obj[2] == i[1]) or (obj[0] == i[0] and obj[1] == i[2] and obj[2] == i[1]):
                return True,index
                
        return False,-1
    for index,i in enumerate(trios):
        val, a = modified_in(i,seen2)
        if not val:
            seen2.append(i)
            real_m2.append(trios_mass[index])
    dic2 = {}
    for index,i in enumerate(seen2):
        if real_m2[index] < threshold:
            dic2[i] = real_m2[index]
    sorted_dic2 = dict(sorted(dic2.items(), key=lambda item: item[1]))
    return sorted_dic,sorted_dic2
"""
def particle_dict():
    """
    Function that creates a dictionary with the particles and their names as keys, it loads all the files in the 
    Particles folder (by hand, could fix this in the future)
    :return: Dictionary with the particles and their names as keys
    """
    particles = read_particles('Particles/particles_unfl.txt')
    charmonmium = read_particles('Particles/charmonium.txt')
    Ds = read_particles('Particles/Ds.txt')
    all_particles = charmonmium + particles + Ds
    particle_dict = {}
    for p in all_particles:
        particle_dict[p.name] = p
    return particle_dict   
def possible_partial_waves_per_irrep_rest(channel,Jmax,threshold):
    """
    Function that calculates the possible partial waves for each irrep at rest given a channel (set of quantum numbers) and a threshold energy
    :param channel: Dictionary containing the quantum numbers of the channel
    :param Jmax: Maximum angular momentum being considered
    :param threshold: Threshold energy
    :return: Dictionary containing the possible partial waves for each irrep at rest, this is a dictionary of dictionaries.
    The keys of the first one being the irreps and the inner one keys being the pair of particles which subduce to that irrep, the final values
    are the partial waves of that pair of partciles which subduce to that irrep.
    """
    
    # Obtain dictionary of dictionaries with the first one having keys J^{PC} in a given channel and the inner one having pair of particles 
    # as keys and lists of partial waves as values
    partial_waves_dic = partial_waves(channel,Jmax,threshold)

    # Repeat for particle pairs with additional symmetry constraints
    Ds_partial_waves_dic = Ds_partial_waves(channel,Jmax,threshold)
    for i in Ds_partial_waves_dic.keys():
        if i in partial_waves_dic.keys():
            partial_waves_dic[i].update(Ds_partial_waves_dic[i])
        else:
            partial_waves_dic[i] = Ds_partial_waves_dic[i]
    # Create a dictionary with the possible angular momentum for each irrep
    irrep_Js = gt.J_in_irrep_rest()
    dic = {}
    # Iterate over the possible angular momentum for each irrep and create a dictionary with the possible partial waves for each irrep
    for key in irrep_Js.keys():
        for J in irrep_Js[key]:

            if J in partial_waves_dic.keys():
                for chan in partial_waves_dic[J].keys():
                    if key in dic.keys():
                        if chan in dic[key].keys():
                            dic[key][chan] += partial_waves_dic[J][chan].copy()
                        else:
                            dic[key][chan] = partial_waves_dic[J][chan].copy()
                        
                    else:
                        dic[key] = {chan:partial_waves_dic[J][chan].copy()}
                
    return dic
def possible_partial_waves_per_irrep_D4(channel,Jmax,threshold):
    """
    Function that calculates the possible partial waves for each irrep in a moving frame with one unit of momentum (100) given a channel
    (set of quantum numbers) and a threshold energy
    :param channel: Dictionary containing the quantum numbers of the channel
    :param Jmax: Maximum angular momentum being considered
    :param threshold: Threshold energy
    :return: Dictionary containing the possible partial waves for each irrep at rest, this is a dictionary of dictionaries.
    The keys of the first one being the irreps and the inner one keys being the pair of particles which subduce to that irrep, the final values
    are the partial waves of that pair of partciles which subduce to that irrep.
    """
    # Obtain dictionary of dictionaries with the first one having keys J^{PC} in a given channel and the inner one having pair of particles 
    # as keys and lists of partial waves as values
    partial_waves_dic = partial_waves(channel,Jmax,threshold)

    # Repeat for particle pairs with additional symmetry constraints
    Ds_partial_waves_dic = Ds_partial_waves(channel,Jmax,threshold)
    for i in Ds_partial_waves_dic.keys():
        if i in partial_waves_dic.keys():
            partial_waves_dic[i].update(Ds_partial_waves_dic[i])
        else:
            partial_waves_dic[i] = Ds_partial_waves_dic[i]
    # Create a dictionary with the possible angular momentum for each irrep with one unit of momenutm (100)
    irrep_Js = gt.Js_in_irrep_D4(Jmax)
    dic = {}

    # Iterate over the possible angular momentum for each irrep and create a dictionary with the possible partial waves for each irrep
    for key in irrep_Js.keys():
        
        for J in irrep_Js[key]:
           

            if J in partial_waves_dic.keys():
                for chan in partial_waves_dic[J].keys():
                    m = partial_waves_dic[J][chan]
                    if key in dic.keys():
                        if chan in dic[key].keys():
                            
                            for partial_wave in partial_waves_dic[J][chan]:
                               
                                    
                                
                                dic[key][chan] += [partial_wave]
                            
                               
                        else:
                            dic[key][chan] = partial_waves_dic[J][chan].copy()
                        
                    else:
                        dic[key] = {chan:partial_waves_dic[J][chan].copy()}
                
    return dic


def partial_waves(channel,Jmax,threshold):
    particles_dict = particle_dict()
    two_channel, _ = channels(channel,threshold)
    JP = {}

    for key in two_channel.keys():
        p1 = particles_dict[key[0]]
        p2 = particles_dict[key[1]]
        J1 = float(p1.J)
        J2 = float(p2.J)
        parity1 = p1.Parity
        parity2 = p2.Parity
    
        Jcombschannel = angular_mom.joinpartial_waves(J1,J2,Jmax,parity1,parity2)

        for kez in Jcombschannel.keys():
            key2 = '$' +p1.name+ " " + p2.name + '$'

            if kez in JP.keys():
                JP[kez].update({key2:format_1(kez,key2,Jcombschannel[kez])})
            else:
                JP[kez] = {key2:format_1(kez,key2,Jcombschannel[kez])}
    

    return JP
def create_table_partial_wave_irrep(channel,Jmax,threshold,mom):
    names = {0:'S',1:'P',2:'D',3:'F',4:'G',5:'H',6:'I',7:'J',8:'K',9:'L',10:'M',11:'N',12:'O'}
    max_l =4
    reversed_names = {v: k for k, v in names.items()}
    channelss = all_possible_channels()
    flag = False
    if channel['C_parity'] == 1:
        cp = '+'
    elif channel['C_parity'] == -1:
        cp = '-'
    if mom == "000":
        dic = possible_partial_waves_per_irrep_rest(channel,Jmax,threshold)
    elif mom == "100":
        dic = possible_partial_waves_per_irrep_D4(channel,Jmax,threshold)
        flag = True
    else:
        raise ValueError('Invalid momentum for now only 000 and 100 are supported')
    for key in dic.keys():
        irrep_name = key
        if not flag:
            ir = key.split('^')[0]
            par = key.split('^')[1]
            irrep_name = ir +'^{'+par+cp+'}'
        
        print( ' \\begin{table}[H]\\begin{tabularx}{\\textwidth}{YYY}  \\toprule Channel  &Partial Waves & $E_{th}$\\\\ \midrule')

        for ch in dic[key].keys():
           
            all_waves = ""
            for i,wave in enumerate(dic[key][ch]):
                letter = wave.split('_')[0][len(wave.split('_')[0])-1]
                number = reversed_names[letter]
                wavey = ""
                if number >= max_l:

                    wavey +="\\textcolor{gray}{"+ str(wave)+"}"
                else:
                    wavey += str(wave)

                if i == len(dic[key][ch])-1:
                    all_waves += wavey
                else: 
                    all_waves += wavey+ ','

            print(ch + " & " + all_waves + "&"+ str(round(channelss[ch],3))+'\\\\')
        print( '\\bottomrule\end{tabularx}')
        print('\caption{Partial Waves for $' + irrep_name +'$ with maximum threshold energy $E_t$ = $' +str(threshold)+'$} \end{table}')

def format_1(J,name,projection):
    texts = []
    names = {0:'S',1:'P',2:'D',3:'F',4:'G',5:'H',6:'I',7:'J',8:'K',9:'L',10:'M',11:'N',12:'O'}
    for proj in projection:
        S = proj[1]
        L = proj[0]
        J = J.split('.')[0]
        text = "$^"+str(int(2*S+1)) +names[L] + "_{"+str(J)+"}$"
        texts.append(text)
    return texts
def all_possible_channels():
    all_particles = read_particles('Particles/particles_unfl.txt') + read_particles('Particles/charmonium.txt') + read_particles('Particles/Ds.txt')
    channels = {}
    for p1 in all_particles:
        for p2 in all_particles:
            key =  '$' +p1.name+ " " + p2.name + '$'
            channels[key] = (p1.Mass+p2.Mass)
    return channels
def printing(channel,Jmax,threshold):
    print('\documentclass{article}\\usepackage{float}\\usepackage{graphicx}\\usepackage{xcolor}\\usepackage{amsmath}\\usepackage{tabularx,booktabs} \\newcolumntype{Y}{>{\\centering\\arraybackslash}X}\\title{channels}\\author{Pablo Morande}\\date{November 2024}\\begin{document}')
    JP = partial_waves(channel,Jmax,threshold)
    #JP = Ds_partial_waves(channel,Jmax,threshold)
    JP2 = Ds_partial_waves(channel,Jmax,threshold)
    #JP2 = {}
    channelss = all_possible_channels()
    C_p_dic = {1:'+',-1:'-',0:''}
    C_p = C_p_dic[channel['C_parity']]
    max_l = 4
    names = {0:'S',1:'P',2:'D',3:'F',4:'G',5:'H',6:'I',7:'J',8:'K',9:'L',10:'M',11:'N',12:'O'}

    reversed_names = {v: k for k, v in names.items()}
    ch = channels(channel,threshold)[0]
    ch2 = Ds(channel,threshold)
    for i in ch2.keys():
        ch[i] = ch2[i]
    print( ' \\begin{table}[H]\\begin{tabularx}{\\textwidth}{YY}  \\toprule Channel  & Threshold E\\\\ \midrule')

    for key in ch.keys():
        name ="$"+ key[0]+' '+key[1]+"$"
        print(name + " & " + str(round(ch[key],2)) +'\\\\')
    print( '\\bottomrule\end{tabularx}')
    print('\caption{Channels for Max threshold $E_t$ = $' +str(threshold)+'$} \end{table}')
    print('\n')

        
        


    for i in JP2.keys():
        JP[i].update(JP2[i]) 
    for key in JP.keys():
        print( ' \\begin{table}[H]\\begin{tabularx}{\\textwidth}{YYY}  \\toprule Channel &Partial Waves & Threshold E\\\\ \midrule')
        for key2 in JP[key].keys():
            string = ''
            for index,a in enumerate(JP[key][key2]):
                if index == len(JP[key][key2])-1:
                    letter = a.split('_')[0][len(a.split('_')[0])-1]
                    number = reversed_names[letter]
                    if number >= max_l:

                        string +="\\textcolor{gray}{"+ str(a)+"}"
                    else:
                        string += str(a)
                else:
                    letter = a.split('_')[0][len(a.split('_')[0])-1]
                    number = reversed_names[letter]
                    if number >= max_l:

                        string +="\\textcolor{gray}{"+ str(a)+"}" +','
                    else:
                        string += str(a) +','
            print(key2 + " & " + string + "&"+ str(round(channelss[key2],3))+'\\\\')
        print( '\\bottomrule\end{tabularx}')
        print('\caption{Partial Waves for $' + key.strip('.')[0]+'^{'+key.strip('.')[3] +C_p+ '}$ with maximum threshold energy $E_t$ = $' +str(threshold)+'$} \end{table}')
        print('\n')
    print('\end{document}')
def bar(string):
    if "\\bar{" in string:
        b = string.split('\\bar{')[1]
        l = ''
        count = 0
        for let in b:
            if let == '{':
                count += 1
            if let == '}' and count == 0:
                continue
            if let == '}' and count != 0:
                count -= 1
           
            l += let

        return l
    else:
        b = ''
        count = 0
        for let in string:
            if (let == '_' or let == '^' or let == '*') and count == 0:
                b+= '}'
                count += 1
            b += let
        if count == 0:
            b += '}'
        return '\\bar{'+b 
def Ds(channel,threshold):
    names = []
    masses = []
    Channel_charm = channel['Charm']
    Channel_strange = channel['Strange']
    Channel_isospin = channel['Isospin']
    Channel_charm_isospin = channel['Charm_Isospin']
    Channel_C_parity = channel['C_parity']
    all_particles =   read_particles('Particles/Ds.txt')
    for p1 in all_particles:
        for p2 in all_particles:
            if two_particles_inchannel(p1,p2,Channel_charm,Channel_strange,Channel_isospin,Channel_charm_isospin,Channel_C_parity):
                names.append((p1.name,p2.name))
                masses.append((p1.Mass+p2.Mass))
    seen = []
    real_m = []
    def modified_in(obj,list):
        for index,i in enumerate(list):
            if (obj[0] == i[0] and obj[1] == i[1]) or (obj[0] == i[1] and obj[1] == i[0]):
                return True,index
                
        return False,-1
    for index,i in enumerate(names):
        val, a = modified_in(i,seen)
        if not val:
            seen.append(i)
            real_m.append(masses[index])
    dic = {}
    for index,i in enumerate(seen):
        if real_m[index] < threshold:
            dic[i] = real_m[index]
    seen_keys = []
    final_dic = {}
    for index,i in enumerate(dic.keys()):
        if i not in seen_keys:
            seen_keys.append(i)
            alternative_key = (bar(i[1]),bar(i[0]))
            alternative_key_2 = (bar(i[0]),bar(i[1]))
            seen_keys.append(alternative_key)
            seen_keys.append(alternative_key_2)
            final_dic[i] = dic[i]
        
    sorted_dic = dict(sorted(final_dic.items(), key=lambda item: item[1]))
    return sorted_dic
def Ds_partial_waves(channel,Jmax,threshold):
    particles_dict = particle_dict()
    two_channel = Ds(channel,threshold)
    JP = {}
    for key in two_channel.keys():
        p1 = particles_dict[key[0]]
        p2 = particles_dict[key[1]]
        J1 = float(p1.J)
        J2 = float(p2.J)
        parity1 = p1.Parity
        parity2 = p2.Parity
        name1 = p1.name
        name2 = bar(p2.name)
        C_parity = channel['C_parity']
        Jcombschannel = angular_mom.joinpartial_waves_Ds(J1,J2,Jmax,parity1,parity2,C_parity,name1,name2)
        for kez in Jcombschannel.keys():
            key2 = '$' +p1.name+ " " + p2.name + '$'

            if kez in JP.keys():
                if kez == '4.0+':
                    pass
                JP[kez].update({key2:format_1(kez,key2,Jcombschannel[kez])})
            else:
                JP[kez] = {key2:format_1(kez,key2,Jcombschannel[kez])}
    

    return JP




channel = {'Charm':0,'Strange':0,'Isospin':0,'Charm_Isospin':1,'C_parity':1}
threshold = 0.74


#(Ds_partial_waves(channel,4))
#print(channels(channel,threshold))
#printing(channel,4,threshold)
#print(angular_mom.possible_Js_partial_waves(1,1,3)[1][2][1])
#print(bar('D'))
def possible_irreps_rest(p_1,p_2,channel,Jmax,threshold):
    key = '$'+p_1+' '+p_2+'$'
    alt_key = '$'+p_2+' '+ p_1+'$'
    JP = partial_waves(channel,Jmax,threshold)
    #JP = Ds_partial_waves(channel,Jmax,threshold)
    JP2 = Ds_partial_waves(channel,Jmax,threshold)
    #JP2 = {}
    channelss = all_possible_channels()
    C_p_dic = {1:'+',-1:'-',0:''}
    C_p = C_p_dic[channel['C_parity']]
    max_l = 4
    names = {0:'S',1:'P',2:'D',3:'F',4:'G',5:'H',6:'I',7:'J',8:'K',9:'L',10:'M',11:'N',12:'O'}

    reversed_names = {v: k for k, v in names.items()}
    for i in JP2.keys():
        if i in JP.keys():
            JP[i].update(JP2[i])
        else:
            JP[i] = JP2[i]
    possible_angular_momentum = []
    m_l = 3
    for keys in JP.keys():
        for key2 in JP[keys].keys():
            if key2 == key or alt_key == key2:
                if key2 == key:
                    for i in JP[keys][key2]:
                        stri = i.split('_')[0][len(i.split('_')[0])-1]
                        num = reversed_names[stri]
                        if num <= m_l:
                            possible_angular_momentum.append(keys)
                            break
                else:
                    for i in JP[keys][alt_key]:
                        stri = i.split('_')[0][len(i.split('_')[0])-1]
                        num = reversed_names[stri]
                        if num <= m_l:
                            possible_angular_momentum.append(keys)
                            break
                break
    possible_irreps = []
    dict = {}
    dict['0.0+'] = ['A_1^+']
    dict['0.0-'] = ['A_1^-']
    dict['4.0+'] = ['A_1^+', 'T_1^+', 'T_2^+','E^+']
    dict['4.0-'] = ['A_1^-', 'T_1^-', 'T_2^-','E^-']
    dict['3.0+'] = ['A_2^+','T_1^+', 'T_2^+']
    dict['3.0-']= ['A_2^-','T_1^-', 'T_2^-']
    dict['1.0+'] = ['T_1^+']
    dict['1.0-'] = ['T_1^-']
    dict['2.0+'] = ['T_2^+', 'E^+']
    dict['2.0-'] = ['T_2^-', 'E^-']
    for i in possible_angular_momentum:
        possible_irreps += dict[i]
    return possible_irreps

def all_particles_table():
    all_particles = read_particles('Particles/particles_unfl.txt') + read_particles('Particles/charmonium.txt') + read_particles('Particles/Ds.txt')
    print( ' \\begin{table}[H]\\begin{tabularx}{\\textwidth}{YYY}  \\toprule Particle  & $J^{P(C)}$& $a_t m$\\\\ \midrule')
    for p in all_particles:
        symbol = '+' if p.Parity == 1 else '-'
        if p.C_parity == 1:
            C_parity = '+'
        elif p.C_parity == -1:
            C_parity = '-'
        else:
            C_parity = ''
        print('$'+p.name +'$'+ " & " +'$'+ str(p.J) +'^{'+str(symbol)+C_parity +'}' +'$'+'& $'+ str(p.Mass)+'$\\\\')
    print( '\\bottomrule\end{tabularx}\end{table}')
#all_particles_table()


def check_partial_waves_in_irrep(channel,Jmax,threshold):
    dic_irreps = possible_partial_waves_per_irrep_rest(channel,Jmax,threshold)
    particles_dic = particle_dict()
    names = {0:'S',1:'P',2:'D',3:'F',4:'G',5:'H',6:'I',7:'J',8:'K',9:'L',10:'M',11:'N',12:'O'}
    possible_Jp_in_irrep = gt.J_in_irrep_rest()

    reversed_names = {v: k for k, v in names.items()}
    for key in dic_irreps.keys():
        for chan in dic_irreps[key].keys():
            for wave in dic_irreps[key][chan]:
                l = reversed_names[wave[3]]
                J = float(wave[6])
                def remove_dol(str):
                    r = ''
                    for s in str:
                        if s != '$':
                            r += s
                    return r
                p_1 = remove_dol(chan.split(' ')[0])
                p_2 = remove_dol(chan.split(' ')[1])
                for particle in particle_dict().keys():
                    if particle == p_1:
                        p1 = particle_dict()[particle]
                    if particle == p_2:
                        p2 = particle_dict()[particle]
                parity_1 = p1.Parity
                parity_2 = p2.Parity
                parity_t = parity_1*parity_2*(-1)**(l)
                if parity_t == 1:
                    parity = '+'
                else:
                    parity = '-'
                J_p = str(J) + parity
                if J_p in possible_Jp_in_irrep[key]:
                    print('The partial wave '+ wave + ' is in the irrep ' + key + ' for the channel ' + p_1 + ' ' + p_2)
                else:
                    raise ValueError('The partial wave '+ wave + ' is not in the irrep ' + key + ' for the channel ' + p_1 + ' ' + p_2)
                

def check_partial_waves_in_irrep_D4(channel,Jmax,threshold):
    dic_irreps = possible_partial_waves_per_irrep_D4(channel,Jmax,threshold)
    particles_dic = particle_dict()
    names = {0:'S',1:'P',2:'D',3:'F',4:'G',5:'H',6:'I',7:'J',8:'K',9:'L',10:'M',11:'N',12:'O'}
    possible_Jp_in_irrep = gt.Js_in_irrep_D4(Jmax)

    reversed_names = {v: k for k, v in names.items()}
    for key in dic_irreps.keys():
        for chan in dic_irreps[key].keys():
            for wave in dic_irreps[key][chan]:
                l = reversed_names[wave[3]]
                J = float(wave[6])
                def remove_dol(str):
                    r = ''
                    for s in str:
                        if s != '$':
                            r += s
                    return r
                p_1 = remove_dol(chan.split(' ')[0])
                p_2 = remove_dol(chan.split(' ')[1])
                for particle in particle_dict().keys():
                    if particle == p_1:
                        p1 = particle_dict()[particle]
                    if particle == p_2:
                        p2 = particle_dict()[particle]
                parity_1 = p1.Parity
                parity_2 = p2.Parity
                parity_t = parity_1*parity_2*(-1)**(l)
                if parity_t == 1:
                    parity = '+'
                else:
                    parity = '-'
                J_p = str(J) + parity
                if J_p in possible_Jp_in_irrep[key]:
                    print('The partial wave '+ wave + ' is in the irrep ' + key + ' for the channel ' + p_1 + ' ' + p_2)
                else:
                    raise ValueError('The partial wave '+ wave + ' is not in the irrep ' + key + ' for the channel ' + p_1 + ' ' + p_2)
                


                
           