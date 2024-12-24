import quantum_numbers
import angular_mom
import matplotlib.pyplot as plt
class Particle(object):
    def __init__(self,name,Isospin,Isospin_charm, C_parity,Charm,Strange,Mass,Parity,J):
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
        Parity = '+' if self.Parity == 1 else '-'
        if self.C_parity == 1:
            C_parity = '+'
        elif self.C_parity == -1:
            C_parity = '-'
        else:
            C_parity = ''
        return self.name + ' ' + str(self.Isospin) + ' ' +  str(self.J) + '^{'+str(Parity)+str(C_parity) + '} '
def two_particles_inchannel(particel_one,particle_two,channel_charm,channel_strange,channel_isospin,channel_charm_isospin,channel_C_parity):
    Isospins = quantum_numbers.Isospin(particel_one.Isospin,particle_two.Isospin)
    charm_Isospins = quantum_numbers.Isospin(particel_one.Isospin_charm,particle_two.Isospin_charm)
    Charm = quantum_numbers.Charm(particel_one.Charm,particle_two.Charm)
    Strange = quantum_numbers.Strangeness(particel_one.Strange,particle_two.Strange)
    if particel_one.Charm == 0 and particle_two.Charm == 0 and particel_one.Strange == 0 and particle_two.Strange == 0:
        C_parity =quantum_numbers.Charge_conjugation(particel_one.C_parity,particle_two.C_parity)

        if channel_charm== Charm and channel_strange == Strange and channel_isospin in Isospins and channel_charm_isospin in charm_Isospins and channel_C_parity == C_parity:
            return True
        else:
            return False
    else:
        if channel_charm== Charm and channel_strange == Strange and channel_isospin in Isospins and channel_charm_isospin in charm_Isospins:
            return True
        else:
            return False
def three_particles_inchannel(particel_one,particle_two,particle_three,channel_charm,channel_strange,channel_isospin,channel_charm_isospin,channel_C_parity):
    Isospins = quantum_numbers.Isospin_three_particles(particel_one.Isospin,particle_two.Isospin,particle_three.Isospin)
    charm_Isospins = quantum_numbers.Isospin_three_particles(particel_one.Isospin_charm,particle_two.Isospin_charm,particle_three.Isospin_charm)
    Charm = quantum_numbers.Charm_three_particles(particel_one.Charm,particle_two.Charm,particle_three.Charm)
    Strange = quantum_numbers.Strangeness_three_particles(particel_one.Strange,particle_two.Strange,particle_three.Strange)
    C_parity = quantum_numbers.Charge_conjugation_three_particles(particel_one.C_parity,particle_two.C_parity,particle_three.C_parity)
    if particel_one.Charm == 0 and particle_two.Charm == 0 and particle_three.Charm == 0 and particel_one.Strange == 0 and particle_two.Strange == 0 and particle_three.Strange == 0:
        if channel_charm== Charm and channel_strange == Strange and channel_isospin in Isospins and channel_charm_isospin in charm_Isospins and channel_C_parity == C_parity:
            return True
        else:
            return False
    else:
        if channel_charm== Charm and channel_strange == Strange and channel_isospin in Isospins and channel_charm_isospin in charm_Isospins:
            return True
        else:
            return False

def read_particles(filename):
    particles = []
    with open(filename) as f:
        for i,line in enumerate(f):
            if i == 0:
                continue
            data = line.split()
            particles.append(Particle(data[0],float(data[1]),float(data[2]),int(data[3]),int(data[4]),float(data[5]),data[6],int(data[7]),data[8]))
    return particles
def channels(channel,threshold):
    
    particles = read_particles('particles_unfl.txt')
    charmonmium = read_particles('charmonium.txt')
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
def all_channels():
    
    particles = read_particles('particles_unfl.txt')
    charmonmium = read_particles('charmonium.txt')
    Ds = read_particles('Ds.txt')
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
def particle_dict():
    particles = read_particles('particles_unfl.txt')
    charmonmium = read_particles('charmonium.txt')
    Ds = read_particles('Ds.txt')
    all_particles = charmonmium + particles + Ds
    particle_dict = {}
    for p in all_particles:
        particle_dict[p.name] = p
    return particle_dict   
                 
def partial_waves(channel,Jmax,threshold):
    particles_dict = particle_dict()
    two_channel, three_channel = channels(channel,threshold)
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
    all_particles = read_particles('particles_unfl.txt') + read_particles('charmonium.txt') + read_particles('Ds.txt')
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
    if string[0] == '\\':
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
    all_particles =   read_particles('Ds.txt')
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
    all_particles = read_particles('particles_unfl.txt') + read_particles('charmonium.txt') + read_particles('Ds.txt')
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


