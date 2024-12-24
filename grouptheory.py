import particle as p
import numpy as np
def identify_irreps_particle_rest(name):
    all_particles = p.read_particles('Particles/particles_unfl.txt') + p.read_particles('Particles/charmonium.txt')+ p.read_particles('Particles/Ds.txt') 
    for particle in all_particles:
        if particle.name == name:
            J = particle.J
            P = '+' if particle.Parity == 1 else '-'
    filename = 'GroupTheory/O_D^h.txt'
    irreps = []
    with open(filename) as f:
        lines = f.readlines()
        for line in lines:
            data = line.split()
            if data[0][len(data[0])-1] != P:
                continue
            for i in range(1,len(data)):
                if data[i] == J:
                    irreps.append(data[0])
                    break
    return irreps


def idetintify_irreps_two_particles_rest(name1,name2):
    irreps1 = identify_irreps_particle_rest(name1)
    irreps2 = identify_irreps_particle_rest(name2)
    combinations = []
    for irrep1 in irreps1:
        for irrep2 in irreps2:
            combinations.append((irrep1,irrep2))
    return combinations
def subductions():
    dict = {'A_1':0,'A_2':1,'E':2,'T_1':3,'T_2':4}
    reversed_dict = {0:'A_1',1:'A_2',2:'E',3:'T_1',4:'T_2'}
    table = [[0 for i in range(5)] for j in range(10)]
    table[0][0] = [0]
    table[0][1] = [1]
    table[0][2] = [2]
    table[0][3] = [3]
    table[0][4] = [4]
    table[1][0] = table[0][1]
    table[1][1] = [0]
    table[1][2] = [2]
    table[1][3] = [4]
    table[1][4] = [3]
    table[2][0] = table[0][2]
    table[2][1] = table[1][2]
    table[2][2] = [0,1,2]
    table[2][3] = [3,4]
    table[2][4] = [3,4]
    table[3][0] = table[0][3]
    table[3][1] = table[1][3]
    table[3][2] = table[2][3]
    table[3][3] = [0,2,3,4]
    table[3][4] = [1,2,3,4]
    table[4][0] = table[0][4]
    table[4][1] = table[1][4]
    table[4][2] = table[2][4]
    table[4][3] = table[3][4]
    table[4][4] = [0,2,3,4]
    real_table = [[[reversed_dict[table[i][j][z]] for z in range(len(table[i][j]))] for i in range(5)] for j in range(5)]
    return real_table
#print(idetintify_irreps_two_particles_rest('\eta_c','\eta_c'))
def subductions_rest(irrep1,irrep2):
    subduction = subductions()
    dict = {'A_1':0,'A_2':1,'E':2,'T_1':3,'T_2':4}

    reversed_dict = {0:'A_1',1:'A_2',2:'E',3:'T_1',4:'T_2'}
    index1 = dict[irrep1]
    index2 = dict[irrep2]
    return subduction[index1][index2]


def identify_ni_at_rest_levels_subduction(name1,name2):
    combinations = idetintify_irreps_two_particles_rest(name1,name2)
    irreps = []
    for combination in combinations:
        parity_1 = combination[0][len(combination[0])-1]
        parity_2 = combination[1][len(combination[1])-1]
        name_irrep1 = combination[0][0:len(combination[0])-2]
        name_irrep2 = combination[1][0:len(combination[1])-2]
        if parity_1 != parity_2:
            P = '-'
        else:
            P = '+'
        ls =  subductions_rest(name_irrep1,name_irrep2)
        for l in range(len(ls)):
            ls[l]+= '^'+P
        irreps += ls
        
        #irreps.append(subductions_rest(combination[0],combination[1]))
    return irreps

#print(identify_ni_at_rest_levels_subduction('\psi','w'))
def identify_irreps_particle_Dic4(name):
    all_particles = p.read_particles('Particles/particles_unfl.txt') + p.read_particles('Particles/charmonium.txt')+ p.read_particles('Particles/Ds.txt') 
    for particle in all_particles:
        if particle.name == name:
            J = particle.J
            P = '+' if particle.Parity == 1 else '-'
            break
    filename = 'GroupTheory/Dic_4.txt'
    irreps = []
    eta = '+' if particle.Parity*(-1)**int(J) == 1 else '-'
    possible_lambdas = np.arange(0,int(J)+1)
    possible_lambdas = [str(i)+eta if i ==0 else str(i) for i in possible_lambdas]
    with open(filename) as f:
        lines = f.readlines()
        for line in lines:
            data = line.split()
            for i in range(1,len(data)):
                if data[i] in possible_lambdas :
                    irreps.append(data[0])
                    break
    return irreps

def idetintify_irreps_two_particles_Dic4(name1,name2):
    irreps1 = identify_irreps_particle_Dic4(name1)
    irreps2 = identify_irreps_particle_Dic4(name2)
    combinations = []
    for irrep1 in irreps1:
        for irrep2 in irreps2:
            combinations.append((irrep1,irrep2))
    return combinations

#print(idetintify_irreps_two_particles_Dic4('\psi','w'))

def subductions_Dic4():
    dict1 = {'A_1':0,'A_2':1,'B_1':2,'B_2':3,'E_2':4}
    dict2 = {'A_1^+':0,'A_2^+':1,'E^+':2,'T_1^+':3,'T_2^+':4,'A_1^-':5,'A_2^-':6,'E^-':7,'T_1^-':8,'T_2^-':9}
    reversed_dict1 = {0:'A_1',1:'A_2',2:'B_1',3:'B_2',4:'E_2'}
    reversed_dict2 = {0:'A_1^+',1:'A_2^+',2:'E^+',3:'T_1^+',4:'T_2^+',5:'A_1^-',6:'A_2^-',7:'E^-',8:'T_1^-',9:'T_2^-'}
    table = [[0 for i in range(5)] for j in range(5)]
    table[0][0] = ['A_1^+','E^+','T_1^-']
    table[0][1] = ['T1^+','E^-','A_1^-']
    table[0][2] = ['A_2^+','T_2^-','E^+']
    table[0][3] =['T2^+','E^-','A_2^-']
    table[0][4] = ['T_1^-','T_2^-','T_1^+','T_2^+']
    table[1][0] = table[0][1]
    table[1][1] = ['A_1^+','E^+','T_1^-']
    table[1][2] = ['T_2^+','E^-','A_2^-']
    table[1][3] = ['A_2^+','T_2^-','E^+']
    table[1][4] = ['T_1^-','T_2^-','T_1^+','T_2^+']
    table[2][0] = table[0][2]
    table[2][1] = table[1][2]
    table[2][2] = ['A_1^+','E^+','T_1^-']
    table[2][3] = ['T_1^+','E^-','A_1^-']
    table[2][4] = ['T_1^-','T_2^-','T_1^+','T_2^+']
    table[3][0] = table[0][3]
    table[3][1] = table[1][3]
    table[3][2] = table[2][3]
    table[3][3] = ['A_1^+','T_1^-','E^+']
    table[3][4] = ['T_1^-','T_2^-','T_1^+','T_2^+']
    table[4][0] = table[0][4]
    table[4][1] = table[1][4]
    table[4][2] = table[2][4]
    table[4][3] = table[3][4]
    table[4][4] = ['A_1^+','A_2^+','E^+','E^+','T_1^+','T_2^+','A_1^-','A_2^-','E^-','E^-','T_1^-','T_2^-']

    return table
def subductions_Dic4_from_irreps(irrep1,irrep2):
    subduction = subductions_Dic4()
    dict1 = {'A_1':0,'A_2':1,'B_1':2,'B_2':3,'E_2':4}
    reversed_dict1 = {0:'A_1',1:'A_2',2:'B_1',3:'B_2',4:'E_2'}
    index1 = dict1[irrep1]
    index2 = dict1[irrep2]

    return subduction[index1][index2]
#print(subductions_Dic4_from_irreps('E_2','E_2'))
def identify_ni_Dic4_levels_subduction(name1,name2):
    combinations = idetintify_irreps_two_particles_Dic4(name1,name2)
    irreps = []
    for combination in combinations:
        name_irrep1 = combination[0]
        name_irrep2 = combination[1]

        irreps += subductions_Dic4_from_irreps(name_irrep1,name_irrep2)
        #irreps.append(subductions_rest(combination[0],combination[1]))
    dict = {}
    for irrep in irreps:
        if irrep in dict.keys():
            dict[irrep] += 1
        else:
            dict[irrep] = 1
    irrepss = []
    for key in dict.keys():
        
        irrepss.append(str(dict[key])+key)
    return irrepss
def identify_irreps_particle_norest(name,group):
    all_particles = p.read_particles('Particles/particles_unfl.txt') + p.read_particles('Particles/charmonium.txt')+ p.read_particles('Particles/Ds.txt') 
    for particle in all_particles:
        if particle.name == name:
            J = particle.J
            P = '+' if particle.Parity == 1 else '-'
    filename = group
    irreps = []
    eta = '+' if particle.Parity*(-1)**int(J) == 1 else '-'
    possible_lambdas = np.arange(0,int(J)+1)
    possible_lambdas = [str(i)+eta if i ==0 else str(i) for i in possible_lambdas]
    with open(filename) as f:
        lines = f.readlines()
        for line in lines:
            data = line.split()
            for i in range(1,len(data)):
                if data[i] in possible_lambdas :
                    irreps.append(data[0])
                    break
    return irreps
def idetintify_irreps_two_particles_norest(name1,name2,group):
    irreps1 = identify_irreps_particle_norest(name1,group)
    irreps2 = identify_irreps_particle_norest(name2,group)
    combinations = []
    for irrep1 in irreps1:
        for irrep2 in irreps2:
            combinations.append((irrep1,irrep2))
    return combinations
#print(identify_irreps_particle_norest('w','Dic_2.txt'))
def subductions_Dic2():
    dict1 = {'A_1':0,'A_2':1,'B_1':2,'B_2':3,}
    dict2 = {'A_1^+':0,'A_2^+':1,'E^+':2,'T_1^+':3,'T_2^+':4,'A_1^-':5,'A_2^-':6,'E^-':7,'T_1^-':8,'T_2^-':9}
    reversed_dict1 = {0:'A_1',1:'A_2',2:'B_1',3:'B_2',4:'E_2'}
    reversed_dict2 = {0:'A_1^+',1:'A_2^+',2:'E^+',3:'T_1^+',4:'T_2^+',5:'A_1^-',6:'A_2^-',7:'E^-',8:'T_1^-',9:'T_2^-'}
    table = [[0 for i in range(4)] for j in range(4)]
    table[0][0] = ['A_1^+','E^+','T_2^+','T_1^-','T_2^-']
    table[0][1] = ['A_1^-','E^-','T_2^-','T_1^+','T_2^+']
    table[0][2] = ['A_2^+','E^+','T_1^+','T_1^-','T_2^-']
    table[0][3] =['A_2^-','E^-','T_2^+','T_1^+','T_1^-']
    table[1][0] = table[0][1]
    table[1][1] = ['A_1^+','E^+','T_2^+','T_1^-','T_2^-']
    table[1][2] = ['A_2^-','E^-','T_2^+','T_1^-','T_1^+']
    table[1][3] = ['A_2^+','E^+','T_1^-','T_1^+','T_2^-']
    table[2][0] = table[0][2]
    table[2][1] = table[1][2]
    table[2][2] = ['A_1^+','E^+','T_2^+','T_1^-','T_2^-']
    table[2][3] = ['A_1^-','E^-','T_2^-','T_1^+','T_2^+']
    table[3][0] = table[0][3]
    table[3][1] = table[1][3]
    table[3][2] = table[2][3]
    table[3][3] = ['A_1^+','E^+','T_2^+','T_1^-','T_2^-']
   
    return table
def subductions_Dic2_from_irreps(irrep1,irrep2):
    subduction = subductions_Dic2()
    dict1 = {'A_1':0,'A_2':1,'B_1':2,'B_2':3}
    reversed_dict1 = {0:'A_1',1:'A_2',2:'B_1',3:'B_2',4:'E_2'}
    index1 = dict1[irrep1]
    index2 = dict1[irrep2]

    return subduction[index1][index2]
#print(subductions_Dic4_from_irreps('E_2','E_2'))
def identify_ni_Dic2_levels_subduction(name1,name2):
    combinations = idetintify_irreps_two_particles_norest(name1,name2,'GroupTheory/Dic_2.txt')
    irreps = []
    for combination in combinations:
        name_irrep1 = combination[0]
        name_irrep2 = combination[1]

        irreps += subductions_Dic2_from_irreps(name_irrep1,name_irrep2)
        #irreps.append(subductions_rest(combination[0],combination[1]))
    dict = {}
    for irrep in irreps:
        if irrep in dict.keys():
            dict[irrep] += 1
        else:
            dict[irrep] = 1
    irrepss = []
    for key in dict.keys():
        
        irrepss.append(str(dict[key])+key)
    return irrepss
#print(identify_ni_Dic4_levels_subduction('\psi','w'))
def subductions_Dic3():
    dict1 = {'A_1':0,'A_2':1,'E_2':2}
    dict2 = {'A_1^+':0,'A_2^+':1,'E^+':2,'T_1^+':3,'T_2^+':4,'A_1^-':5,'A_2^-':6,'E^-':7,'T_1^-':8,'T_2^-':9}
    reversed_dict1 = {0:'A_1',1:'A_2',2:'E_2'}
    reversed_dict2 = {0:'A_1^+',1:'A_2^+',2:'E^+',3:'T_1^+',4:'T_2^+',5:'A_1^-',6:'A_2^-',7:'E^-',8:'T_1^-',9:'T_2^-'}
    table = [[0 for i in range(4)] for j in range(4)]
    table[0][0] = ['A_1^+','T_2^+','A_2^-','T_1^-']
    table[0][1] = ['A_2^+','T_1^+','A_1^-','T_2^-']
    table[0][2] = ['T_2^+','E^+','T_1^+','E^-','T_1^-','T_2^-']
    table[1][0] = table[0][1]
    table[1][1] = ['A_1^+','T_2^+','A_2^-','T_1^-']
    table[1][2] = ['T_2^+','E^+','T_1^+','E^-','T_1^-','T_2^-']
    table[2][0] = table[0][2]
    table[2][1] = table[1][2]
    table[2][2] = ['A_1^+','A_2^+','E^+','T_1^+','T_1^+','T_2^+','T_2^+','A_1^-','A_2^-','E^-','T_1^-','T_2^-','T_1^-','T_2^-']
   
    return table
def subductions_Dic3_from_irreps(irrep1,irrep2):
    subduction = subductions_Dic3()
    dict1 = {'A_1':0,'A_2':1,'E_2':2}
    reversed_dict1 = {0:'A_1',1:'A_2',2:'B_1',3:'B_2',4:'E_2'}
    index1 = dict1[irrep1]
    index2 = dict1[irrep2]

    return subduction[index1][index2]
def identify_ni_Dic3_levels_subduction(name1,name2):
    combinations = idetintify_irreps_two_particles_norest(name1,name2,'GroupTheory/Dic_3.txt')
    irreps = []
    for combination in combinations:
        name_irrep1 = combination[0]
        name_irrep2 = combination[1]

        irreps += subductions_Dic3_from_irreps(name_irrep1,name_irrep2)
        #irreps.append(subductions_rest(combination[0],combination[1]))
    dict = {}
    for irrep in irreps:
        if irrep in dict.keys():
            dict[irrep] += 1
        else:
            dict[irrep] = 1
    irrepss = []
    for key in dict.keys():
        
        irrepss.append(str(dict[key])+key)
    return irrepss
#print(identify_ni_at_rest_levels_subduction('\psi','\phi'))