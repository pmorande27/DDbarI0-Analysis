import numpy as np
def possible_s_twoParticle(s1,s2):
    return np.arange(abs(s1-s2),s1+s2+1)

def angular_momentum(l,s):
    return np.arange(abs(l-s),l+s+1)
def possible_Js(s,J_max):
    J_max = J_max+1
    J = {}
    Jcombs = {}
    J[J_max] = 0
    l = 0
    while True:

        possible_Js = angular_momentum(l,s)
        for j in possible_Js:
            if j in Jcombs.keys():
                Jcombs[j] += [(l,s)]
            else:
                Jcombs[j] = [(l,s)]

            if j in J.keys():
                J[j] += 1
                
            else:
                J[j] = 1
        l+= 1
        if J[J_max] == 2*s + 1:
            break
    sorted_J = dict(sorted(J.items()))
    Jcombs = dict(sorted(Jcombs.items()))

    return sorted_J,Jcombs
def possible_Js_partial_waves(s1,s2,J_max):
    ss = possible_s_twoParticle(s1,s2)
    Js = []
    Jcombs =[]
    for s in ss:
        J,Jcomb = possible_Js(s,J_max)
        Js.append(J)
        Jcombs.append(Jcomb)
    return Js,Jcombs
def joinpartial_waves(s1,s2,J_max,p1,p2):
    Js,Jcombs = possible_Js_partial_waves(s1,s2,J_max)
    total_Jcombs = {}
    for index,j in enumerate(Jcombs):
        for key in j.keys():
            if key > J_max:
                continue
                
            if key in total_Jcombs.keys():
                total_Jcombs[key] += Jcombs[index][key]
            else:
                total_Jcombs[key] = Jcombs[index][key]
    parity_factor = p1*p2
    Jparity = {}
    for key in total_Jcombs.keys():
        combs = total_Jcombs[key]
        for index,comb in enumerate(combs):
            pkey = str(key)+ str("+")
            mkey = str(key)+ str("-")
            comb_p = parity_factor*(-1)**comb[0]
            
            if comb_p >0:
                if pkey in Jparity.keys():
                    Jparity[pkey] += [comb]
                else:
                    Jparity[pkey] = [comb]
            else:
                if mkey in Jparity.keys():
                    Jparity[mkey] += [comb]
                else:
                    Jparity[mkey] = [comb]
    return Jparity

def possible_Js_Ds(s,J_max,C_parity,name1,name2):

    J_max = J_max+1
    J = {}
    Jcombs = {}
    J[J_max] = 0
    l = 0
    
    if name1 == name2:
        f = 1
    else:
        f = -1
    while True:

        possible_Js = angular_momentum(l,s)
        for j in possible_Js:
            if j in Jcombs.keys():

                if name1 == name2:
                    if (-1)**(l+s) == C_parity:
                        Jcombs[j] += [(l,s)]
                else:
                    Jcombs[j] += [(l,s)]
            else:
                if name1 == name2:

                    if (-1)**(l+s) == C_parity:

                        Jcombs[j] = [(l,s)]
                else:
                    Jcombs[j] = [(l,s)]

            if j in J.keys():
                J[j] += 1
                
            else:
                J[j] = 1
        l+= 1
        if J[J_max] == 2*s + 1:
            break
    sorted_J = dict(sorted(J.items()))
    Jcombs = dict(sorted(Jcombs.items()))

    return sorted_J,Jcombs
def possible_Js_partial_wavesDs(s1,s2,J_max,C_parity,name1,name2):
    ss = possible_s_twoParticle(s1,s2)
    Js = []
    Jcombs =[]
    for s in ss:
        J,Jcomb = possible_Js_Ds(s,J_max,C_parity,name1,name2)
        Js.append(J)
        Jcombs.append(Jcomb)
    return Js,Jcombs
def joinpartial_waves_Ds(s1,s2,J_max,p1,p2,C_parity,name1,name2):
    Js,Jcombs = possible_Js_partial_wavesDs(s1,s2,J_max,C_parity,name1,name2)
    total_Jcombs = {}
    for index,j in enumerate(Jcombs):
        for key in j.keys():
            if key > J_max:
                continue
                
            if key in total_Jcombs.keys():
                total_Jcombs[key] += Jcombs[index][key]
            else:
                total_Jcombs[key] = Jcombs[index][key]
    parity_factor = p1*p2
    Jparity = {}
    for key in total_Jcombs.keys():
        combs = total_Jcombs[key]
        for index,comb in enumerate(combs):
            pkey = str(key)+ str("+")
            mkey = str(key)+ str("-")
            comb_p = parity_factor*(-1)**comb[0]
            
            if comb_p >0:
                if pkey in Jparity.keys():
                    Jparity[pkey] += [comb]
                else:
                    Jparity[pkey] = [comb]
            else:
                if mkey in Jparity.keys():
                    Jparity[mkey] += [comb]
                else:
                    Jparity[mkey] = [comb]
    return Jparity



        

    

