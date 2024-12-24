import numpy as np
def Isospin(isospin_one,isospin_two):
    limit_1 = abs(isospin_one - isospin_two)
    limit_2 = abs(isospin_one + isospin_two)
    isospin = np.arange(limit_1,limit_2+1)
    return isospin
def Charge_conjugation(charge_conguation_1,charge_conjugation_2):
    return charge_conguation_1*charge_conjugation_2

def Charm(charm_1,charm_2):
    return charm_1+charm_2
def Strangeness(strangeness_1,strangeness_2):
    return strangeness_1+strangeness_2

def Isospin_three_particles(isospin_one,isospin_two,isospin_three):
    list_1 = Isospin(isospin_one,isospin_two)
    isospin = []
    for i in list_1:
        isospin = Isospin(i,isospin_three)
    return isospin
def Charge_conjugation_three_particles(charge_conguation_1,charge_conjugation_2,charge_conjugation_3):
    return charge_conguation_1*charge_conjugation_2*charge_conjugation_3
def Charm_three_particles(charm_1,charm_2,charm_3):
    return charm_1+charm_2+charm_3
def Strangeness_three_particles(strangeness_1,strangeness_2,strangeness_3):
    return strangeness_1+strangeness_2+strangeness_3
    
