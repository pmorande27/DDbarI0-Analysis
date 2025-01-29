import numpy as np
import colors
import spectrum_analysis
import os
import no_int_E_levels
import grouptheory as gt
import particle
def preliminary_analysis_T1pM_corrected():
    irrep = "T1pM-corrected"
    create_structure(24,[8,9,10,11,12],irrep)
    create_structure(16,[8,9,10,11,12],irrep)
    create_structure(20,[8,9,10,11,12,13],irrep)
    figures_save, create_files,save_hist = False,False,False
    spectrumV16t09 = spectrum_analysis.Spectrum(16,9,irrep,create_files,figures_save,save_hist)
    
    spectrumV20t010 = spectrum_analysis.Spectrum(20,10,irrep,create_files,figures_save,save_hist)
    spectrumV24t010 = spectrum_analysis.Spectrum(24,10,irrep,create_files,figures_save,save_hist)
    #spectrumV20t010.automatic_coloring(30)
    #spectrumV24t010.automatic_coloring(30)
    #spectrumV16t09.automatic_coloring(30)
    #[spectrumV20t010.plot_histogram_of_levels_of_color('blue',i) for i in range(10)]
    spectra = [spectrumV24t010,spectrumV20t010,spectrumV16t09]
    Ls = np.linspace(14,26)
    Ps = [np.array([0,0,0]),np.array([1,0,0]),np.array([1,1,0]),np.array([1,1,1])]
    channel = {'Charm':0,'Strange':0,'Isospin':0,'Charm_Isospin':1,'C_parity':-1}
    #spectrumV24t010.plot_histogram_of_levels_of_color('grey',3)
    spectrum_analysis.Spectrum.plot_irrep_mass(spectra,"T1pM",Ls,"000",30,channel,0.74)
def mom_format(mom):
    i = no_int_E_levels.get_mom_from_text(mom)
    i.sort(reverse=True)
    return "".join([str(j) for j in i])
def obtain_sigma_files(Vs=[16,20,24]):
    accessed = {}
    for V in Vs:
        path = f"sigma/sigma_{V}.txt"
        if os.path.exists(path):
            with open(path) as f:
                for line in f:
                    name = line.split("::")[0].strip()
                    rest = line.split("::")[1].strip()
                    momentum = mom_format(name.split("d")[1].split("_")[0])
                   
                    E = rest.split("+/-")[0].strip()
                    err = rest.split("+/-")[1].strip().split(" ")[0]
                    new_line = E + " +/- " + err +"\n"
                    new_path = f"sigma/sigma_E_levels_{momentum}_{V}.txt"
                    if new_path not in accessed.keys():
                        accessed[new_path] = 1
                        with open(new_path,"w") as f2:
                            f2.write(new_line)
                    else:

                        with open(new_path,"a") as f2:
                            f2.write(new_line)
                    print(momentum,E,err)
def preliminary_analysis_T1pM_nosigma():
    irrep = "T1pM-nosigma"
    create_structure(24,[8,9,10,11,12],irrep)
    create_structure(16,[8,9,10,11,12],irrep)
    create_structure(20,[8,9,10,11,12,13],irrep)
    figures_save, create_files,save_hist = False,False,False
    spectrumV20t010 = spectrum_analysis.Spectrum(20,10,irrep,create_files,figures_save,save_hist)
    spectrumV24t010 = spectrum_analysis.Spectrum(24,10,irrep,create_files,figures_save,save_hist)
    spectrumV16t09 = spectrum_analysis.Spectrum(16,9,irrep,create_files,figures_save,save_hist)
    #spectrumV24t010.automatic_coloring(30)
    #spectrumV20t010.automatic_coloring(30)
    #spectrumV16t09.automatic_coloring(30)
    spectra = [spectrumV24t010,spectrumV20t010,spectrumV16t09]
    Ls = np.linspace(14,26)
    Ps = [np.array([0,0,0]),np.array([1,0,0]),np.array([1,1,0]),np.array([1,1,1])]
    channel = {'Charm':0,'Strange':0,'Isospin':0,'Charm_Isospin':1,'C_parity':-1}
    #spectrumV24t010.plot_histogram_of_levels_of_color('grey',3)
    spectrum_analysis.Spectrum.plot_irrep_mass(spectra,"T1pM",Ls,"000",30,channel,0.74)

def preliminary_analysis_T1pM():
    irrep = "T1pM"
    create_structure(24,[8,9,10,11,12],irrep)
    create_structure(16,[8,9,10,11,12],irrep)
    create_structure(20,[8,9,10,11,12,13],irrep)
    figures_save, create_files,save_hist = False,False,False
    spectrumV20t010 = spectrum_analysis.Spectrum(20,10,irrep,create_files,figures_save,save_hist)
    spectrumV16t09 = spectrum_analysis.Spectrum(16,9,irrep,create_files,figures_save,save_hist)
    spectrumV24t010 = spectrum_analysis.Spectrum(24,10,irrep,create_files,figures_save,save_hist)
    #spectrumV24t010.automatic_coloring(30)
    #spectrumV20t010.automatic_coloring(30)
    #spectrumV16t09.automatic_coloring(30)
    spectra = [spectrumV24t010,spectrumV20t010,spectrumV16t09]
    Ls = np.linspace(14,26)
    Ps = [np.array([0,0,0]),np.array([1,0,0]),np.array([1,1,0]),np.array([1,1,1])]
    channel = {'Charm':0,'Strange':0,'Isospin':0,'Charm_Isospin':1,'C_parity':-1}
    #spectrumV24t010.plot_histogram_of_levels_of_color('grey',3)
    spectrum_analysis.Spectrum.plot_irrep_mass(spectra,"T1pM",Ls,"000",30,channel,0.74)
    #spectrumV16t09.plot_histogram_of_levels_of_color('red',1)

def preliminar_analysis_T1mM_nosigma():
    irrep = "T1mM-nosigma-tmax25"
    create_structure(24,[8,9,10,11,12],irrep)
    create_structure(16,[8,9,10,11,12],irrep)
    create_structure(20,[8,9,10,11,12,13],irrep)
    figures_save, create_files,save_hist = False,False,False
    spectrumV24t011 = spectrum_analysis.Spectrum(24,12,irrep,create_files,figures_save,save_hist)
    spectrumV20t010 = spectrum_analysis.Spectrum(20,9,irrep,create_files,figures_save,save_hist)
    spectrumV16t011 = spectrum_analysis.Spectrum(16,9,irrep,create_files,figures_save,save_hist)
    Ls = np.linspace(14,26)
    Ps = [np.array([0,0,0]),np.array([1,0,0]),np.array([1,1,0]),np.array([1,1,1])]
    channel = {'Charm':0,'Strange':0,'Isospin':0,'Charm_Isospin':1,'C_parity':-1}
    spectrumV24t08full = spectrum_analysis.Spectrum(24,8,"T1mM",create_files,figures_save,save_hist)
    spectrumV20t09full = spectrum_analysis.Spectrum(20,9,"T1mM",create_files,figures_save,save_hist)
    spectrumV16t011full = spectrum_analysis.Spectrum(16,11,"T1mM",create_files,figures_save,save_hist)
    #[spectrumV24t011.plot_histogram_of_levels_of_color('grey',i) for i in range(10)]
    spectra = [spectrumV16t011,spectrumV24t011,spectrumV20t010]
    spectrum_analysis.Spectrum.plot_irrep_mass(spectra,"T1mM",Ls,"000",25,channel,0.74)
    spectra = [spectrumV16t011,spectrumV24t011,spectrumV20t010,spectrumV16t011full,spectrumV24t08full,spectrumV20t09full]
    #spectrum_analysis.Spectrum.plot_irrep_mass_superimposed(spectra,"T1mM",Ls,"000",25,channel,0.74)

def preliminar_analysis_T1mM_fewer_pm():
    irrep = "T1mM-fewer-pm"
    create_structure(24,[8,9,10,11,12],irrep)
    create_structure(16,[8,9,10,11,12],irrep)
    create_structure(20,[8,9,10,11,12,13],irrep)
    figures_save, create_files,save_hist = False,False,False
    spectrumV20t010 = spectrum_analysis.Spectrum(20,10,irrep,create_files,figures_save,save_hist)
    Ls = np.linspace(14,26)
    Ps = [np.array([0,0,0]),np.array([1,0,0]),np.array([1,1,0]),np.array([1,1,1])]
    channel = {'Charm':0,'Strange':0,'Isospin':0,'Charm_Isospin':1,'C_parity':-1}
    #[spectrumV20t010.plot_histogram_of_levels_of_color('grey',i) for i in range(10)]
    spectra = [spectrumV20t010]
    #spectrumV20t010.plot_histogram_of_levels_of_color
    spectrum_analysis.Spectrum.plot_irrep_mass(spectra,"T1mM",Ls,"000",25,channel,0.74)
    
    #spectrum_analysis.Spectrum.plot_irrep_mass_superimposed(spectra,"T1mM",Ls,"000",25,channel,0.74)



def preliminar_analysis_T1mM_fewer_djw():
    irrep = "T1mM-fewer-djw-tmax30"
    create_structure(24,[8,9,10,11,12],irrep)
    create_structure(16,[8,9,10,11,12],irrep)
    create_structure(20,[8,9,10,11,12,13],irrep)
    figures_save, create_files,save_hist = False,False,False
    spectrumV24t011 = spectrum_analysis.Spectrum(24,11,irrep,create_files,figures_save,save_hist)
    spectrumV20t010 = spectrum_analysis.Spectrum(20,10,irrep,create_files,figures_save,save_hist)
    spectrumV16t011 = spectrum_analysis.Spectrum(16,11,irrep,create_files,figures_save,save_hist)
    Ls = np.linspace(14,26)
    Ps = [np.array([0,0,0]),np.array([1,0,0]),np.array([1,1,0]),np.array([1,1,1])]
    channel = {'Charm':0,'Strange':0,'Isospin':0,'Charm_Isospin':1,'C_parity':-1}
    #[spectrumV24t011.plot_histogram_of_levels_of_color('grey',i) for i in range(10)]
    spectra = [spectrumV16t011,spectrumV24t011,spectrumV20t010]
    spectrum_analysis.Spectrum.plot_irrep_mass(spectra,"T1mM",Ls,"000",25,channel,0.75)
       
def preliminary_analysis_T1mM():
    create_structure(24,[8,9,10,11,12],"T1mM")
    create_structure(16,[8,9,10,11,12],"T1mM")
    create_structure(20,[8,9,10,11,12],"T1mM")
    figures_save, create_files,save_hist = False,False,False
    spectrumV24t08 = spectrum_analysis.Spectrum(24,8,"T1mM",create_files,figures_save,save_hist)
    spectrumV20t09 = spectrum_analysis.Spectrum(20,9,"T1mM",create_files,figures_save,save_hist)
    spectrumV16t011 = spectrum_analysis.Spectrum(16,11,"T1mM",create_files,figures_save,save_hist)
    #spectrumV16t011.automatic_coloring(24)
    #spectrumV20t09.automatic_coloring(24)
    #spectrumV24t08.automatic_coloring(24)
    Ls = np.linspace(14,26)
    Ps = [np.array([0,0,0]),np.array([1,0,0]),np.array([1,1,0]),np.array([1,1,1])]
    channel = {'Charm':0,'Strange':0,'Isospin':0,'Charm_Isospin':1,'C_parity':-1}
    spectra = [spectrumV16t011,spectrumV24t08,spectrumV20t09]
    spectrum_analysis.Spectrum.plot_irrep_mass(spectra,"T1mM",Ls,"000",25,channel,0.74)

def create_irrep_ni_plus_C_parity(irrep,mom,th):
    channel = {'Charm':0,'Strange':0,'Isospin':0,'Charm_Isospin':1,'C_parity':1}
    L= np.linspace(14,25,100)
    #no_int_E_levels.generate_xml_pairs_no_symmetry(channel,0.76)


    spectrum_analysis.plot_irrep_mass_in_flight(irrep,L,mom,channel,th)
def create_irrep_ni_minus_C_parity(irrep,mom,th):
    channel = {'Charm':0,'Strange':0,'Isospin':0,'Charm_Isospin':1,'C_parity':-1}
    L= np.linspace(14,25,100)
    #no_int_E_levels.generate_xml_pairs_no_symmetry(channel,th)

    #no_int_E_levels.get_levels_irrep_in_flight_sigma(channel,irrep,"energies_no_sym.txt",mom)
    spectrum_analysis.plot_irrep_mass_in_flight(irrep,L,mom,channel,th)




def main():
    #preliminary_analysis_T1mM()
    #preliminar_analysis_T1mM_fewer_djw()
    #preliminar_analysis_T1mM_fewer_pm()
    #preliminar_analysis_T1mM_nosigma()

    #preliminary_analysis_T1pM()
    #preliminary_analysis_T1pM_nosigma()
    #preliminary_analysis_T1pM_corrected()
    
    th = 0.74
    channel = {'Charm':0,'Strange':0,'Isospin':0,'Charm_Isospin':1,'C_parity':-1}
    irreps_rest = ["T1mM","T1pM","T2mM","T2pM","A1mM","A1pM","A2mM","A2pM","EpM","EmM"]
    irreps_unit_mom = ["A1M","A2M","B1M","B2M","E2M"]
    no_int_E_levels.purge_sym_no_int_levels(channel)
    [create_irrep_ni_minus_C_parity(irrep_rest,"000",th) for irrep_rest in irreps_rest]
    [create_irrep_ni_minus_C_parity(irrep_unit_mom,"100",th) for irrep_unit_mom in irreps_unit_mom]

    th = 0.74
    channelp = {'Charm':0,'Strange':0,'Isospin':0,'Charm_Isospin':1,'C_parity':1}
    irreps_rest = ["T1mP","T1pP","T2mP","T2pP","A1mP","A1pP","A2mP","A2pP","EpP","EmP"]
    irreps_unit_mom = ["A1P","A2P","B1P","B2P","E2P"]
    no_int_E_levels.purge_sym_no_int_levels(channelp)
    [create_irrep_ni_plus_C_parity(irrep_rest,"000",th) for irrep_rest in irreps_rest]
    [create_irrep_ni_plus_C_parity(irrep_unit_mom,"100",th) for irrep_unit_mom in irreps_unit_mom]

    #no_int_E_levels.purge_sym_no_int_levels(channel)
    #no_int_E_levels.create_operator_list(channel,"T1pM","000",0.74)
    
    #print(particle.create_table_partial_wave_irrep(channel=channel,Jmax = 3,threshold=0.74,mom = "100"))
    #print(gt.helicities_in_irrep_D4())
    #print(gt.Js_in_irrep_D4(4))
    #print(particle.possible_partial_waves_per_irrep_D4(channel,4,0.74)['A_1'])
    #print(particle.check_partial_waves_in_irrep_D4(channel,4,0.74))

def create_structure(volume,t0s,irrep):
    if not os.path.exists(irrep):
        os.mkdir(irrep)
    if not os.path.exists(f"{irrep}\\Volume_{volume}"):
        os.mkdir(f"{irrep}\\Volume_{volume}")
    for t0 in t0s:
        path2 = f"{irrep}\\Volume_{volume}\\t0{t0}"
        if not os.path.exists(path2):
            os.mkdir(path2)
        path3 = f"{irrep}\\Volume_{volume}\\t0{t0}\\StateColorFiles"
        if not os.path.exists(path3):
            os.mkdir(path3)
        path4 = f"{irrep}\\Volume_{volume}\\t0{t0}\\ZValues"
        if not os.path.exists(path4):
            os.mkdir(path4)
        path5 = f"{irrep}\\Volume_{volume}\\t0{t0}\\ZvaluesRenormalized"
        if not os.path.exists(path5):
            os.mkdir(path5)
        path6 = f"{irrep}\\Volume_{volume}\\t0{t0}\\MassValues"
        if not os.path.exists(path6):
            os.mkdir(path6) 
        path8 = f"{irrep}\\Volume_{volume}\\t0{t0}\\CorrPlots"
        if not os.path.exists(path8):
            os.mkdir(path8)
        path9 = f"{irrep}\\Volume_{volume}\\t0{t0}\\Annotations"
        if not os.path.exists(path9):
            os.mkdir(path9)
        path10 = f"{irrep}\\Volume_{volume}\\t0{t0}\\HistogramPlots"
        if not os.path.exists(path10):
            os.mkdir(path10)
main()
#[create_structure(i,[8,9,10,11,12],"T1mM-fewer-djw-tmax30") for i in [24,16,20]]
