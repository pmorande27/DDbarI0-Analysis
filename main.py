import numpy as np
import colors
import spectrum_analysis
import os
def main():
    #spectrum = spectrum_analysis.Spectrum(16,11,"T1mM",False)
    #spectrumV20t010 = spectrum_analysis.Spectrum(20,10,"T1mM",True,True)
    figures_save, create_files = False, False
    #create_structure(24,[8,9,10,11,12],"T1mM-nosigma-tmax25")
    #create_structure(20,[8,9,10,11,12],"T1mM-nosigma-tmax25")
    #create_structure(16,[8,9,10,11,12],"T1mM-nosigma-tmax25")
    #spectrumV20t09nosigma = spectrum_analysis.Spectrum(20,9,"T1mM-nosigma-tmax25",create_files,figures_save)
    #spectrumV24t012nosigma = spectrum_analysis.Spectrum(24,12,"T1mM-nosigma-tmax25",create_files,figures_save)
    #spectrumV16t09nosigma = spectrum_analysis.Spectrum(16,9,"T1mM-nosigma-tmax25",create_files,figures_save)
    #spectrumV16t09nosigma.automatic_coloring(24)
    #spectra = [spectrumV20t09nosigma,spectrumV24t012nosigma,spectrumV16t09nosigma]
    
    Ls = np.linspace(14,26)
    Ps = [np.array([0,0,0]),np.array([1,0,0]),np.array([1,1,0]),np.array([1,1,1])]
    #[spectrumV24t012nosigma.plot_histogram_state(i) for i in range(24)]
    #spectrum_analysis.Spectrum.plot_irrep_mass(spectra,"T_1^-",Ls,Ps,24)

    #spectrumV24t012nosigma.automatic_coloring(24)
    #spectrumV20t09 = spectrum_analysis.Spectrum(16,9,"T1mM",create_files,figures_save)
    spectrumV20t010 = spectrum_analysis.Spectrum(20,10,"T1mM",create_files,figures_save)
    #spectrumV20t011 = spectrum_analysis.Spectrum(16,11,"T1mM",create_files,figures_save)
    #spectrumV20t012 = spectrum_analysis.Spectrum(16,12,"T1mM",create_files,figures_save)

    #spectra = [spectrumV20t09,spectrumV20t010,spectrumV20t011,spectrumV20t012]

    spectrumV16t011 = spectrum_analysis.Spectrum(16,11,"T1mM",False,False)
    #spectrumV16t012 = spectrum_analysis.Spectrum(16,12,"T1mM",False,False)
    #spectrumV16t09 = spectrum_analysis.Spectrum(16,9,"T1mM",False,False)
    #spectra = [spectrumV16t09,spectrumV16t011,spectrumV16t012]
    spectrumV20t09 = spectrum_analysis.Spectrum(20,9,"T1mM",False,False)
    spectrumV24t08 = spectrum_analysis.Spectrum(24,8,"T1mM",False,False)
    spectrumV16t010 = spectrum_analysis.Spectrum(16,10,"T1mM",False,False)
    #spectra = [spectrumV20t09,spectrumV24t08,spectrumV16t010]
    
    #[create_structure(i,[8,9,10,11,12],"T1mM") for i in [16,20,24]]

    #colors.save_color_code_state(1,2,spectrumV16t010)
    #spectrumV24t011_fewer = spectrum_analysis.Spectrum(24,11,"T1mM-fewer",False,False)
    #spectrumV24t011_fewer25 = spectrum_analysis.Spectrum(24,11,"T1mM-fewer-tmax25",False,False)
    #spectrumV24t08 = spectrum_analysis.Spectrum(24,8,"T1mM",False,False)
    #spectra = [spectrumV20t09,spectrumV24t08]
    #spectrumV24t010_fewer25 = spectrum_analysis.Spectrum(24,10,"T1mM-fewer-tmax25",False,False)
    #spectra = [spectrumV24t011_fewer,spectrumV24t011_fewer25,spectrumV24t08,spectrumV24t010_fewer25,spectrumV24t012nosigma]
    spectra = [spectrumV16t011,spectrumV20t010,spectrumV24t08]
    #spectrumV20t010.automatic_coloring(20)
    #[spectrum.automatic_coloring(20) for spectrum in spectra]
    #spectrumV16t09.plot_histogram_state(1)
    channel = {'Charm':0,'Strange':0,'Isospin':0,'Charm_Isospin':1,'C_parity':-1}

    spectrum_analysis.Spectrum.plot_irrep_mass(spectra,"T_1^-",Ls,Ps,24,channel)
    #spectrum_analysis.Spectrum.plot_spectrum_multiple(spectra,24)
    #print(spectrum_analysis.plot_sigma_E_levels("\psi",[16],[[0,0,0]],"T_1^-"))
  

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
main()
#[create_structure(i,[8,9,10,11,12],"T1mM-fewer-djw-tmax30") for i in [24,16,20]]
