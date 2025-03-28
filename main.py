import numpy as np
import colors
import spectrum_analysis
import os
import no_int_E_levels
import grouptheory as gt
import particle
import Set_Up_Functions
import xml_functions
import matplotlib.pyplot as plt
import plot_spectrum as ps
import plot_amplitudes as pa
def preliminary_analysis_T1mM_A2mM_A1():
    irrep = "A2mM-nosigma-reduced"
    create_structure(24,[8,9,10,11,12,13,14],irrep)
    create_structure(16,[8,9,10,11,12],irrep)
    create_structure(20,[8,9,10,11,12,13,14,15],irrep)
    figures_save, create_files,save_hist = False,False,False
    #spectrumV24 = spectrum_analysis.Spectrum(24,13,irrep,create_files,figures_save,save_hist)
    spectrumV20 = spectrum_analysis.Spectrum(20,14,irrep,create_files,figures_save,save_hist)
    spectrumV16 = spectrum_analysis.Spectrum(16,12,irrep,create_files,figures_save,save_hist)
    spectrumV24 = spectrum_analysis.Spectrum(24,14,irrep,create_files,figures_save,save_hist)

    
    spectrumV20.automatic_coloring(8)
    spectrumV16.automatic_coloring(8)
    spectrumV24.automatic_coloring(8)
    Ls = np.linspace(14,26)
    channel = {'Charm':0,'Strange':0,'Isospin':0,'Charm_Isospin':1,'C_parity':-1}
    spectra = [spectrumV16,spectrumV20,spectrumV24]
    included = [[2,3,4,6],[2,3,5,4,11,6,8],[2,3,4,13]]
    #spectrum_analysis.Spectrum.plot_irrep_mass_no_sigma(spectra,"A2mM",Ls,"000",7,channel,0.78)
    included = [[1],[1,4],[1,2]]
    #spectrum_analysis.Spectrum.plot_irrep_mass_no_sigma_grey_out(spectra,"A2mM",Ls,"000",7,channel,0.78,included)
    fig,ax = plt.subplots(1,3,sharex=True,sharey=True)

    ps.plot_irrep_mass_with_get_finite_E_vs_L_ax(spectra,"A2mM",Ls,"000",7,channel,0.78,"./Data/levels_A2mM",ax[0])
    #ps.plot_irrep_mass_no_sigma_grey_out(spectra,"A2mM",Ls,"000",7,channel,0.78,included,ax[0])
    irrep = "T1mM-nosigma-reduced"
    create_structure(24,[8,9,10,11,12,13],irrep)
    create_structure(16,[8,9,10,11,12],irrep)
    create_structure(20,[8,9,10,11,12,13],irrep)
    included = [[1,2,3,4,6],[1,2,3,5,4,11,6,8],[1,2,3,4,13,7,8]]

    figures_save, create_files,save_hist = False,False,False
    spectrumV16 = spectrum_analysis.Spectrum(16,12,irrep,create_files,figures_save,save_hist)
    spectrumV20 = spectrum_analysis.Spectrum(20,10,irrep,create_files,figures_save,save_hist)
    #spectrumV24 = spectrum_analysis.Spectrum(24,9,irrep,create_files,figures_save,save_hist)
    spectrumV24 = spectrum_analysis.Spectrum(24,12,irrep,create_files,figures_save,save_hist)

    
    #spectrumV20.automatic_coloring(30)
    #spectrumV16.automatic_coloring(30)
    #spectrumV24.automatic_coloring(30)
    Ls = np.linspace(14,26)
    channel = {'Charm':0,'Strange':0,'Isospin':0,'Charm_Isospin':1,'C_parity':-1}
    spectra = [spectrumV16,spectrumV24,spectrumV20]
    included = [[1,2,4,6],[1,2,3,4,11,6,8],[1,2,3,13]]
    #spectrum_analysis.Spectrum.plot_irrep_mass_no_sigma(spectra,"T1mM",Ls,"000",20,channel,0.74)
    #ps.plot_irrep_mass_no_sigma_grey_out(spectra,"T1mM",Ls,"000",19,channel,0.72,included)
    #fig,axis = plt.subplots()
    
    #spectrumV20.plot_nine_levels([1,2,3,4,5,6,7,8,9],"$T_1^{--} [000]$")
    name = f"prin_corr_fit_t0{spectrumV16.t0}_reorder_state{0}.ax"
    path =  f"{spectrumV16.irrep}\\Volume_{spectrumV16.volume}\\t0{spectrumV16.t0}\\PrinCorrPlots\\{name}"
    #spectrum_analysis.plot_ax_file_in_fig(spectrumV16,path,axis,state=0)
    ps.plot_irrep_mass_with_get_finite_E_vs_L_ax(spectra,"T1mM",Ls,"000",30,channel,0.75,"./Data/levels_T1mM",ax[1])
    #ps.plot_irrep_mass_no_sigma_grey_out(spectra,"T1mM",Ls,"000",19,channel,0.75,included,ax[1])

    irrep = "D4A1M-nosigma"
    create_structure(24,[8,9,10,11,12,13,14,15,16],irrep)
    create_structure(16,[8,9,10,11,12,13,14,15],irrep)
    create_structure(20,[8,9,10,11,12,13,14,15,16],irrep)
    figures_save, create_files,save_hist = False,False,False
    spectrumV16 = spectrum_analysis.Spectrum(16,14,irrep,create_files,figures_save,save_hist)
    spectrumV20 = spectrum_analysis.Spectrum(20,14,irrep,create_files,figures_save,save_hist)
    #spectrumV16.plot_histogram_of_levels_of_color('green',2)
  
    spectrumV24 = spectrum_analysis.Spectrum(24,12,irrep,create_files,figures_save,save_hist)

    #spectrumV20.plot_nine_levels([0,1,2,3,4,5,6,7,8,9],"$A_1^{-} [100]$")
    #spectrumV20.automatic_coloring(30)
    #spectrumV16.automatic_coloring(30)
    #spectrumV24.automatic_coloring(30)
    Ls = np.linspace(14,26)
    channel = {'Charm':0,'Strange':0,'Isospin':0,'Charm_Isospin':1,'C_parity':-1}
    spectra = [spectrumV16,spectrumV20,spectrumV24]
    included = [[1,2,3,4,13,7,8]]
    #spectrum_analysis.Spectrum.plot_irrep_mass(spectra,"A1M",Ls,"100",20,channel,0.74)
    ps.plot_irrep_mass_with_get_finite_E_vs_L_ax(spectra,"A1M",Ls,"100",19,channel,0.78,"./Data/levels_A1",ax[2])
    
    
    plt.show()
    

def preliminary_analysis_T1mM_nosigma_improved():
    irrep = "T1mM-nosigma-improved"
    create_structure(24,[8,9,10,11,12],irrep)
    create_structure(16,[8,9,10,11,12],irrep)
    create_structure(20,[8,9,10,11,12,13],irrep)
    figures_save, create_files,save_hist = False,False,False
    spectrumV16 = spectrum_analysis.Spectrum(16,10,irrep,create_files,figures_save,save_hist)
    spectrumV20 = spectrum_analysis.Spectrum(20,12,irrep,False,False,False)
    #spectrumV24 = spectrum_analysis.Spectrum(24,9,irrep,create_files,figures_save,save_hist)
    spectrumV24 = spectrum_analysis.Spectrum(24,12,irrep,create_files,figures_save,save_hist)

    
    spectrumV20.automatic_coloring(30)
    #spectrumV16.automatic_coloring(30)
    #spectrumV24.automatic_coloring(30)
    
    
    
    spectrumV24.automatic_coloring(30)
    
    Ls = np.linspace(14,26)
    channel = {'Charm':0,'Strange':0,'Isospin':0,'Charm_Isospin':1,'C_parity':-1}
    spectra = [spectrumV16,spectrumV24,spectrumV20]
    spectrum_analysis.Spectrum.plot_irrep_mass_no_sigma(spectra,"T1mM",Ls,"000",30,channel,0.73)


    """

    """
def preliminary_analysis_T1mM_nosigma_reduced():
    irrep = "T1mM-nosigma-reduced"
    fig,ax = plt.subplots()
    figures_save, create_files,save_hist = False,False,False
    spectrumV16 = spectrum_analysis.Spectrum(16,12,irrep,create_files,figures_save,save_hist)
    spectrumV20 = spectrum_analysis.Spectrum(20,10,irrep,create_files,figures_save,save_hist)
    #spectrumV24 = spectrum_analysis.Spectrum(24,9,irrep,create_files,figures_save,save_hist)
    spectrumV24 = spectrum_analysis.Spectrum(24,12,irrep,create_files,figures_save,save_hist)

    
    #spectrumV20.automatic_coloring(30)
    #spectrumV16.automatic_coloring(30)
    #spectrumV24.automatic_coloring(30)
    Ls = np.linspace(14,26)
    channel = {'Charm':0,'Strange':0,'Isospin':0,'Charm_Isospin':1,'C_parity':-1}
    spectra = [spectrumV16,spectrumV24,spectrumV20]
    included = [[1,2,4,6],[1,2,3,4,11,6,8],[1,2,3,13]]
    #spectrum_analysis.Spectrum.plot_irrep_mass_no_sigma(spectra,"T1mM",Ls,"000",20,channel,0.74)
    #ps.plot_irrep_mass_no_sigma_grey_out(spectra,"T1mM",Ls,"000",19,channel,0.72,included)
    #fig,axis = plt.subplots()
    
    #spectrumV20.plot_nine_levels([1,2,3,4,5,6,7,8,9],"$T_1^{--} [000]$")
    name = f"prin_corr_fit_t0{spectrumV16.t0}_reorder_state{0}.ax"
    path =  f"{spectrumV16.irrep}\\Volume_{spectrumV16.volume}\\t0{spectrumV16.t0}\\PrinCorrPlots\\{name}"
    #spectrum_analysis.plot_ax_file_in_fig(spectrumV16,path,axis,state=0)

    ps.plot_irrep_mass_with_get_finite_E_vs_L_ax(spectra,"T1mM",Ls,"000",19,channel,0.75,"./Data/levels_T1mM",ax)
    plt.show()
    #spectrum_analysis.Spectrum.plot_irrep_mass_no_sigma(spectra,"T1mM",Ls,"000",20,channel,0.74)
    #fig,axis = plt.subplots()
    #ps.plot_irrep_mass(spectra,"T1mM",Ls,"000",19,channel,0.75)
    #
    
    #spectrumV20.plot_nine_levels([1,2,3,4,5,6,7,8,9],"$T_1^{--} [000]$")
    name = f"prin_corr_fit_t0{spectrumV16.t0}_reorder_state{0}.ax"
    path =  f"{spectrumV16.irrep}\\Volume_{spectrumV16.volume}\\t0{spectrumV16.t0}\\PrinCorrPlots\\{name}"
    #spectrum_analysis.plot_ax_file_in_fig(spectrumV16,path,axis,state=0)
    #plt.show()
    #spectrum_analysis.Spectrum.plot_irrep_mass_with_get_finite(spectra,"T1mM",Ls,"000",19,channel,0.75,"./Data/output_d000_T1m.spectrum")

    
    """    
    Ls = np.linspace(14,26)
    channel = {'Charm':0,'Strange':0,'Isospin':0,'Charm_Isospin':1,'C_parity':-1}
    spectra = [spectrumV16,spectrumV24,spectrumV20]
    spectrum_analysis.Spectrum.plot_irrep_mass_no_sigma(spectra,"T1mM",Ls,"000",20,channel,0.72)

    """
def preliminary_analysis_T1mM_nosigma_reduced_2():
    irrep = "T1mM-nosigma-reduced-2"
    create_structure(24,[8,9,10,11,12,13],irrep)
    create_structure(16,[8,9,10,11,12],irrep)
    create_structure(20,[8,9,10,11,12,13],irrep)
    figures_save, create_files,save_hist = True,True,True
    #spectrumV16 = spectrum_analysis.Spectrum(16,12,irrep,create_files,figures_save,save_hist)
    spectrumV20 = spectrum_analysis.Spectrum(20,10,irrep,create_files,figures_save,save_hist)
    #spectrumV24 = spectrum_analysis.Spectrum(24,9,irrep,create_files,figures_save,save_hist)
    #spectrumV24 = spectrum_analysis.Spectrum(24,12,irrep,create_files,figures_save,save_hist)

    
    spectrumV20.automatic_coloring(30)
    #spectrumV16.automatic_coloring(30)
    #spectrumV24.automatic_coloring(30)
    Ls = np.linspace(14,26)
    channel = {'Charm':0,'Strange':0,'Isospin':0,'Charm_Isospin':1,'C_parity':-1}
    spectra = [spectrumV20]
    included = [[1,2,3,4,13,7,8]]
    spectrum_analysis.Spectrum.plot_irrep_mass_no_sigma(spectra,"T1mM",Ls,"000",20,channel,0.74)
    #spectrum_analysis.Spectrum.plot_irrep_mass_no_sigma_grey_out(spectra,"T1mM",Ls,"000",19,channel,0.72,included)
    #spectrum_analysis.Spectrum.plot_irrep_mass_with_get_finite(spectra,"T1mM",Ls,"000",19,channel,0.73,"./Data/output_d000_T1m.spectrum")

    
    """    
    Ls = np.linspace(14,26)
    channel = {'Charm':0,'Strange':0,'Isospin':0,'Charm_Isospin':1,'C_parity':-1}
    spectra = [spectrumV16,spectrumV24,spectrumV20]
    spectrum_analysis.Spectrum.plot_irrep_mass_no_sigma(spectra,"T1mM",Ls,"000",20,channel,0.72)

    """
def preliminary_analysis_D4A1M():
    irrep = "D4A1M"
    create_structure(24,[8,9,10,11,12,13,14,15,16],irrep)
    create_structure(16,[8,9,10,11,12,13,14,15],irrep)
    create_structure(20,[8,9,10,11,12,13,14,15,16],irrep)
    figures_save, create_files,save_hist = False,False,False
    spectrumV16 = spectrum_analysis.Spectrum(16,13,irrep,create_files,figures_save,save_hist)
    #spectrumV16.plot_histogram_of_levels_of_color('green',0)
  
    #spectrumV24 = spectrum_analysis.Spectrum(24,12,irrep,create_files,figures_save,save_hist)

    
    #spectrumV20.automatic_coloring(30)
    spectrumV16.automatic_coloring(30)
    #spectrumV24.automatic_coloring(30)
    Ls = np.linspace(14,26)
    channel = {'Charm':0,'Strange':0,'Isospin':0,'Charm_Isospin':1,'C_parity':-1}
    spectra = [spectrumV16]
    included = [[1,2,3,4,13,7,8]]
    spectrum_analysis.Spectrum.plot_irrep_mass(spectra,"A1M",Ls,"100",20,channel,0.74)
    #spectrum_analysis.Spectrum.plot_irrep_mass_no_sigma_grey_out(spectra,"T1mM",Ls,"000",19,channel,0.72,included)
    #spectrum_analysis.Spectrum.plot_irrep_mass_with_get_finite(spectra,"T1mM",Ls,"000",19,channel,0.73,"./Data/output_d000_T1m.spectrum")
def preliminary_analysis_D4A1M_nosigma():
    irrep = "D4A1M-nosigma"
    create_structure(24,[8,9,10,11,12,13,14,15,16],irrep)
    create_structure(16,[8,9,10,11,12,13,14,15],irrep)
    create_structure(20,[8,9,10,11,12,13,14,15,16],irrep)
    figures_save, create_files,save_hist = False,False,False
    spectrumV16 = spectrum_analysis.Spectrum(16,14,irrep,create_files,figures_save,save_hist)
    spectrumV20 = spectrum_analysis.Spectrum(20,14,irrep,create_files,figures_save,save_hist)
    #spectrumV16.plot_histogram_of_levels_of_color('green',2)
  
    spectrumV24 = spectrum_analysis.Spectrum(24,12,irrep,create_files,figures_save,save_hist)

    #spectrumV20.plot_nine_levels([0,1,2,3,4,5,6,7,8,9],"$A_1^{-} [100]$")
    #spectrumV20.automatic_coloring(30)
    #spectrumV16.automatic_coloring(30)
    #spectrumV24.automatic_coloring(30)
    Ls = np.linspace(14,26)
    channel = {'Charm':0,'Strange':0,'Isospin':0,'Charm_Isospin':1,'C_parity':-1}
    spectra = [spectrumV16,spectrumV20,spectrumV24]
    included = [[1,2,3,4,13,7,8]]
    #spectrum_analysis.Spectrum.plot_irrep_mass(spectra,"A1M",Ls,"100",20,channel,0.74)
    spectrum_analysis.Spectrum.plot_irrep_mass_with_get_finite(spectra,"A1M",Ls,"100",19,channel,0.72,"./Data/output_d001_A1.spectrum")
    #spectrum_analysis.Spectrum.plot_irrep_mass_with_get_finite(spectra,"T1mM",Ls,"000",19,channel,0.73,"./Data/output_d000_T1m.spectrum")
def preliminary_analysis_D4A1M_nosigma_2():
    irrep = "D4A1M-nosigma-2"
    create_structure(24,[8,9,10,11,12,13,14,15,16],irrep)
    create_structure(16,[8,9,10,11,12,13,14,15],irrep)
    create_structure(20,[8,9,10,11,12,13,14,15,16],irrep)
    figures_save, create_files,save_hist = False,False,False
    #spectrumV16 = spectrum_analysis.Spectrum(16,14,irrep,create_files,figures_save,save_hist)
    #spectrumV20 = spectrum_analysis.Spectrum(20,14,irrep,create_files,figures_save,save_hist)
    #spectrumV16.plot_histogram_of_levels_of_color('green',2)
  
    spectrumV24 = spectrum_analysis.Spectrum(24,14,irrep,create_files,figures_save,save_hist)
    spectrumV24.automatic_coloring(20)


    #spectrumV20.plot_nine_levels([0,1,2,3,4,5,6,7,8,9],"$A_1^{-} [100]$")
    #spectrumV20.automatic_coloring(30)
    #spectrumV16.automatic_coloring(30)
    #spectrumV24.automatic_coloring(30)
    Ls = np.linspace(14,26)
    channel = {'Charm':0,'Strange':0,'Isospin':0,'Charm_Isospin':1,'C_parity':-1}
    spectra = [spectrumV24]
    included = [[1,2,3,4,13,7,8]]
    #spectrum_analysis.Spectrum.plot_irrep_mass(spectra,"A1M",Ls,"100",20,channel,0.74)
    spectrum_analysis.Spectrum.plot_irrep_mass_with_get_finite(spectra,"A1M",Ls,"100",19,channel,0.72,"./Data/output_d001_A1.spectrum")
    #spectrum_analysis.Spectrum.plot_irrep_mass_with_get_finite(spectra,"T1mM",Ls,"000",19,channel,0.73,"./Data/output_d000_T1m.spectrum")
def preliminary_analysis_D4A1M_nosigma_reduced():
    irrep = "D4A1M-nosigma-reduced"
    create_structure(24,[8,9,10,11,12,13,14,15,16],irrep)
    create_structure(16,[8,9,10,11,12,13,14,15],irrep)
    create_structure(20,[8,9,10,11,12,13,14,15,16],irrep)
    figures_save, create_files,save_hist = True,True,True
    spectrumV16 = spectrum_analysis.Spectrum(16,12,irrep,create_files,figures_save,save_hist)
    #spectrumV16.plot_histogram_of_levels_of_color('green',2)
  
    #spectrumV24 = spectrum_analysis.Spectrum(24,12,irrep,create_files,figures_save,save_hist)

    
    #spectrumV20.automatic_coloring(30)
    spectrumV16.automatic_coloring(20)
    #spectrumV24.automatic_coloring(30)
    Ls = np.linspace(14,26)
    channel = {'Charm':0,'Strange':0,'Isospin':0,'Charm_Isospin':1,'C_parity':-1}
    spectra = [spectrumV16]
    included = [[1,2,3,4,13,7,8]]
    #spectrum_analysis.Spectrum.plot_irrep_mass(spectra,"A1M",Ls,"100",20,channel,0.74)
    spectrum_analysis.Spectrum.plot_irrep_mass_no_sigma(spectra,"A1M",Ls,"100",10,channel,0.77)
    #spectrum_analysis.Spectrum.plot_irrep_mass_with_get_finite(spectra,"T1mM",Ls,"000",19,channel,0.73,"./Data/output_d000_T1m.spectrum")

    
def preliminary_analysis_T1mM_nosigma_reduced_3():
    irrep = "T1mM-nosigma-reduced-3"
    irrep2 = "T1mM-nosigma-reduced"
    create_structure(24,[8,9,10,11,12,13,14,15,16],irrep)
    create_structure(16,[8,9,10,11,12],irrep)
    create_structure(20,[8,9,10,11,12,13,14,15,16],irrep)
    figures_save, create_files,save_hist = False,False,False
    spectrumV16 = spectrum_analysis.Spectrum(16,12,irrep2,create_files,figures_save,save_hist)
    spectrumV20 = spectrum_analysis.Spectrum(20,11,irrep,create_files,figures_save,save_hist)
    spectrumV24 = spectrum_analysis.Spectrum(24,12,irrep2,create_files,figures_save,save_hist)
    #spectrumV24 = spectrum_analysis.Spectrum(24,12,irrep,create_files,figures_save,save_hist)

    
    spectrumV20.automatic_coloring(30)
    #spectrumV16.automatic_coloring(30)
    #spectrumV24.automatic_coloring(30)
    Ls = np.linspace(14,26)
    channel = {'Charm':0,'Strange':0,'Isospin':0,'Charm_Isospin':1,'C_parity':-1}
    spectra = [spectrumV20,spectrumV16,spectrumV24]
    included = [[1,2,3,4,13,7,8]]
    spectrum_analysis.Spectrum.plot_irrep_mass_with_get_finite(spectra,"T1mM",Ls,"000",19,channel,0.75,"./Data/output_d000_T1m.spectrum")
    #spectrum_analysis.Spectrum.plot_irrep_mass_no_sigma_grey_out(spectra,"T1mM",Ls,"000",19,channel,0.72,included)
    #spectrum_analysis.Spectrum.plot_irrep_mass_with_get_finite(spectra,"T1mM",Ls,"000",19,channel,0.73,"./Data/output_d000_T1m.spectrum")

    
    """    
    Ls = np.linspace(14,26)
    channel = {'Charm':0,'Strange':0,'Isospin':0,'Charm_Isospin':1,'C_parity':-1}
    spectra = [spectrumV16,spectrumV24,spectrumV20]
    spectrum_analysis.Spectrum.plot_irrep_mass_no_sigma(spectra,"T1mM",Ls,"000",20,channel,0.72)

    """

    
    
def preliminary_analysis_T1pM_nosigma_reduced():
    irrep = "T1pM-nosigma-reduced"
    create_structure(24,[8,9,10,11,12,13],irrep)
    create_structure(16,[8,9,10,11,12],irrep)
    create_structure(20,[8,9,10,11,12,13],irrep)
    figures_save, create_files,save_hist = True,True,True
    spectrumV16 = spectrum_analysis.Spectrum(16,11,irrep,create_files,figures_save,save_hist)
    spectrumV20 = spectrum_analysis.Spectrum(20,10,irrep,create_files,figures_save,save_hist)
    #spectrumV24 = spectrum_analysis.Spectrum(24,9,irrep,create_files,figures_save,save_hist)
    #spectrumV24 = spectrum_analysis.Spectrum(24,12,irrep,create_files,figures_save,save_hist)

    
    spectrumV20.automatic_coloring(30)
    spectrumV16.automatic_coloring(30)
    #spectrumV24.automatic_coloring(30)
    Ls = np.linspace(14,26)
    channel = {'Charm':0,'Strange':0,'Isospin':0,'Charm_Isospin':1,'C_parity':-1}
    spectra = [spectrumV16,spectrumV20]
    included = [[2,3,4,6],[2,3,5,4,11,6,8],[2,3,4,13]]
    spectrum_analysis.Spectrum.plot_irrep_mass_no_sigma(spectra,"T1pM",Ls,"000",20,channel,0.74)
    #spectrum_analysis.Spectrum.plot_irrep_mass_with_get_finite(spectra,"T1mM",Ls,"000",19,channel,0.72,"./Data/output_d000_T1m.spectrum")
    
    
    """    
    Ls = np.linspace(14,26)
    channel = {'Charm':0,'Strange':0,'Isospin':0,'Charm_Isospin':1,'C_parity':-1}
    spectra = [spectrumV16,spectrumV24,spectrumV20]
    spectrum_analysis.Spectrum.plot_irrep_mass_no_sigma(spectra,"T1mM",Ls,"000",20,channel,0.72)

    """
def preliminary_analysis_A2mM_nosigma_reduced():
    irrep = "A2mM-nosigma-reduced"
    volumes = [24,20,16]
    t0s = [14,14,12]
    PathDatas = []
    PathsOpss = []
    for i in range(3):
        PathDatas.append(f"../A2mM-nosigma-reduced/Volume_{volumes[i]}/t0{t0s[i]}")
        PathsOpss.append(f"../A2mM-nosigma-reduced/Volume_{volumes[i]}/ops.txt")

    create_structure(24,[8,9,10,11,12,13,14],irrep)
    create_structure(16,[8,9,10,11,12],irrep)
    create_structure(20,[8,9,10,11,12,13,14,15],irrep)
    figures_save, create_files,save_hist = True,True,True
    #spectrumV24 = spectrum_analysis.Spectrum(24,13,irrep,create_files,figures_save,save_hist)
    spectrumV20 = spectrum_analysis.Spectrum(20,14,irrep,create_files,figures_save,save_hist,PathDatas[1],PathsOpss[1])
    spectrumV16 = spectrum_analysis.Spectrum(16,12,irrep,create_files,figures_save,save_hist,PathDatas[2],PathsOpss[2])
    spectrumV24 = spectrum_analysis.Spectrum(24,14,irrep,create_files,figures_save,save_hist,PathDatas[0],PathsOpss[0])

    
    spectrumV20.automatic_coloring(8)
    spectrumV16.automatic_coloring(8)
    spectrumV24.automatic_coloring(8)
    Ls = np.linspace(14,26)
    channel = {'Charm':0,'Strange':0,'Isospin':0,'Charm_Isospin':1,'C_parity':-1}
    spectra = [spectrumV16,spectrumV20,spectrumV24]
    included = [[2,3,4,6],[2,3,5,4,11,6,8],[2,3,4,13]]
    #spectrum_analysis.Spectrum.plot_irrep_mass_no_sigma(spectra,"A2mM",Ls,"000",7,channel,0.78)
    included = [[0],[0,1],[0,1]]
    #spectrum_analysis.Spectrum.plot_irrep_mass_no_sigma_grey_out(spectra,"A2mM",Ls,"000",7,channel,0.78,included)
    fig,ax = plt.subplots()
    ps.plot_irrep_mass_with_get_finite_ax(spectra,"A2mM",Ls,"000",7,channel,0.78,"./Data/output_d000_A2m.spectrum",ax)
    plt.show()
    
    """    
    Ls = np.linspace(14,26)
    channel = {'Charm':0,'Strange':0,'Isospin':0,'Charm_Isospin':1,'C_parity':-1}
    spectra = [spectrumV16,spectrumV24,spectrumV20]
    spectrum_analysis.Spectrum.plot_irrep_mass_no_sigma(spectra,"T1mM",Ls,"000",20,channel,0.72)

    """
def preliminary_analysis_A2mM_nosigma():
    irrep = "A2mM-nosigma"
    create_structure(24,[8,9,10,11,12,13],irrep)
    create_structure(16,[8,9,10,11,12],irrep)
    create_structure(20,[8,9,10,11,12,13],irrep)
    figures_save, create_files,save_hist = False,False,False
    spectrumV24 = spectrum_analysis.Spectrum(24,13,irrep,create_files,figures_save,save_hist)
    spectrumV20 = spectrum_analysis.Spectrum(20,10,irrep,create_files,figures_save,save_hist)
    spectrumV16 = spectrum_analysis.Spectrum(16,12,irrep,create_files,figures_save,save_hist)
    #spectrumV24 = spectrum_analysis.Spectrum(24,12,irrep,create_files,figures_save,save_hist)

    
    spectrumV20.automatic_coloring(8)
    spectrumV16.automatic_coloring(8)
    spectrumV24.automatic_coloring(8)
    Ls = np.linspace(14,26)
    channel = {'Charm':0,'Strange':0,'Isospin':0,'Charm_Isospin':1,'C_parity':-1}
    spectra = [spectrumV16,spectrumV24,spectrumV20]
    included = [[2,3,4,6],[2,3,5,4,11,6,8],[2,3,4,13]]
    spectrum_analysis.Spectrum.plot_irrep_mass_with_get_finite(spectra,"A2mM",Ls,"000",8,channel,0.78,"./Data/output_d000_A2m.spectrum")
    #spectrum_analysis.Spectrum.plot_irrep_mass_with_get_finite(spectra,"T1mM",Ls,"000",19,channel,0.72,"./Data/output_d000_T1m.spectrum")
    
    
    """    
    Ls = np.linspace(14,26)
    channel = {'Charm':0,'Strange':0,'Isospin':0,'Charm_Isospin':1,'C_parity':-1}
    spectra = [spectrumV16,spectrumV24,spectrumV20]
    spectrum_analysis.Spectrum.plot_irrep_mass_no_sigma(spectra,"T1mM",Ls,"000",20,channel,0.72)

    """
def preliminary_analysis_A2mM():
    irrep = "A2mM"
    create_structure(24,[8,9,10,11,12,13],irrep)
    create_structure(16,[8,9,10,11,12],irrep)
    create_structure(20,[8,9,10,11,12,13],irrep)
    figures_save, create_files,save_hist = True,True,True
    spectrumV16 = spectrum_analysis.Spectrum(20,11,irrep,create_files,figures_save,save_hist)
    #spectrumV20 = spectrum_analysis.Spectrum(20,11,irrep,create_files,figures_save,save_hist)
    #spectrumV24 = spectrum_analysis.Spectrum(24,12,irrep,create_files,figures_save,save_hist)
    #spectrumV24 = spectrum_analysis.Spectrum(24,12,irrep,create_files,figures_save,save_hist)

    
    #spectrumV20.automatic_coloring(8)
    spectrumV16.automatic_coloring(8)
    #spectrumV24.automatic_coloring(8)
    Ls = np.linspace(14,26)
    channel = {'Charm':0,'Strange':0,'Isospin':0,'Charm_Isospin':1,'C_parity':-1}
    spectra = [spectrumV16]
    included = [[2,3,4,6],[2,3,5,4,11,6,8],[2,3,4,13]]
    spectrum_analysis.Spectrum.plot_irrep_mass_with_get_finite(spectra,"A2mM",Ls,"000",8,channel,0.78,"./Data/output_d000_A2m.spectrum")
    #spectrum_analysis.Spectrum.plot_irrep_mass_with_get_finite(spectra,"T1mM",Ls,"000",19,channel,0.72,"./Data/output_d000_T1m.spectrum")
    
    
    """    
    Ls = np.linspace(14,26)
    channel = {'Charm':0,'Strange':0,'Isospin':0,'Charm_Isospin':1,'C_parity':-1}
    spectra = [spectrumV16,spectrumV24,spectrumV20]
    spectrum_analysis.Spectrum.plot_irrep_mass_no_sigma(spectra,"T1mM",Ls,"000",20,channel,0.72)

    """

    
    
    
def preliminary_analysis_T1mM_nosigma_fewer():
    irrep = "T1mM-nosigma-fewer"
    create_structure(24,[8,9,10,11,12],irrep)
    create_structure(16,[8,9,10,11,12],irrep)
    create_structure(20,[8,9,10,11,12,13],irrep)
    figures_save, create_files,save_hist = False,False,False
    
    spectrumV24 = spectrum_analysis.Spectrum(24,10,irrep,create_files,figures_save,save_hist)
    colors.create_table_operators("T1mM-nosigma-fewer",24)
    colors.create_table_operators("T1mM-nosigma-fewer",20)
    colors.create_table_operators("T1mM-nosigma-fewer",16)
    spectrumV20 = spectrum_analysis.Spectrum(20,11,irrep,create_files,figures_save,save_hist)
    spectrumV16 = spectrum_analysis.Spectrum(16,11,irrep,create_files,figures_save,save_hist)
    #spectrumV20.automatic_coloring(30)
    #spectrumV16.automatic_coloring(30)
    #spectrumV24.automatic_coloring(30)
    
    
    
    #spectrumV24.automatic_coloring(30)
    
    Ls = np.linspace(14,26)
    channel = {'Charm':0,'Strange':0,'Isospin':0,'Charm_Isospin':1,'C_parity':-1}
    spectra = [spectrumV16,spectrumV24,spectrumV20]
    spectrum_analysis.Spectrum.plot_irrep_mass_with_get_finite(spectra,"T1mM",Ls,"000",14,channel,0.73,"./Data/output_d000_T1m.spectrum")
    included = [[2,7,4,3],[3,4,2,5,6,7],[2,7,3,4]]
    #spectrum_analysis.Spectrum.plot_irrep_mass_no_sigma_grey_out(spectra,"T1mM",Ls,"000",14,channel,0.71,included)

    
def create_E_file_no_sigma_fewer():
    mom = "000"
    irrep = "T1mM-nosigma-fewer"
    figures_save, create_files,save_hist = False,False,False
    spectrumV24 = spectrum_analysis.Spectrum(24,10,irrep,create_files,False,False)

    spectrumV20 = spectrum_analysis.Spectrum(20,11,irrep,create_files,figures_save,save_hist)
    spectrumV16 = spectrum_analysis.Spectrum(16,11,irrep,create_files,figures_save,save_hist)
    
    prefix = "/CETQCD3/cet34/HadSpec"
    prefix_volume_16 =  "/szscl21_16_128_b1p50_t_x4p300_um0p0840_sm0p0743_n1p265_per"
    prefix_volume_20 = "/szscl21_20_256_b1p50_t_x4p300_um0p0840_sm0p0743_n1p265_per"
    prefix_volume_24 = "/szscl21_24_128_b1p50_t_x4p300_um0p0840_sm0p0743_n1p265"
    prefix_redstar = "/redstar/DDbar_I0/fits_pm"
    irrep_prefix = "/p000-T1mM-nosigma-fewer"
    path_24 = prefix + prefix_volume_24 + prefix_redstar +irrep_prefix
    path_20 = prefix + prefix_volume_20 + prefix_redstar + irrep_prefix
    path_16 = prefix + prefix_volume_16 + prefix_redstar + irrep_prefix
    spectra_per_irrep_per_volume = {"T1m":{16:spectrumV16,20:spectrumV20,24:spectrumV24}}
    paths_per_irrep_per_volume =  {"T1m":{16:path_16,20:path_20,24:path_24}}
    excluded_min_16 = 7
    excluded_min_20 = 7
    excluded_min_24 = 7
    excluded_levels_16 = []
    excluded_levels_20 = []
    excluded_levels_24 = []
    for j in range(60):
        if j>excluded_min_16:
            excluded_levels_16.append(j)
    for j in range(60):
        if j>excluded_min_20:
            excluded_levels_20.append(j)
    for j in range(60):
        if j>excluded_min_24:
            excluded_levels_24.append(j)

    excluded_levels_per_irrep_per_volume =  {"T1m":{16:excluded_levels_16,20:excluded_levels_20,24:excluded_levels_24}}


    Set_Up_Functions.SpectrumInfo(spectra_per_irrep_per_volume,[16,20,24],["T1m"],excluded_levels_per_irrep_per_volume,0.71,paths_per_irrep_per_volume)

def create_E_file_no_sigma():
    mom = "000"
    irrep = "T1mM-nosigma-tmax25"
    create_structure(24,[8,9,10,11,12],irrep)
    create_structure(16,[8,9,10,11,12],irrep)
    create_structure(20,[8,9,10,11,12,13],irrep)
    figures_save, create_files,save_hist = False,False,False
    spectrumV24t011 = spectrum_analysis.Spectrum(24,12,irrep,create_files,figures_save,save_hist)
    spectrumV20t010 = spectrum_analysis.Spectrum(20,9,irrep,create_files,figures_save,save_hist)
    spectrumV16t011 = spectrum_analysis.Spectrum(16,9,irrep,create_files,figures_save,save_hist)
    prefix = "/CETQCD3/cet34/HadSpec"
    prefix_volume_16 =  "/szscl21_16_128_b1p50_t_x4p300_um0p0840_sm0p0743_n1p265_per"
    prefix_volume_20 = "/szscl21_20_256_b1p50_t_x4p300_um0p0840_sm0p0743_n1p265_per"
    prefix_volume_24 = "/szscl21_24_128_b1p50_t_x4p300_um0p0840_sm0p0743_n1p265"
    prefix_redstar = "/redstar/DDbar_I0/fits_pm"
    irrep_prefix = "/p000-T1mM-nosigma"
    #path_24 = "/t012/MassJackFiles"
    #path_20 = "/t09/MassJackFiles"
    #path_16 = "/t09/MassJackFiles"
    path_24 = prefix + prefix_volume_24 + prefix_redstar +irrep_prefix
    path_20 = prefix + prefix_volume_20 + prefix_redstar + irrep_prefix
    path_16 = prefix + prefix_volume_16 + prefix_redstar + irrep_prefix
    spectra_per_irrep_per_volume = {"T1m":{16:spectrumV16t011,20:spectrumV20t010,24:spectrumV24t011}}
    paths_per_irrep_per_volume =  {"T1m":{16:path_16,20:path_20,24:path_24}}
    excluded_levels = []
    excluded_min_16 = 8
    excluded_min_20 = 5
    excluded_min_24 = 8
    excluded_levels_16 = []
    excluded_levels_20 = []
    excluded_levels_24 = []
    for j in range(60):
        if j>excluded_min_16:
            excluded_levels_16.append(j)
    for j in range(60):
        if j>excluded_min_20:
            excluded_levels_20.append(j)
    for j in range(60):
        if j>excluded_min_24:
            excluded_levels_24.append(j)

    excluded_levels_per_irrep_per_volume =  {"T1m":{16:excluded_levels_16+[3,4],20:excluded_levels_20+[4,3],24:excluded_levels_24+[7,5,3,8,6]}}


    Set_Up_Functions.SpectrumInfo(spectra_per_irrep_per_volume,[16,20,24],["T1m"],excluded_levels_per_irrep_per_volume,0.71,paths_per_irrep_per_volume)


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
    #spectrum_analysis.Spectrum.plot_irrep_mass(spectra,"T1pM",Ls,"000",30,channel,0.74)
    spectrum_analysis.Spectrum.plot_irrep_mass_with_get_finite(spectra,"T1pM",Ls,"000",30,channel,0.74,"./Data/output_d000_T1m.spectrum ")
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
    #colors.create_table_operators("T1mM-nosigma-tmax25",24)
    #colors.create_table_operators("T1mM-nosigma-tmax25",20)
    #colors.create_table_operators("T1mM-nosigma-tmax25",16)
    Ls = np.linspace(14,26)
    Ps = [np.array([0,0,0]),np.array([1,0,0]),np.array([1,1,0]),np.array([1,1,1])]
    channel = {'Charm':0,'Strange':0,'Isospin':0,'Charm_Isospin':1,'C_parity':-1}
    spectrumV24t08full = spectrum_analysis.Spectrum(24,8,"T1mM",create_files,figures_save,save_hist)
    spectrumV20t09full = spectrum_analysis.Spectrum(20,9,"T1mM",create_files,figures_save,save_hist)
    spectrumV16t011full = spectrum_analysis.Spectrum(16,11,"T1mM",create_files,figures_save,save_hist)
    #[spectrumV24t011.plot_histogram_of_levels_of_color('grey',i) for i in range(10)]
    spectra = [spectrumV16t011,spectrumV24t011,spectrumV20t010]
    spectrum_analysis.Spectrum.plot_irrep_mass_no_sigma(spectra,"T1mM",Ls,"000",25,channel,0.74)

    #spectrum_analysis.Spectrum.plot_irrep_mass_with_get_finite(spectra,"T1mM",Ls,"000",25,channel,0.74,"./Data/output_d000_T1m.spectrum ")
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
    #spectrumV24t08.plot_nine_levels([0,1,2,3,4,5,6,7,8],"$T_1^{--}$",0.9,1.2,0,30)
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

import Set_Up_Functions as SUF


def main():
    #xml_functions.create_hadron_xml("./Data/Hadrons.xml","./Data/Hadrons.dat")
    #create_E_file_no_sigma_fewer()
    #rig,ax = plt.subplots()
    #ps.plot_elastic_phase_shift_discrete("./Data/discrete.data",ax)
    #labels = ["$D\\bar{D}\\to D\\bar{D}$"]
    #paths = ['Amplitudes_plot_data/phase_shift.dat']
    #pa.plot_amplitudes_phase_shift(paths,labels,1,-1,ax)
    #ax.set_xlim(0.64,0.73)
    #plt.show()
    #spectrum_analysis.bias_analysis_create_structure("T1mM-nosigma-reduced",20,[25,26,27,28,29,30],10,14)
    #spectrum_analysis.bias_analysis_load_and_create_plots("T1mM-nosigma-improved",20,[27,28,29],10,14)
    #spectrum_analysis.Spectrum(20,10,"T1mM-nosigma-reduced",False,False,False).automatic_coloring(30)
    markers = ['o','s','^','v','<','>','D','p','*']
    #spectrum_analysis.Spectrum(20,10,"T1mM-nosigma-improved_tmax_30",False,False,False).automatic_coloring(16)
    #spectrum_analysis.plot_nine_levels_bias_analysis_ref_irrep("T1mM-nosigma-reduced",["T1mM-nosigma-reduced","T1mM-nosigma-improved"],20,[27,28,29,30],10,11,30,10,30,[1,2,3,13,15,6,7,8,9],markers)
    #plt.show()
    #preliminary_analysis_T1mM_A2mM_A1()
    #preliminary_analysis_T1mM_nosigma_reduced()
    #preliminary_analysis_T1mM()
    preliminary_analysis_A2mM_nosigma_reduced()
    #preliminar_analysis_T1mM_fewer_pm()
    #preliminary_analysis_D4A1M_nosigma_2()
    #preliminary_analysis_T1mM_nosigma_reduced_3()
    #preliminary_analysis_T1mM_A2mM_A1()
    #preliminary_analysis_T1pM_corrected()
    #preliminary_analysis_T1pM_nosigma()
    #preliminary_analysis_T1pM_corrected()
    #xml_functions.create_splines_xml("./Data/make_splines.xml","./Data/PairsOfHadrons.dat","./Data/LatticeParameters.dat",4)
    #preliminary_analysis_T1mM_nosigma_fewer()
    th = 0.74
    channel = {'Charm':0,'Strange':0,'Isospin':0,'Charm_Isospin':1,'C_parity':-1}
    #SUF.pairsofHadrons(channel,[0,1],th)
    irreps_rest = ["T1mM","T1pM","T2mM","T2pM","A1mM","A1pM","A2mM","A2pM","EpM","EmM"]
    irreps_unit_mom = ["A1M","A2M","B1M","B2M","E2M"]
    #no_int_E_levels.purge_sym_no_int_levels(channel)
    #[create_irrep_ni_minus_C_parity(irrep_rest,"000",th) for irrep_rest in irreps_rest]
    #[create_irrep_ni_minus_C_parity(irrep_unit_mom,"100",th) for irrep_unit_mom in irreps_unit_mom]

    th = 0.74
    channelp = {'Charm':0,'Strange':0,'Isospin':0,'Charm_Isospin':1,'C_parity':1}
    irreps_rest = ["T1mP","T1pP","T2mP","T2pP","A1mP","A1pP","A2mP","A2pP","EpP","EmP"]
    irreps_unit_mom = ["A1P","A2P","B1P","B2P","E2P"]
    #no_int_E_levels.purge_sym_no_int_levels(channelp)
    #[create_irrep_ni_plus_C_parity(irrep_rest,"000",th) for irrep_rest in irreps_rest]
    #[create_irrep_ni_plus_C_parity(irrep_unit_mom,"100",th) for irrep_unit_mom in irreps_unit_mom]

    #no_int_E_levels.purge_sym_no_int_levels(channel)
    #no_int_E_levels.create_operator_list(channel,"T1mM","000",0.74)
    
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
