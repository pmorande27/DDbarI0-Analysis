from Set_Up_Functions import *
import spectrum_analysis
def create_E_file_no_sigma_fewer_coupled_1_3_reduced_A2():
    mom = "000"
    irrep = "T1mM-nosigma-reduced"
    irrep20 = "T1mM-nosigma-reduced-3"
    figures_save, create_files,save_hist = False,False,False
    spectrumV24 = spectrum_analysis.Spectrum(24,12,irrep,create_files,figures_save,save_hist)

    spectrumV20 = spectrum_analysis.Spectrum(20,11,irrep20,create_files,figures_save,save_hist)
    spectrumV16 = spectrum_analysis.Spectrum(16,12,irrep,create_files,figures_save,save_hist)
    irrep = "A2mM-nosigma-reduced"
    spectrumV16_A2mM  = spectrum_analysis.Spectrum(16,12,irrep,create_files,figures_save,save_hist)
    spectrumV20_A2mM = spectrum_analysis.Spectrum(20,14,irrep,create_files,figures_save,save_hist)
    spectrumV24_A2mM = spectrum_analysis.Spectrum(24,14,irrep,create_files,figures_save,save_hist)

    
    prefix = "/CETQCD3/cet34/HadSpec"
    prefix_volume_16 =  "/szscl21_16_128_b1p50_t_x4p300_um0p0840_sm0p0743_n1p265_per"
    prefix_volume_20 = "/szscl21_20_256_b1p50_t_x4p300_um0p0840_sm0p0743_n1p265_per"
    prefix_volume_24 = "/szscl21_24_128_b1p50_t_x4p300_um0p0840_sm0p0743_n1p265"
    prefix_redstar = "/redstar/DDbar_I0/fits_pm"
    irrep_prefix = "/p000-T1mM-nosigma-reduced"
    irrep_prefix_20 = "/p000-T1mM-nosigma-reduced-2"
    irrep_prefix_A2 = "/p000-A2mM-nosigma-reduced"
    irrep_prefix_A2_16 = "/p000-A2mM-nosigma-reduced"
    path_24 = prefix + prefix_volume_24 + prefix_redstar +irrep_prefix
    path_20 = prefix + prefix_volume_20 + prefix_redstar + irrep_prefix_20
    path_16 = prefix + prefix_volume_16 + prefix_redstar + irrep_prefix
    path_16_A2 = prefix + prefix_volume_16 + prefix_redstar + irrep_prefix_A2_16
    path_20_A2 = prefix + prefix_volume_20 + prefix_redstar + irrep_prefix_A2
    path_24_A2 = prefix + prefix_volume_24 + prefix_redstar + irrep_prefix_A2
    spectra_per_irrep_per_volume = {"T1m":{16:spectrumV16,20:spectrumV20,24:spectrumV24},"A2m":{16:spectrumV16_A2mM,20:spectrumV20_A2mM,24:spectrumV24_A2mM}}
    paths_per_irrep_per_volume =  {"T1m":{16:path_16,20:path_20,24:path_24},"A2m":{16:path_16_A2,20:path_20_A2,24:path_24_A2}}
    excluded_min_16 = 6
    excluded_min_20 = 10
    excluded_min_24 = 11
    excluded_levels_16 = [5]
    excluded_levels_20 = [6,10,7]
    excluded_levels_24 = [7,9,10]
    excluded_levels_16_A2 = [5]
    excluded_levels_20_A2 = [5,6,9,10,11,12,14]
    excluded_levels_24_A2 = [7,9,10]
    excluded_min_16_A2 = 1
    excluded_min_20_A2 = 0
    excluded_min_24_A2 = 2
    for j in range(60):
        if j>excluded_min_16:
            excluded_levels_16.append(j)
    for j in range(60):
        if j>excluded_min_20:
            excluded_levels_20.append(j)
    for j in range(60):
        if j>excluded_min_24:
            excluded_levels_24.append(j)
    for j in range(60):
        if j>excluded_min_16_A2:
            excluded_levels_16_A2.append(j)
    for j in range(60):
        if j>excluded_min_20_A2:
            excluded_levels_20_A2.append(j)
    for j in range(60):
        if j>excluded_min_24_A2:
            excluded_levels_24_A2.append(j)
    excluded_levels_per_irrep_per_volume =  {"T1m":{16:excluded_levels_16,20:excluded_levels_20,24:excluded_levels_24},"A2m":{16:excluded_levels_16_A2,20:excluded_levels_20_A2,24:excluded_levels_24_A2}}


    SpectrumInfo(spectra_per_irrep_per_volume,[16,20,24],["T1m","A2m"],excluded_levels_per_irrep_per_volume,0.76,paths_per_irrep_per_volume)
def create_E_file_no_sigma_fewer_coupled_1_3_reduced():
    mom = "000"
    irrep = "T1mM-nosigma-reduced"
    figures_save, create_files,save_hist = False,False,False
    spectrumV24 = spectrum_analysis.Spectrum(24,12,irrep,create_files,figures_save,save_hist)

    spectrumV20 = spectrum_analysis.Spectrum(20,10,irrep,create_files,figures_save,save_hist)
    spectrumV16 = spectrum_analysis.Spectrum(16,12,irrep,create_files,figures_save,save_hist)
    
    prefix = "/CETQCD3/cet34/HadSpec"
    prefix_volume_16 =  "/szscl21_16_128_b1p50_t_x4p300_um0p0840_sm0p0743_n1p265_per"
    prefix_volume_20 = "/szscl21_20_256_b1p50_t_x4p300_um0p0840_sm0p0743_n1p265_per"
    prefix_volume_24 = "/szscl21_24_128_b1p50_t_x4p300_um0p0840_sm0p0743_n1p265"
    prefix_redstar = "/redstar/DDbar_I0/fits_pm"
    irrep_prefix = "/p000-T1mM-nosigma-reduced"
    path_24 = prefix + prefix_volume_24 + prefix_redstar +irrep_prefix
    path_20 = prefix + prefix_volume_20 + prefix_redstar + irrep_prefix
    path_16 = prefix + prefix_volume_16 + prefix_redstar + irrep_prefix
    spectra_per_irrep_per_volume = {"T1m":{16:spectrumV16,20:spectrumV20,24:spectrumV24}}
    paths_per_irrep_per_volume =  {"T1m":{16:path_16,20:path_20,24:path_24}}
    excluded_min_16 = 6
    excluded_min_20 = 15
    excluded_min_24 = 11
    excluded_levels_16 = [5]
    excluded_levels_20 = [5,6,9,10,11,12,14]
    excluded_levels_24 = [7,9,10]
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


    SpectrumInfo(spectra_per_irrep_per_volume,[16,20,24],["T1m"],excluded_levels_per_irrep_per_volume,0.72,paths_per_irrep_per_volume)
def create_E_file_no_sigma_fewer_coupled_1_3():
    mom = "000"
    irrep = "T1mM-nosigma-improved"
    figures_save, create_files,save_hist = False,False,False
    spectrumV24 = spectrum_analysis.Spectrum(24,12,irrep,create_files,figures_save,save_hist)

    spectrumV20 = spectrum_analysis.Spectrum(20,12,irrep,create_files,figures_save,save_hist)
    spectrumV16 = spectrum_analysis.Spectrum(16,10,irrep,create_files,figures_save,save_hist)
    
    prefix = "/CETQCD3/cet34/HadSpec"
    prefix_volume_16 =  "/szscl21_16_128_b1p50_t_x4p300_um0p0840_sm0p0743_n1p265_per"
    prefix_volume_20 = "/szscl21_20_256_b1p50_t_x4p300_um0p0840_sm0p0743_n1p265_per"
    prefix_volume_24 = "/szscl21_24_128_b1p50_t_x4p300_um0p0840_sm0p0743_n1p265"
    prefix_redstar = "/redstar/DDbar_I0/fits_pm"
    irrep_prefix = "/p000-T1mM-nosigma-improved"
    path_24 = prefix + prefix_volume_24 + prefix_redstar +irrep_prefix
    path_20 = prefix + prefix_volume_20 + prefix_redstar + irrep_prefix
    path_16 = prefix + prefix_volume_16 + prefix_redstar + irrep_prefix
    spectra_per_irrep_per_volume = {"T1m":{16:spectrumV16,20:spectrumV20,24:spectrumV24}}
    paths_per_irrep_per_volume =  {"T1m":{16:path_16,20:path_20,24:path_24}}
    excluded_min_16 = 9
    excluded_min_20 = 5
    excluded_min_24 = 11
    excluded_levels_16 = [5,6,7]
    excluded_levels_20 = []
    excluded_levels_24 = [7,9,10]
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


    SpectrumInfo(spectra_per_irrep_per_volume,[16,20,24],["T1m"],excluded_levels_per_irrep_per_volume,0.71,paths_per_irrep_per_volume)
def create_E_file_no_sigma_improved_coupled_1():
    mom = "000"
    irrep = "T1mM-nosigma-improved"
    figures_save, create_files,save_hist = False,False,False
    spectrumV24 = spectrum_analysis.Spectrum(24,12,irrep,create_files,figures_save,save_hist)

    spectrumV20 = spectrum_analysis.Spectrum(20,12,irrep,create_files,figures_save,save_hist)
    spectrumV16 = spectrum_analysis.Spectrum(16,10,irrep,create_files,figures_save,save_hist)
    
    prefix = "/CETQCD3/cet34/HadSpec"
    prefix_volume_16 =  "/szscl21_16_128_b1p50_t_x4p300_um0p0840_sm0p0743_n1p265_per"
    prefix_volume_20 = "/szscl21_20_256_b1p50_t_x4p300_um0p0840_sm0p0743_n1p265_per"
    prefix_volume_24 = "/szscl21_24_128_b1p50_t_x4p300_um0p0840_sm0p0743_n1p265"
    prefix_redstar = "/redstar/DDbar_I0/fits_pm"
    irrep_prefix = "/p000-T1mM-nosigma-improved"
    path_24 = prefix + prefix_volume_24 + prefix_redstar +irrep_prefix
    path_20 = prefix + prefix_volume_20 + prefix_redstar + irrep_prefix
    path_16 = prefix + prefix_volume_16 + prefix_redstar + irrep_prefix
    spectra_per_irrep_per_volume = {"T1m":{16:spectrumV16,20:spectrumV20,24:spectrumV24}}
    paths_per_irrep_per_volume =  {"T1m":{16:path_16,20:path_20,24:path_24}}
    excluded_min_16 = 8
    excluded_min_20 = 7
    excluded_min_24 = 7
    excluded_levels_16 = [3]
    excluded_levels_20 = [4]
    excluded_levels_24 = [5,6,8,7]
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


    SpectrumInfo(spectra_per_irrep_per_volume,[16,20,24],["T1m"],excluded_levels_per_irrep_per_volume,0.71,paths_per_irrep_per_volume)

def create_E_file_no_sigma_fewer_coupled_1():
    mom = "000"
    irrep = "T1mM-nosigma-fewer"
    figures_save, create_files,save_hist = False,False,False
    spectrumV24 = spectrum_analysis.Spectrum(24,10,irrep,create_files,figures_save,save_hist)

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
    excluded_levels_16 = [3]
    excluded_levels_20 = [4]
    excluded_levels_24 = [5,6,7]
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


    SpectrumInfo(spectra_per_irrep_per_volume,[16,20,24],["T1m"],excluded_levels_per_irrep_per_volume,0.71,paths_per_irrep_per_volume)

def create_E_file_no_sigma_fewer_elastic():
    mom = "000"
    irrep = "T1mM-nosigma-fewer"
    figures_save, create_files,save_hist = False,False,False
    spectrumV24 = spectrum_analysis.Spectrum(24,10,irrep,create_files,figures_save,save_hist)

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
    excluded_levels_16 = [4,3]
    excluded_levels_20 = [3,4]
    excluded_levels_24 = [2,7,6,5]
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


    SpectrumInfo(spectra_per_irrep_per_volume,[16,20,24],["T1m"],excluded_levels_per_irrep_per_volume,0.71,paths_per_irrep_per_volume)

def data_nosigma_fewer_elastic_simple():
    channel = {'Charm':0,'Strange':0,'Isospin':0,'Charm_Isospin':1,'C_parity':-1}
    target_pairs = "./Data/PairsOfHadrons.dat"
    target_JS = ["1.0-"]
    create_E_file_no_sigma_fewer_elastic()
    Hadrons()
    pairsofHadrons(channel,[1],0.74)
    set_up_partial_waves(channel,1,0.74,target_JS,target_pairs)
    
    model_parameters({"1.0-":{'n_poles':2,'poly_order':-1,"chew_man":"true"
    ,"chew_man_sub_point": "pole0"}})
    #set_up_initial_values("./Data/ModelParameters.dat")
    #copy initial values to the right place
    """
    ('m_pole0_1.0_-', '0.64566')
    ('m_pole0_err_1.0_-', '0.00111')
    ('m_pole0_fix_1.0_-', 'false')
    ('m_pole0_limits_1.0_-', '0.64 0.65')
    ('g_D:\\bar{D}/1^P_1_pole0', '0.39627')
    ('g_D:\\bar{D}/1^P_1_pole0_err', '0.03')
    ('g_D:\\bar{D}/1^P_1_pole0_fix', 'false')
    ('m_pole1_1.0_-', '0.671994')
    ('m_pole1_err', '0.000596')
    ('m_pole1_fix', 'false')
    ('m_pole1_limits', '0.665 0.673')
    ('g_D:\\bar{D}/1^P_1_pole1', '1.49362')
    ('g_D:\\bar{D}/1^P_1_pole1_err', '0.149')
    ('g_D:\\bar{D}/1^P_1_pole1_fix', 'false')
    """
def data_nosigma_fewer_elastic_simple_one_pole():
    channel = {'Charm':0,'Strange':0,'Isospin':0,'Charm_Isospin':1,'C_parity':-1}
    target_pairs = "./Data/PairsOfHadrons.dat"
    target_JS = ["1.0-"]
    create_E_file_no_sigma_fewer_elastic()
    Hadrons()
    pairsofHadrons(channel,[1],0.74)
    set_up_partial_waves(channel,1,0.74,target_JS,target_pairs)
    
    model_parameters({"1.0-":{'n_poles':1,'poly_order':0,"chew_man":"true"
    ,"chew_man_sub_point": "pole0"}})
    #set_up_initial_values("./Data/ModelParameters.dat")
    #copy initial values to the right place
    """
    ('m_pole0_1.0_-', '0.64566')
    ('m_pole0_err_1.0_-', '0.00111')
    ('m_pole0_fix_1.0_-', 'false')
    ('m_pole0_limits_1.0_-', '0.64 0.65')
    ('g_D:\\bar{D}/1^P_1_pole0', '0.39627')
    ('g_D:\\bar{D}/1^P_1_pole0_err', '0.03')
    ('g_D:\\bar{D}/1^P_1_pole0_fix', 'false')
    ('m_pole1_1.0_-', '0.671994')
    ('m_pole1_err', '0.000596')
    ('m_pole1_fix', 'false')
    ('m_pole1_limits', '0.665 0.673')
    ('g_D:\\bar{D}/1^P_1_pole1', '1.49362')
    ('g_D:\\bar{D}/1^P_1_pole1_err', '0.149')
    ('g_D:\\bar{D}/1^P_1_pole1_fix', 'false')
    ('gamma_D:\\bar{D}/1^P_1|D:\\bar{D}/1^P_1_order0', '37.452')
    ('gamma_D:\\bar{D}/1^P_1|D:\\bar{D}/1^P_1_order0_err', '20.412')
    ('gamma_D:\\bar{D}/1^P_1|D:\\bar{D}/1^P_1_order0_fix', 'false')
    """
def data_nosigma_fewer_elastic_simple_two_poles_poly0():
    channel = {'Charm':0,'Strange':0,'Isospin':0,'Charm_Isospin':1,'C_parity':-1}
    target_pairs = "./Data/PairsOfHadrons.dat"
    target_JS = ["1.0-"]
    create_E_file_no_sigma_fewer_elastic()
    Hadrons()
    pairsofHadrons(channel,[1],0.74)
    set_up_partial_waves(channel,1,0.74,target_JS,target_pairs)
    
    model_parameters({"1.0-":{'n_poles':2,'poly_order':0,"chew_man":"true"
    ,"chew_man_sub_point": "pole0"}})
    #set_up_initial_values("./Data/ModelParameters.dat")
    #copy initial values to the right place
    """
    ('m_pole0_1.0_-', '0.671994')
    ('m_pole0_err_1.0_-', '0.000596')
    ('m_pole0_fix_1.0_-', 'false')
    ('m_pole0_limits_1.0_-', '0.665 0.673')
    ('g_D:\\bar{D}/1^P_1_pole0', '1.49362')
    ('g_D:\\bar{D}/1^P_1_pole0_err', '0.149')
    ('g_D:\\bar{D}/1^P_1_pole0_fix', 'false')
    ('gamma_D:\\bar{D}/1^P_1|D:\\bar{D}/1^P_1_order0', '37.452')
    ('gamma_D:\\bar{D}/1^P_1|D:\\bar{D}/1^P_1_order0_err', '20.412')
    ('gamma_D:\\bar{D}/1^P_1|D:\\bar{D}/1^P_1_order0_fix', 'false')
    """
def data_nosigma_fewer_coupled_1_one_pole_poly_0():
    channel = {'Charm':0,'Strange':0,'Isospin':0,'Charm_Isospin':1,'C_parity':-1}
    target_pairs = "./Data/PairsOfHadrons.dat"
    target_JS = ["1.0-"]
    create_E_file_no_sigma_fewer_coupled_1()
    Hadrons()
    pairsofHadrons(channel,[0,1],0.74)
    set_up_partial_waves(channel,1,0.74,target_JS,target_pairs)
    
    model_parameters({"1.0-":{'n_poles':1,'poly_order':0,"chew_man":"true"
    ,"chew_man_sub_point": "pole0"}})
    #set_up_initial_values("./Data/ModelParameters.dat")
    #copy initial values to the right place
    """
    ('m_pole0_1.0_-', '0.66813')
    ('m_pole0_err_1.0_-', '0.00073')
    ('m_pole0_fix_1.0_-', 'false')
    ('m_pole0_limits_1.0_-', '0.667 0.673')
    ('g_\\psi:\\eta/3^P_1_pole0', '0.066641')
    ('g_\\psi:\\eta/3^P_1_pole0_err', '0.17336')
    ('g_\\psi:\\eta/3^P_1_pole0_fix', 'false')
    ('g_D:\\bar{D}/1^P_1_pole0', '1.7639')
    ('g_D:\\bar{D}/1^P_1_pole0_err', '0.19549')
    ('g_D:\\bar{D}/1^P_1_pole0_fix', 'false')
    ('gamma_\\psi:\\eta/3^P_1|\\psi:\\eta/3^P_1_order0', '-9.0297')
    ('gamma_\\psi:\\eta/3^P_1|\\psi:\\eta/3^P_1_order0_err', '4.995')
    ('gamma_\\psi:\\eta/3^P_1|\\psi:\\eta/3^P_1_order0_fix', 'false')
    ('gamma_D:\\bar{D}/1^P_1|\\psi:\\eta/3^P_1_order0', '-6.6344')
    ('gamma_D:\\bar{D}/1^P_1|\\psi:\\eta/3^P_1_order0_err', '7.5993')
    ('gamma_D:\\bar{D}/1^P_1|\\psi:\\eta/3^P_1_order0_fix', 'false')
    ('gamma_D:\\bar{D}/1^P_1|D:\\bar{D}/1^P_1_order0', '37.452')
    ('gamma_D:\\bar{D}/1^P_1|D:\\bar{D}/1^P_1_order0_err', '20.412')
    ('gamma_D:\\bar{D}/1^P_1|D:\\bar{D}/1^P_1_order0_fix', 'false')

    """
def data_nosigma_fewer_coupled_1_one_pole_poly_0_3_one_pole_poly_m1():
    channel = {'Charm':0,'Strange':0,'Isospin':0,'Charm_Isospin':1,'C_parity':-1}
    target_pairs = "./Data/PairsOfHadrons.dat"
    target_JS = ["1.0-","3.0-"]
    create_E_file_no_sigma_fewer_coupled_1_3()
    Hadrons()
    pairsofHadrons(channel,[0,1],0.74)
    set_up_partial_waves(channel,3,0.74,target_JS,target_pairs)
    
    model_parameters({"1.0-":{'n_poles':1,'poly_order':0,"chew_man":"true"
    ,"chew_man_sub_point": "pole0"},"3.0-":{'n_poles':1,'poly_order':-1,"chew_man":"true"
    ,"chew_man_sub_point": "pole0"}})
    set_up_initial_values("./Data/ModelParameters.dat")
    #copy initial values to the right place
    """
    ('m_pole0_1.0_-', '0.67070')
    ('m_pole0_err_1.0_-', '0.00049')
    ('m_pole0_fix_1.0_-', 'false')
    ('m_pole0_limits_1.0_-', '0.668 0.673')
    ('g_\\psi:\\eta/3^P_1_pole0', '0')
    ('g_\\psi:\\eta/3^P_1_pole0_err', '0')
    ('g_\\psi:\\eta/3^P_1_pole0_fix', 'true')
    ('g_D:\\bar{D}/1^P_1_pole0', '1.7639')
    ('g_D:\\bar{D}/1^P_1_pole0_err', '0.19549')
    ('g_D:\\bar{D}/1^P_1_pole0_fix', 'false')
    ('gamma_\\psi:\\eta/3^P_1|\\psi:\\eta/3^P_1_order0', '0.6')
    ('gamma_\\psi:\\eta/3^P_1|\\psi:\\eta/3^P_1_order0_err', '2')
    ('gamma_\\psi:\\eta/3^P_1|\\psi:\\eta/3^P_1_order0_fix', 'false')
    ('gamma_D:\\bar{D}/1^P_1|\\psi:\\eta/3^P_1_order0', '0')
    ('gamma_D:\\bar{D}/1^P_1|\\psi:\\eta/3^P_1_order0_err', '0')
    ('gamma_D:\\bar{D}/1^P_1|\\psi:\\eta/3^P_1_order0_fix', 'true')
    ('gamma_D:\\bar{D}/1^P_1|D:\\bar{D}/1^P_1_order0', '115.61')
    ('gamma_D:\\bar{D}/1^P_1|D:\\bar{D}/1^P_1_order0_err', '26.368')
    ('gamma_D:\\bar{D}/1^P_1|D:\\bar{D}/1^P_1_order0_fix', 'false')
    ('m_pole0_3.0_-', '0.67853')
    ('m_pole0_err_3.0_-', '0.00058')
    ('m_pole0_fix_3.0_-', 'false')
    ('m_pole0_limits_3.0_-', '0.672 0.682')
    ('g_\\psi:\\eta/3^F_3_pole0', '0.1')
    ('g_\\psi:\\eta/3^F_3_pole0_err', '0.2')
    ('g_\\psi:\\eta/3^F_3_pole0_fix', 'false')
    ('g_D:\\bar{D}/1^F_3_pole0', '13.816')
    ('g_D:\\bar{D}/1^F_3_pole0_err', '1.5122')
    ('g_D:\\bar{D}/1^F_3_pole0_fix', 'false')
    """
    
def data_nosigma_fewer_coupled_1_two_pole_poly_0():
    channel = {'Charm':0,'Strange':0,'Isospin':0,'Charm_Isospin':1,'C_parity':-1}
    target_pairs = "./Data/PairsOfHadrons.dat"
    target_JS = ["1.0-"]
    create_E_file_no_sigma_fewer_coupled_1()
    Hadrons()
    pairsofHadrons(channel,[0,1],0.74)
    set_up_partial_waves(channel,1,0.74,target_JS,target_pairs)
    
    model_parameters({"1.0-":{'n_poles':2,'poly_order':0,"chew_man":"true"
    ,"chew_man_sub_point": "pole0"}})
    #set_up_initial_values("./Data/ModelParameters.dat")
    #copy initial values to the right place
    """
    ('m_pole0_1.0_-', '0.64535')
    ('m_pole0_err_1.0_-', '0.0018')
    ('m_pole0_fix_1.0_-', 'false')
    ('m_pole0_limits', '0.64 0.65')
    ('g_\\psi:\\eta/3^P_1_pole0', '0.1')
    ('g_\\psi:\\eta/3^P_1_pole0_err', '0.17336')
    ('g_\\psi:\\eta/3^P_1_pole0_fix', 'false')
    ('g_D:\\bar{D}/1^P_1_pole0', '3.0119')
    ('g_D:\\bar{D}/1^P_1_pole0_err', '0.3813')
    ('g_D:\\bar{D}/1^P_1_pole0_fix', 'false')
    ('m_pole1_1.0_-', '0.66813')
    ('m_pole1_err_1.0_-', '0.00073')
    ('m_pole1_fix_1.0_-', 'false')
    ('m_pole1_limits_1.0_-', '0.667 0.673')
    ('g_\\psi:\\eta/3^P_1_pole1', '0.0')
    ('g_\\psi:\\eta/3^P_1_pole1_err', '0.0')
    ('g_\\psi:\\eta/3^P_1_pole1_fix', 'true')
    ('g_D:\\bar{D}/1^P_1_pole1', '1.7639')
    ('g_D:\\bar{D}/1^P_1_pole1_err', '0.19549')
    ('g_D:\\bar{D}/1^P_1_pole1_fix', 'false')
    ('gamma_\\psi:\\eta/3^P_1|\\psi:\\eta/3^P_1_order0', '-9.0297')
    ('gamma_\\psi:\\eta/3^P_1|\\psi:\\eta/3^P_1_order0_err', '4.995')
    ('gamma_\\psi:\\eta/3^P_1|\\psi:\\eta/3^P_1_order0_fix', 'false')
    ('gamma_D:\\bar{D}/1^P_1|\\psi:\\eta/3^P_1_order0', '0')
    ('gamma_D:\\bar{D}/1^P_1|\\psi:\\eta/3^P_1_order0_err', '0')
    ('gamma_D:\\bar{D}/1^P_1|\\psi:\\eta/3^P_1_order0_fix', 'true')
    ('gamma_D:\\bar{D}/1^P_1|D:\\bar{D}/1^P_1_order0', '37.452')
    ('gamma_D:\\bar{D}/1^P_1|D:\\bar{D}/1^P_1_order0_err', '20.412')
    ('gamma_D:\\bar{D}/1^P_1|D:\\bar{D}/1^P_1_order0_fix', 'false')

    """
def data_nosigma_fewer_coupled_1_one_pole_poly_0_3_one_pole_poly_0():
    channel = {'Charm':0,'Strange':0,'Isospin':0,'Charm_Isospin':1,'C_parity':-1}
    target_pairs = "./Data/PairsOfHadrons.dat"
    target_JS = ["1.0-","3.0-"]
    create_E_file_no_sigma_fewer_coupled_1_3()
    Hadrons()
    pairsofHadrons(channel,[0,1],0.74)
    set_up_partial_waves(channel,3,0.74,target_JS,target_pairs)
    
    model_parameters({"1.0-":{'n_poles':1,'poly_order':0,"chew_man":"true"
    ,"chew_man_sub_point": "pole0"},"3.0-":{'n_poles':1,'poly_order':0,"chew_man":"true"
    ,"chew_man_sub_point": "pole0"}})
    #set_up_initial_values("./Data/ModelParameters.dat")
    #copy initial values to the right place
    """
    ('m_pole0_1.0_-', '0.67070')
    ('m_pole0_err_1.0_-', '0.00049')
    ('m_pole0_fix_1.0_-', 'false')
    ('m_pole0_limits_1.0_-', '0.668 0.673')
    ('g_\\psi:\\eta/3^P_1_pole0', '0')
    ('g_\\psi:\\eta/3^P_1_pole0_err', '0')
    ('g_\\psi:\\eta/3^P_1_pole0_fix', 'true')
    ('g_D:\\bar{D}/1^P_1_pole0', '1.7639')
    ('g_D:\\bar{D}/1^P_1_pole0_err', '0.19549')
    ('g_D:\\bar{D}/1^P_1_pole0_fix', 'false')
    ('gamma_\\psi:\\eta/3^P_1|\\psi:\\eta/3^P_1_order0', '0.6')
    ('gamma_\\psi:\\eta/3^P_1|\\psi:\\eta/3^P_1_order0_err', '2')
    ('gamma_\\psi:\\eta/3^P_1|\\psi:\\eta/3^P_1_order0_fix', 'false')
    ('gamma_D:\\bar{D}/1^P_1|\\psi:\\eta/3^P_1_order0', '0')
    ('gamma_D:\\bar{D}/1^P_1|\\psi:\\eta/3^P_1_order0_err', '0')
    ('gamma_D:\\bar{D}/1^P_1|\\psi:\\eta/3^P_1_order0_fix', 'true')
    ('gamma_D:\\bar{D}/1^P_1|D:\\bar{D}/1^P_1_order0', '115.61')
    ('gamma_D:\\bar{D}/1^P_1|D:\\bar{D}/1^P_1_order0_err', '26.368')
    ('gamma_D:\\bar{D}/1^P_1|D:\\bar{D}/1^P_1_order0_fix', 'false')
    ('m_pole0_3.0_-', '0.67853')
    ('m_pole0_err_3.0_-', '0.00058')
    ('m_pole0_fix_3.0_-', 'false')
    ('m_pole0_limits_3.0_-', '0.672 0.682')
    ('g_\\psi:\\eta/3^F_3_pole0', '0.1')
    ('g_\\psi:\\eta/3^F_3_pole0_err', '0.2')
    ('g_\\psi:\\eta/3^F_3_pole0_fix', 'false')
    ('g_D:\\bar{D}/1^F_3_pole0', '13.816')
    ('g_D:\\bar{D}/1^F_3_pole0_err', '1.5122')
    ('g_D:\\bar{D}/1^F_3_pole0_fix', 'false')
    ('gamma_D:\\bar{D}/1^F_3|D:\\bar{D}/1^F_3_order0', '10')
    ('gamma_D:\\bar{D}/1^F_3|D:\\bar{D}/1^F_3_order0_err', '2')
    ('gamma_D:\\bar{D}/1^F_3|D:\\bar{D}/1^F_3_order0_fix', 'false')
    ('gamma_D:\\bar{D}/1^F_3|\\psi:\\eta/3^F_3_order0', '0')
    ('gamma_D:\\bar{D}/1^F_3|\\psi:\\eta/3^F_3_order0_err', '0')
    ('gamma_D:\\bar{D}/1^F_3|\\psi:\\eta/3^F_3_order0_fix', 'true')
    ('gamma_\\psi:\\eta/3^F_3|\\psi:\\eta/3^F_3_order0', '1')
    ('gamma_\\psi:\\eta/3^F_3|\\psi:\\eta/3^F_3_order0_err', '0.1')
    ('gamma_\\psi:\\eta/3^F_3|\\psi:\\eta/3^F_3_order0_fix', 'false')
    """
#data_nosigma_fewer_coupled_1_one_pole_poly_0_3_one_pole_poly_0()
#data_nosigma_fewer_elastic_simple_two_poles_poly0()
#create_E_file_no_sigma_improved_coupled_1()
#create_E_file_no_sigma_fewer_coupled_1_3()
create_E_file_no_sigma_fewer_coupled_1_3_reduced_A2()