def create_amplitude_xml_tsq(pair_one,pair_two,path_fit,path_amplitude,twoJ,P,Ecm_min,Ecm_max):
    with open(path_amplitude, 'w') as f:
        f.write('<?xml version="1.0"?>\n')
        f.write('<PlotControl>\n')
        f.write('  <fit_minimum_filename>'+path_fit+'</fit_minimum_filename>\n')
        f.write('  <plot>plot_amplitude</plot>\n')
        f.write('  <nbins>100</nbins>\n')
        f.write('  <twoJ>'+str(twoJ)+'</twoJ>\n')
        f.write('  <P>'+str(P)+'</P>\n')
        f.write('  <Ecm_min>'+str(Ecm_min)+'</Ecm_min>\n')
        f.write('  <Ecm_max>'+str(Ecm_max)+'</Ecm_max>\n')
        f.write('  <npts>1000</npts>\n')
        f.write('  <Channels>\n')
        f.write('    <elem>\n')
        hadron_one_pair_one = pair_one['hardon_one']
        hadron_two_pair_one = pair_one['hardon_two']
        twoS = pair_one['twoS']
        ell = pair_one['ell']
        f.write('      <hadron1>'+hadron_one_pair_one+'</hadron1>\n')
        f.write('      <hadron2>'+hadron_two_pair_one+'</hadron2>\n')
        f.write('      <twoS>'+str(twoS)+'</twoS>\n')
        f.write('      <ell>'+str(ell)+'</ell>\n')
        f.write('      <twoJ>'+str(twoJ)+'</twoJ>\n')
        f.write('    </elem>\n')
        hadron_one_pair_two = pair_two['hardon_one']
        hadron_two_pair_two = pair_two['hardon_two']
        twoS = pair_two['twoS']
        ell = pair_two['ell']
        f.write('    <elem>\n')
        f.write('      <hadron1>'+hadron_one_pair_two+'</hadron1>\n')
        f.write('      <hadron2>'+hadron_two_pair_two+'</hadron2>\n')
        f.write('      <twoS>'+str(twoS)+'</twoS>\n')
        f.write('      <ell>'+str(ell)+'</ell>\n')
        f.write('      <twoJ>'+str(twoJ)+'</twoJ>\n')
        f.write('    </elem>\n')
        f.write('  </Channels>\n')
        f.write('</PlotControl>\n')
def read_PartwialWaves(path_read):
    partial_waves = {}
    with open(path_read, 'r') as r:
        for i,line in enumerate(r.readlines()):
            if line[0] == '#':
                string = line.split()[1]
                Jstring = ''
                Pstring = ''
                for let in string:
                    if let.isdigit() or let == ".":
                        Jstring += let
                    else:
                        Pstring += let
                twoJ = int(2*float(Jstring))
                P = 1 if Pstring == '+' else -1
                partial_waves[(twoJ,P)] = {}
            else:
                string = line.split()
                name_1 = string[0]
                name_2 = string[1]
                waves = string[2].replace("[","").replace("]","").replace("'","").replace("{","").replace("}","").split(",")
                
                for i,wave in enumerate(waves):
                    upper = False
                    spin = ""
                    other = ""
                    for let in wave:
                        if let == "^":
                            continue
                        if let.isalpha():
                            upper = True
                        if not upper:
                            spin += let
                        if upper:
                            other += let
                    waves[i] = spin+"^"+other
                partial_waves[(twoJ,P)][(name_1,name_2)] = waves
    return partial_waves
def create_combinations_partial_wave(partial_waves_given_JP):
    possible_combinations = []
    total_single_channels = []
    for channel in partial_waves_given_JP.keys():
        for partial_wave in partial_waves_given_JP[channel]:
            name_1 = channel[0]
            name_2 = channel[1]
            chan = name_1 + ":" + name_2 + "/" + partial_wave
            total_single_channels.append(chan)
    for i in range(len(total_single_channels)):
        for j in range(i,len(total_single_channels)):
            

            possible_combinations.append((total_single_channels[j]+"|"+total_single_channels[i]))


    print(possible_combinations)

        
    return total_single_channels,possible_combinations

def iterate_over_pairs(path_fit,path_partial_waves,Ecm_min=0.66,Ecm_max=0.71):

    
    dict_waves = read_PartwialWaves(path_partial_waves)
    names = {0:'S',1:'P',2:'D',3:'F',4:'G',5:'H',6:'I',7:'J',8:'K',9:'L',10:'M',11:'N',12:'O'}
    max_l =4
    reversed_names = {v: k for k, v in names.items()}
    for JP in dict_waves.keys():
        two_J = JP[0]
        P = JP[1]
        _,combinations = create_combinations_partial_wave(dict_waves[JP])
        for combi in combinations:
            elems = combi.split("|")
            pair_one = elems[0].split("/")[0]
            pair_two = elems[1].split("/")[0]
            pair_one_hadron_name_one = pair_one.split(":")[0]
            pair_one_hadron_name_two = pair_one.split(":")[1]
            ell_pair_one = reversed_names[elems[0].split("/")[1].split("^")[1].split("_")[0]]
            twoS_pair_one = int(elems[0].split("/")[1].split("^")[0])-1
            pair_two_hadron_name_one = pair_two.split(":")[0]
            pair_two_hadron_name_two = pair_two.split(":")[1]
            ell_pair_two = reversed_names[elems[1].split("/")[1].split("^")[1].split("_")[0]]
            twoS_pair_two = int(elems[1].split("/")[1].split("^")[0])-1
            p_one_dict = {'hardon_one':pair_one_hadron_name_one,'hardon_two':pair_one_hadron_name_two,'ell':ell_pair_one,'twoS':twoS_pair_one}
            p_two_dict = {'hardon_one':pair_two_hadron_name_one,'hardon_two':pair_two_hadron_name_two,'ell':ell_pair_two,'twoS':twoS_pair_two}
            combo = combi.replace(":","").replace("|","--").replace("/","-").replace("\\","")
            path_amplitude = "./amplitude/Amplitude_"+combo+".xml"


            create_amplitude_xml_tsq(p_one_dict,p_two_dict,path_fit,path_amplitude,two_J,P,Ecm_min,Ecm_max)

def main(args):
    path_fit = "./fit/out.xml"
    path_waves = "./Data/PartialWaves.dat"
    iterate_over_pairs(path_fit,path_waves)
