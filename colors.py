def color_coding_dict():
    dict_nums_to_channel = {}

    dict_nums_to_channel[4] = ('D','\\bar{D}')
    dict_nums_to_channel[5] = ('D','\\bar{D}^*')
    dict_nums_to_channel[6] = ('D^*','\\bar{D}^*')
    dict_nums_to_channel[7] = ('D_s','\\bar{D}_s')
    dict_nums_to_channel[8] = ('D_s','\\bar{D}^*_s')
    dict_nums_to_channel[9] = ('D_s^*','\\bar{D}_s^*')


    dict_nums_to_channel[10] = ('\\psi','\\sigma')
    dict_nums_to_channel[11] = ('\\psi','\\eta')
    dict_nums_to_channel[12] = ('\\eta_c','w')
    dict_nums_to_channel[13] = ('\\eta_c','\\phi')
    dict_nums_to_channel[14] = ('h_c', '\\eta')
    dict_nums_to_channel[15] =  ('\\psi',"\\eta'")

    color_code_dict = {}


    color_code_dict[4] = 'green'
    color_code_dict[5] = 'orange'
    color_code_dict[6] = 'blue'
    color_code_dict[7] = 'lime'
    color_code_dict[8] = 'gold'
    color_code_dict[9] = 'cyan'

    color_code_dict[10] = 'grey'
    color_code_dict[11] = 'red'
    color_code_dict[12] = 'mediumslateblue'
    color_code_dict[13] = 'darkviolet'
    color_code_dict[14] = 'pink'
    color_code_dict[15] = 'darkred'

    color_code_dict_new = {}
    for key in color_code_dict.keys():
        color_code_dict_new[dict_nums_to_channel[key]] = color_code_dict[key]
    return color_code_dict_new
def identify_particle(particle):
    if 'Dneg_proj0' in particle:
        return 'D'
    if 'Dneg_proj1' in particle:
        return 'D^*'
    if 'etace_proj0' in particle:
        return '\\eta_c'
    if 'etace_proj1' in particle:
        return "\\eta_c'"
    if 'psice_proj0' in particle:
        return '\\psi'
    if 'psice_proj1' in particle:
        return "\\psi'"
    if 'Dsneg_proj0' in particle:
        return 'D_s'
    if 'Dsneg_proj1' in particle:
        return 'D_s^*'
    if 'eta_proj0' in particle:
        return '\\eta'
    if 'eta_proj1' in particle:
        return "\\eta'"
    if 'omega_proj0' in particle:
        return 'w'
    if 'omega_proj1' in particle:
        return "\\phi"
    if 'f_' in particle:
        return 'f_0'
    if 'hce_proj0' in particle:
        return 'h_c'
    if 'Ebarneg_proj0' in particle:
        return '\\bar{D}'
    if 'Ebarneg_proj1' in particle:
        return '\\bar{D}^*'
    if 'Esbarneg_proj0' in particle:
        return '\\bar{D}_s'
    if 'Esbarneg_proj1' in particle:
        return '\\bar{D}_s^*'
    raise ValueError(f"Particle not found p: {particle}")
def operator_identification(file_name):
    with open(file_name) as f:
        lines = f.readlines()
        operators_full_name = []
        op_short_name = []
        for i,line in enumerate(lines):
            operators_full_name.append(line)
            op_short_name.append(f"op{i}")
    dict_ops = {}
    for i, op in enumerate(op_short_name):
        if 'xx' in operators_full_name[i]:
            p1 = operators_full_name[i].split('xx')[0]
            p2 = operators_full_name[i].split('xx')[1] 
            p1 = identify_particle(p1)
            p2 = identify_particle(p2)
            channel = "$"+p1+' '+p2 + "$"
            dict_ops[op] = channel
        else:
            J = operators_full_name[i].split('__J')[1].split('_')[0]
            
            irrep = operators_full_name[i].split('_')[1].split(' ')[0]
            C = irrep[len(irrep)-1]
            P = irrep[len(irrep)-2]
            if P == 'm':
                P = '-'
            else:
                P = '+'
            if C == 'P':
                C = '+'
            else:
                C = '-'
            dict_ops[op] = '$q\\Gamma \\bar{q} ~'+J+'^{'+P+C+'}$'
    return dict_ops
def color_coding_file(filename):
    new_dict = operator_identification(filename)
    dict_nums_to_channel = {}
    dict_nums_to_channel[1] = '$q\\Gamma \\bar{q} ~1^{--}$'
    dict_nums_to_channel[2] = '$q\\Gamma \\bar{q} ~3^{--}$'
    dict_nums_to_channel[3] = '$q\\Gamma \\bar{q} ~4^{--}$'


    dict_nums_to_channel[4] = '$D \\bar{D}$'
    dict_nums_to_channel[5] = '$D \\bar{D}^*$'
    dict_nums_to_channel[6] = '$D^* \\bar{D}^*$'
    dict_nums_to_channel[7] = '$D_s \\bar{D}_s$'
    dict_nums_to_channel[8] = '$D_s \\bar{D}_s^*$'
    dict_nums_to_channel[9] = '$D_s^* \\bar{D}_s^*$'


    dict_nums_to_channel[10] = '$\\psi f_0$'
    dict_nums_to_channel[11] = '$\\psi \\eta$'
    dict_nums_to_channel[12] = '$\\eta_c w$'
    dict_nums_to_channel[13] = '$\\eta_c \\phi$'
    dict_nums_to_channel[14] = '$h_c \\eta$'
    dict_nums_to_channel[15] = "$\\psi \\eta'$"

    color_code_dict = {}
    color_code_dict[1] = 'aquamarine'
    color_code_dict[2] = 'teal'
    color_code_dict[3] = 'skyblue'

    color_code_dict[4] = 'green'
    color_code_dict[5] = 'orange'
    color_code_dict[6] = 'blue'
    color_code_dict[7] = 'lime'
    color_code_dict[8] = 'gold'
    color_code_dict[9] = 'cyan'

    color_code_dict[10] = 'grey'
    color_code_dict[11] = 'red'
    color_code_dict[12] = 'mediumslateblue'
    color_code_dict[13] = 'darkviolet'
    color_code_dict[14] = 'pink'
    color_code_dict[15] = 'darkred'

    color_code_dict_new = {}
    for key in color_code_dict.keys():
        color_code_dict_new[dict_nums_to_channel[key]] = color_code_dict[key]
    for k in new_dict.keys():
        if not new_dict[k] in color_code_dict_new.keys():
            raise ValueError(f"Key not found {new_dict[k]}, op{k} without color")
    return new_dict,color_code_dict_new

def color_coding_file_simple():
    dict_nums_to_channel = {}
    dict_nums_to_channel[1] = '$q\\Gamma \\bar{q} ~1^{--}$'
    dict_nums_to_channel[2] = '$q\\Gamma \\bar{q} ~3^{--}$'
    dict_nums_to_channel[3] = '$q\\Gamma \\bar{q} ~4^{--}$'


    dict_nums_to_channel[4] = '$D \\bar{D}$'
    dict_nums_to_channel[5] = '$D \\bar{D}^*$'
    dict_nums_to_channel[6] = '$D^* \\bar{D}^*$'
    dict_nums_to_channel[7] = '$D_s \\bar{D}_s$'
    dict_nums_to_channel[8] = '$D_s \\bar{D}_s^*$'
    dict_nums_to_channel[9] = '$D_s^* \\bar{D}_s^*$'


    dict_nums_to_channel[10] = '$\\psi f_0$'
    dict_nums_to_channel[11] = '$\\psi \\eta$'
    dict_nums_to_channel[12] = '$\\eta_c w$'
    dict_nums_to_channel[13] = '$\\eta_c \\phi$'
    dict_nums_to_channel[14] = '$h_c \\eta$'
    dict_nums_to_channel[15] = "$\\psi \\eta'$"

    color_code_dict = {}
    color_code_dict[1] = 'aquamarine'
    color_code_dict[2] = 'teal'
    color_code_dict[3] = 'skyblue'

    color_code_dict[4] = 'green'
    color_code_dict[5] = 'orange'
    color_code_dict[6] = 'blue'
    color_code_dict[7] = 'lime'
    color_code_dict[8] = 'gold'
    color_code_dict[9] = 'cyan'

    color_code_dict[10] = 'grey'
    color_code_dict[11] = 'red'
    color_code_dict[12] = 'mediumslateblue'
    color_code_dict[13] = 'darkviolet'
    color_code_dict[14] = 'pink'
    color_code_dict[15] = 'darkred'

    
    return color_code_dict

def save_color_code_state(color_code,state,spectrum):
    name_file = f"{spectrum.irrep}\\Volume_{spectrum.volume}\\t0{spectrum.t0}\\StateColorFiles\\state{state}.txt"
    with open(name_file, "w") as f:
        f.write(str(color_code))

def color_code():
    dict_operators = {}
    dict_nums_to_channel = {}
    dict_nums_to_channel[1] = '$q\\Gamma \\bar{q} ~1^{--}$'
    dict_nums_to_channel[2] = '$q\\Gamma \\bar{q} ~3^{--}$'
    dict_nums_to_channel[3] = '$D \\bar{D}$'
    dict_nums_to_channel[4] = '$D \\bar{D}^*$'
    dict_nums_to_channel[5] = '$D^* \\bar{D}^*$'
    dict_nums_to_channel[6] = '$D_s \\bar{D}_s$'
    dict_nums_to_channel[7] = '$D_s \\bar{D}^*_s$'
    dict_nums_to_channel[8] = '$D_s^* \\bar{D}_s^*$'
    dict_nums_to_channel[9] = '$\\psi f_0$'
    dict_nums_to_channel[10] = '$\\psi \\eta$'
    dict_nums_to_channel[11] = '$\\eta_c w$'
    dict_nums_to_channel[13] = '$\\eta_c \\phi$'
    dict_nums_to_channel[12] = '$q\\Gamma \\bar{q} ~4^{--}$'
    dict_nums_to_channel[14] = '$h_c \\eta$'




    for i in range(0,19):
        dict_operators[f'op{i}'] = 1
    for i in range(19,25):
        dict_operators[f'op{i}'] = 2
    dict_operators['op25'] = 12
    for i in range(26,29):
        dict_operators[f'op{i}'] = 3
    for i in range(29,31):
        dict_operators[f'op{i}'] = 6
    for i in range(31,34):
        dict_operators[f'op{i}'] = 4
    dict_operators['op34'] = 9
    dict_operators['op35'] = 10
    for i in range(36,40):
        dict_operators[f'op{i}'] = 9
    for i in range(40,42):
        dict_operators[f'op{i}'] = 10
    for i in range(42,49):
        dict_operators[f'op{i}'] = 9
    dict_operators['op49'] = 10
    dict_operators['op50'] = 11
    dict_operators['op51'] = 11
    dict_operators['op52'] = 11
    dict_operators['op53'] = 13
    dict_operators['op54'] = 14
    new_dict = {}
    keys = dict_operators.keys()
    for key in keys:
        new_dict[key] = dict_nums_to_channel[dict_operators[key]]
    color_code_dict = {}
    color_code_dict[1] = 'red'
    color_code_dict[2] = 'pink'
    color_code_dict[3] = 'green'
    color_code_dict[4] = 'lime'
    color_code_dict[5] = 'blue'
    color_code_dict[6] = 'magenta'
    color_code_dict[7] = 'violet'
    color_code_dict[8] = 'cyan'
    color_code_dict[9] = 'grey'
    color_code_dict[10] = 'coral'
    color_code_dict[11] = 'olive'
    color_code_dict[12] = 'brown'
    color_code_dict[13] = 'purple'
    color_code_dict[14] = 'teal'
    color_code_dict_new = {}
    for key in color_code_dict.keys():
        color_code_dict_new[dict_nums_to_channel[key]] = color_code_dict[key]
    
    return new_dict,color_code_dict_new
#print(color_code())


    