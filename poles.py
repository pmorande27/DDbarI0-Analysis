import os
import particle
import matplotlib.pyplot as plt
import numpy as np

def pole_extraction_from_xml(path_xml):
    """
    Extracts the poles from the xml file
    :param path_xml: path to the xml file
    :return: list of poles
    """
    poles = []
    pole = False
    pole_find = False
    in_values = False
    sheets = []
    sqrt_ss_r = []
    sqrt_ss_im = []
    sqrt_ss_r_err = []
    sqrt_ss_im_err = []
    sqrt_s = False
    with open(path_xml, 'r') as f:
        for line in f:
            if '<Poles>' in line:
                pole = True

            if pole:
                if "<elem>" in line and pole_find == False:
                    pole_find = True
                    continue
                if pole_find and "<Val>" in line:
                    in_values = True
                    continue
                if pole_find and "</elem>" in line and in_values == False:
                    pole_find = False
                    continue
                if "</Val>" in line:
                    in_values = False
                    continue
                if "<First>" in line:
                    twoJ = line.split("<First>")[1].split("</First>")[0]
                if "<Second>" in line:
                    P = line.split("<Second>")[1].split("</Second>")[0]
                if in_values:
                    
                        
                        

                    
                    if "sheet" in line:
                        sheets.append((line,twoJ,P))
                    if "<sqrt_s_pole>" in line:
                        sqrt_s = True
                        continue
                    if "</sqrt_s_pole>" in line:
                        sqrt_s = False
                    if sqrt_s:
                        if "<elem>" in line:
                            continue
                        if "</elem>" in line:
                            sqrt_s = False
                            continue
                        if "<re_mean>" in line:
                            value = line.split("<re_mean>")[1].split("</re_mean>")[0]
                            sqrt_ss_r.append(value)
                        if "<im_mean>" in line:
                            value = line.split("<im_mean>")[1].split("</im_mean>")[0]
                            sqrt_ss_im.append(value)
                        if "<re_err>" in line:
                            value = line.split("<re_err>")[1].split("</re_err>")[0]
                            sqrt_ss_r_err.append(value)
                        if "<im_err>" in line:
                            value = line.split("<im_err>")[1].split("</im_err>")[0]
                            sqrt_ss_im_err.append(value)
    poles = {}
    for j in range(len(sheets)):
        l =sheets[j][0].split("<sheet>")[1].split("</sheet>")[0]
        sheet = l
        twoJ = sheets[j][1]
        P = sheets[j][2]
        re = sqrt_ss_r[j]
        im = sqrt_ss_im[j]
        re_err = sqrt_ss_r_err[j]
        im_err = sqrt_ss_im_err[j]
        if (sheet,twoJ,P) not in poles.keys():
            poles[(sheet,twoJ,P)] = [(re,im,re_err,im_err)]
        else:
            poles[(sheet,twoJ,P)].append((re,im,re_err,im_err))
    example_sheet =sheets[0][0].split("<sheet>")[1].split("</sheet>")[0]
    channels = example_sheet.split(" ")
    all_Particles =  particle.read_particles('Particles/particles_unfl.txt') + particle.read_particles('Particles/charmonium.txt')+ particle.read_particles('Particles/Ds.txt')
    ms = []
    mchan = {}
    for chan in channels:
        particles = chan.split(":")
        last = len(particles)-1
        new_str  = ""
        for let in particles[last]:
            if let != "[" and let != "]" and let != "+" and let != "-":
                new_str += let
        particles[last] = new_str
        for p in all_Particles:
            if p.name == particles[0]:
                m_1 = p.Mass
            if p.name == particles[1]:
                m_2 = p.Mass
        m_th = m_1 + m_2
        chan_w = ""
        for st in chan:
            if st != "+" and st != "-" and st != "[" and st != "]":
                chan_w += st
        mchan[chan_w] = m_th
        ms.append(m_th)
    new_dict = {}
    for key in poles.keys():
        sheet,twoJ,P = key
        new_key = (sheet,twoJ)
        chans = sheet.split(" ")
        masses = []
        for chan in chans:
            new_ch =""
            for st in chan:
                if st != "+" and st != "-" and st != "[" and st != "]":
                    new_ch += st
            masses.append(mchan[new_ch])
        "order chans by mass in ascending order"

        for i in range(len(masses)):
            for j in range(i+1,len(masses)):
                if masses[j] < masses[i]:
                    masses[i],masses[j] = masses[j],masses[i]
                    chans[i],chans[j] = chans[j],chans[i]

        new_sheet =" ".join(chans)
        new_key = (new_sheet,twoJ,P)

        
        new_dict[new_key] = poles[key]
     
    return new_dict,ms

def plot_poles_from_xml(path_xml):
    """
    Plots the poles from the xml file
    :param path_xml: path to the xml file
    :return: None
    """
    colors = ['red','blue','green','orange','purple','yellow','pink','brown','cyan','magenta','grey','lime','olive','teal','navy','maroon','aqua','fuchsia','silver']

    sheets = []
    twoJsP = []
    

    poles,ms = pole_extraction_from_xml(path_xml)

    print(poles)
    for key in poles.keys():
        sheet,twoJ,P = key
        sheet_reformed = ""
        for i in range(len(sheet)):
            if sheet[i] != "-" and sheet[i] != "+":
                continue
            else:
                sheet_reformed += sheet[i]
        if sheet_reformed not in sheets:

            sheets.append(sheet_reformed)
        if (twoJ,P) not in twoJsP:
            twoJsP.append((twoJ,P))
    fix,axes = plt.subplots(len(twoJsP),1,figsize=(10,10))
    if len(twoJsP) == 1:
        axes = [axes]
    ax_dict = {}
    sheet_color_dict = {}
    #print(sheets)
    for i,sheet in enumerate(sheets):
        sheet_color_dict[sheet] = colors[i]
    labels_dict ={}
    for twoJ,P in twoJsP:
        ax_dict[(twoJ,P)] = axes[twoJsP.index((twoJ,P))]
        
        labels_dict[(twoJ,P)] = []

            
    ys  = [0 for i in range(len(ms))]
    for ax in axes:
        ax.set_xlabel("$\Re\sqrt{s}$")
        ax.set_ylabel("$\Im 2\sqrt{s}$")

        ax.plot(ms,ys,'o')
        ax.set_xlim(min(ms)-0.2,max(ms)+0.2)
        ns = np.linspace(min(ms)-0.2,max(ms)+0.2,100)
        ax.plot(ns,ns*0,'k--')
    for key in poles.keys():
        for j in range(len(poles[key])):
            sheet,twoJ,P = key
            sheet_reformed = ""
            for i in range(len(sheet)):
                if sheet[i] != "-" and sheet[i] != "+":
                    continue
                else:
                    sheet_reformed += sheet[i]
            color = sheet_color_dict[sheet_reformed]
            label = sheet_reformed
            ax = ax_dict[(twoJ,P)]
            re,im,re_err,im_err = poles[key][j]
            im = 2*float(im)
            im_err = 2*float(im_err)
            labels = labels_dict[(twoJ,P)]
            if label not in labels:
                labels.append(label)
                ax.errorbar(float(re),float(im),xerr=float(re_err),yerr=float(im_err),fmt='o',color=color,label=label)
                ax.text(float(re),float(im),f"({label})",color=color)
            else:
                ax.errorbar(float(re),float(im),xerr=float(re_err),yerr=float(im_err),fmt='o',color=color)
                ax.text(float(re),float(im),f"({label})",color=color)
    ns = np.linspace(min(ms)-0.2,max(ms)+0.2,100)
    for ax in axes:
        ax.legend()
    plt.show()

def plot_poles_from_xml_list(path_xml_list):
    """
    Plots the poles from the xml file
    :param path_xml: path to the xml file
    :return: None
    """
    cond = False

    colors = ['red','blue','green','orange','purple','yellow','pink','brown','cyan','magenta','grey','lime','olive','teal','navy','maroon','aqua','fuchsia','silver']
    for k,path_xml in enumerate(path_xml_list):
        sheets = []
        twoJsP = []
        

        poles,ms = pole_extraction_from_xml(path_xml)
        #print(poles)
        for key in poles.keys():
            sheet,twoJ,P = key
            sheet_reformed = ""
            for i in range(len(sheet)):
                if sheet[i] != "-" and sheet[i] != "+":
                    continue
                else:
                    sheet_reformed += sheet[i]
            if sheet_reformed not in sheets:

                sheets.append(sheet_reformed)
            if (twoJ,P) not in twoJsP:
                twoJsP.append((twoJ,P))
        if cond == False:
            print("cond")
            fix,axes = plt.subplots(len(twoJsP),1,figsize=(10,10))
            
            ax_dict = {}
            sheet_color_dict = {}
            #print(sheets)
            for i,sheet in enumerate(sheets):
                sheet_color_dict[sheet] = colors[i]
            labels_dict ={}
            for twoJ,P in twoJsP:
                ax_dict[(twoJ,P)] = axes[twoJsP.index((twoJ,P))]
                
                labels_dict[(twoJ,P)] = []
            cond = True
                
        ys  = [0 for i in range(len(ms))]
        for ax in axes:
            ax.set_xlabel("$\Re\sqrt{s}$")
            ax.set_ylabel("$\Im 2\sqrt{s}$")

            ax.plot(ms,ys,'o')
            ax.set_xlim(min(ms)-0.2,max(ms)+0.2)
            ns = np.linspace(min(ms)-0.2,max(ms)+0.2,100)
            ax.plot(ns,ns*0,'k--')
        for key in poles.keys():
            for j in range(len(poles[key])):
                sheet,twoJ,P = key
                sheet_reformed = ""
                for i in range(len(sheet)):
                    if sheet[i] != "-" and sheet[i] != "+":
                        continue
                    else:
                        sheet_reformed += sheet[i]
                color = sheet_color_dict[sheet_reformed]
                label = sheet_reformed
                ax = ax_dict[(twoJ,P)]
                re,im,re_err,im_err = poles[key][j]
                im = 2*float(im)
                im_err = 2*float(im_err)
                labels = labels_dict[(twoJ,P)]
                alpha = (k+1)/len(path_xml_list)
                if im <=1:
                    if label not in labels:
                        labels.append(label)
                        ax.errorbar(float(re),float(im),xerr=float(re_err),yerr=float(im_err),fmt='o',color=color,label=label,alpha=alpha)
                        ax.text(float(re),float(im),f"({label})",color=color,alpha=alpha)
                    else:
                        ax.errorbar(float(re),float(im),xerr=float(re_err),yerr=float(im_err),fmt='o',color=color,alpha=alpha)
                        ax.text(float(re),float(im),f"({label})",color=color,alpha=alpha)
        ns = np.linspace(min(ms)-0.2,max(ms)+0.2,100)
        for ax in axes:
            ax.legend()
    plt.show()

    

                    
                        



                

      
#plot_poles_from_xml_list(['Data/poles.xml',])
plot_poles_from_xml('Data/poles2.xml',)