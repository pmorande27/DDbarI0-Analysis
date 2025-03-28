#!/usr/bin/python
import os
import sys
def init_xml_create(path_init_xml_inital,path_init_xml_final,tmax,t0_max,t0_min):
    with open(path_init_xml_inital) as f:
        lines = f.readlines()
        with open(path_init_xml_final,"w") as f:
            for line in lines:
                if "tmax" in line:
                    f.write("  <tmax>"+str(tmax)+"</tmax>\n")
                if "t0low" in line:
                    f.write("  <t0low>"+str(t0_min)+"</t0low>\n")
                if "t0high" in line:
                    f.write("  <t0high>"+str(t0_max)+"</t0high>\n")
                else:
                    f.write(line)
def create_tmax_spectrum_files(irrep,tmax_max,tmax_min,t0_max,t0_min):
    tmaxs = range(tmax_min,tmax_max+1)
    for t in tmaxs:
        os.mkdir(irrep + "_tmax_" + str(t))
        path_init_xml_inital = "./sfit.ini.xml"
        path_init_xml_final = irrep + "_tmax_" + str(t) + "/sfit.ini.xml"
        init_xml_create(path_init_xml_inital,path_init_xml_final,t,t0_max,t0_min)
def main(irrep,tmax_max,tmax_min,t0_max,t0_min):

    create_tmax_spectrum_files(irrep,tmax_max,tmax_min,t0_max,t0_min)
if __name__ == "__main__":
    tmax_max = int(sys.argv[1])
    tmax_min = int(sys.argv[2])
    irrep = sys.argv[3]
    t0_max = int(sys.argv[4])
    t0_min = int(sys.argv[5])

    main(irrep,tmax_max,tmax_min)

