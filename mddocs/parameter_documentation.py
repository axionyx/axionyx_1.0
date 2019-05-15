import sys
import re
import os
import numpy as np

BASEPATH = os.path.dirname(__file__)
dirs = [BASEPATH+"/../Source/",BASEPATH+"/../Exec/Test_Only_Axions/"]
our_dict = {}
params = []
def scan_files(dirs,recursive=False):
    if(recursive):
        for root, dirs, files in os.walk("../Source/"):
                for file in files:
                    if(r".cpp" in file or r".H" in file):
                            scan_file(os.path.join(root,file))

    else:
        for D in dirs:
            files = os. listdir(D)
            for f in files:
                if(r".cpp" in f or r".H" in f):
                    scan_file(f)
def scan_file(fname):
    with open(fname) as f:
        lines = f.readlines()
    pname=None
    for i, line in enumerate(lines):
#        print line 
        r=re.search(r"ParmParse (.*?)(\(\"(.*?)\"\)|);",line)
        if(r is not None):
            pname = r.group(1)
            psubdir = r.group(3)
        
        
        if(pname is not None):
            r=re.search(r"({}\.get|{}\.query|{}\.getarr|{}\.queryarr)\(\"(.*?)\"".format(pname,pname,pname,pname),line)
            if (r is not None):
                param = r.group(2)
                if(psubdir!=None):
                    param = psubdir + "." +  param
                if("get" in r.group(1)):
                        type_string = ("_(obligatory)_")
                else:
                        type_string = ("_(optional)_")
                if(not param in params):
                    params.append(param)
                #print "found parameter",param
                if(i==0):
                    continue
                myline = lines[i-1]
                r=re.search(r"(\\pparam) (.*)", myline) 
                if (r is not None):
                    docstring = r.group(2)
                    #print "found docstring:",docstring
                    entry = {"location":[fname,i],"docstring":type_string+" "+docstring}
                    if(param in our_dict):
                        if (our_dict[param]["docstring"] != ""):
                            print "WARNING: parameter",param,"has two doc strings. Originally in",our_dict[param]["location"],
                            "again found in", entry["location"]
                    else:
                        our_dict[param] = entry
def check_dict():
    for p in params:
        if(p not in our_dict.keys()):
            print "WARNING: parameter",p,"is not documented!"
def get_doc_line(key,entry):
    s = ""
    s +=r"`"+key+"`- "+entry["docstring"]+"  \n"
    return s
def get_header():
    s = "## Input Parameters\nThis section briefly describes the parameter available in the inputs file.\n\n"
    return s
def write_params_doc(fname):
    ks = our_dict.keys()
    ks.sort()
    with open(fname,"w") as f:
        f.write(get_header())
        for key in ks:
            f.write(get_doc_line(key,our_dict[key]))
scan_files(dirs,recursive=True)    
check_dict()
write_params_doc(BASEPATH+"/input_parameters.md")
