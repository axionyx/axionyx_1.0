import sys
import re
import os
import numpy as np

BASEPATH = os.path.dirname(os.path.abspath(__file__))
our_dict = {}
params = []
def scan_files():
    for root, dirs, files in os.walk(BASEPATH+"/../Exec/"):
        for file in files:
            if(r".cpp" in file or r".H" in file or r".f90" in file):
                    full_path = os.path.join(root,file)
                    print full_path
                    scan_file(full_path)

def scan_file(fname):
    with open(fname) as f:
        lines = f.readlines()
    pname=None
    for i, line in enumerate(lines):
#        print line 
        r=re.search(r"(\\testproblem) (.*)", line) 
        if(r is not None):
            docstring = r.group(2)
            entry = {"location":fname,"docstring":docstring}
            t = os.path.split(fname)[0]
            key = os.path.split(t)[1]
            if(key in our_dict):
                print "WARNING: "+key+" is documented more than once. First occurence in "+our_dict[key]["location"]+" found again in "+entry["location"]
            else:
                our_dict[key] = entry

def get_doc_line(key,entry):
    s = ""
    s +=r"`"+key+"`- "+entry["docstring"]+"  \n"
    return s
def get_header():
    s = "## Available Test Problems\nThis section briefly describes the FDM-related test problems available in `Exec/`\n\n"
    return s
def write_testproblem_doc(fname):
    ks = our_dict.keys()
    ks.sort()
    with open(fname,"w") as f:
        f.write(get_header())
        for key in ks:
            f.write(get_doc_line(key,our_dict[key]))
        
scan_files()    
write_testproblem_doc("test_problems.md")
