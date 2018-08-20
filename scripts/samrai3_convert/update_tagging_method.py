import sys
import os
from shutil import copyfile

def find_database(file_name):
    file_handle = open(file_name, "r")
    fl = file_handle.readlines()
    strs = []
    found_db = -1
    num_bracks = 0
    for x in fl:
        if "StandardTagAndInitialize" in x:
            found_db = 1
            if "{" in x:
                num_bracks = num_bracks+1
            strs.append(x)
            continue
        # Read until we find matching }
        if ((found_db == 1) & ("{" in x)):
            num_bracks = num_bracks+1
        if ((found_db == 1) & ("}" in x)):
            num_bracks = num_bracks-1
        if found_db == 1:
            # Get rid of comments
            if x.find("//") > -1:
                x = x[:x.find("//")]+"\n"
                if x == "\n":
                    continue
            strs.append(x)
        if ((found_db == 1) & (num_bracks == 0)):
            found_db = 0
            break
    file_handle.close()
    return strs

def find_in_list(strs, to_find):
    # Search in list of strings for particular string, then return that line.
    # Note only finds first occurence.
    for x in strs:
        if to_find in x:
            return x
    return ""

def parse_db(strs):
    # Create dictionary of relevant info.
    parsed = dict()
    # First check for tagging method
    x = find_in_list(strs, "tagging_method")
    parsed['tagging_method'] = []
    if ("REFINE_BOXES" in x):
        parsed['tagging_method'].append("REFINE_BOXES")
    if ("GRADIENT_DETECTOR" in x):
        parsed['tagging_method'].append("GRADIENT_DETECTOR")
    # Now we check for level info. Let's max out on 10 levels...
    l = 0
    while True:
        level = "level_"+str(l)
        x = find_in_list(strs, level)
        if x == "":
            break
        rhs = x.split(" = ")
        parsed[level] = rhs[1]
        l = l+1
    return parsed

def print_database(file_name, strs_to_print):
    # Similar to find_database, but this one replaces the
    # StandardTagAndInitialize{} database
    fh = open(file_name, 'r')
    fh_temp = open(file_name+".tmp", 'w')
    write_old = 1
    done_writing_new = 0
    num_brackets = 0
    for line in fh:
        if ("StandardTagAndInitialize" in line):
            write_old = 0
        if ((write_old == 0) & ("{" in line)):
            num_brackets = num_brackets+1
        if ((write_old == 0) & ("}" in line)):
            num_brackets = num_brackets-1
        if(write_old == 1):
            fh_temp.write(line)
        elif(done_writing_new == 0):
            for new_lines in strs_to_print:
                fh_temp.write(new_lines)
            done_writing_new = 1
        if((done_writing_new == 1) & (num_brackets == 0)):
            write_old = 1
    fh.close()
    fh_temp.close()
    copyfile(file_name+".tmp", file_name)
    os.remove(file_name+".tmp")


def generate_new_db(db):
    types = db["tagging_method"]
    i = 0
    strs = ["StandardTagAndInitialize {\n", 
            "    at_0 {\n",
            "        time = 0.0\n"]
    for type in types:
        if (type == "REFINE_BOXES"):
            strs.extend(generateRefineBoxes(db, i))
            i = i+1
        elif(type == "GRADIENT_DETECTOR"):
            strs.extend(generateGradientDetector(db, i))
            i = i+1
    strs.append("    }\n")
    strs.append("}\n")
    return strs

def generateRefineBoxes(db, i):
    strs = ["        tag_"+str(i)+" {\n",
            "            tagging_method = \"REFINE_BOXES\"\n"]
    l = 0

    while(("level_"+str(l)) in db):
        strs.append("            level_"+str(l)+" {\n")
        strs.append("                boxes = "+db["level_"+str(l)])
        strs.append("            }\n")
        l = l+1
    strs.append("        }\n")
    return strs

def generateGradientDetector(db, i):
    strs = ["        tag_"+str(i)+" {\n",
            "            tagging_method = \"GRADIENT_DETECTOR\"\n",
            "        }\n"]
    return strs


file_dir = sys.argv[1]
file_list = []
for subdir, dirs, files in os.walk(file_dir):
    for file in files:
        if ("input" in file):
            file_list.append(os.path.join(subdir, file))

for file in file_list:
    print "Working on file: " + file
    str_list = find_database(file)
    db = parse_db(str_list)
    str_list = generate_new_db(db)
    print_database(file, str_list)
