""" check.py: helper to maintain code option consistency """
import re 
import sys
import os

def parseIf(string, defines, fin):
    """ Parse pre-processor line from a source file. """
    s = string
    
    while len(s) > 0:
        if s.startswith("defined"):
            s = s[7:]
            continue
            
        if s.startswith("//"):
            return
            
        m = re.search("[a-zA-Z_][a-zA-Z_0-9]*",s)
        if m is not None and m.start() == 0:
            defines.update([m.group()])
            #print "%s : %s"%(string,m.group())
            s = s[m.span()[1]:]
            continue
        
        if s.startswith("/*"):
            m = re.search("\*/",s)
            if m is not None:
                s = s[m.span()[1]:]
                continue
            else:
                return
            
        if s[0] == "\n":
            return
            
        if s == "\\\n":
            string = fin.readline()
            s = string
            continue
            
        if s[0] in '<>+-/=!&|() 0123456789\t':
            s = s[1:]
            continue
            
        print("Strange character in '%s' detected: '%s', skipping it."%(string, s[0]))
        s = s[1:]

def load(fin):
    """ Wrapper for loading a list of strings from file. """
    s = set()
    
    for i in fin:
        if not i.startswith("#"):
            s.update([i.strip()])
        
    return s
    
def write(s, fout):
    """ Wrapper for writing out a sorted list of strings to file. """
    fout = open(fout,'w')
    for i in sorted(s):
        fout.write(i + os.linesep)
        
    fout.close()

# ----------------------------

def filter_code(fin):
    """ Find all macros used  in a .c or .h files. """
    defines = set()
    
    line = fin.readline()
    while line != "":
        s = line.lstrip()

        if s.startswith("#if "):
            parseIf(s[4:],defines,fin)
        elif s.startswith("#elseif "):
            parseIf(s[8:],defines,fin)
        elif s.startswith("#elif "):
            parseIf(s[6:],defines,fin)    
            
        elif s.startswith("#ifdef ") or s.startswith("#ifndef "):
            if s.startswith("#ifdef "):
                s = s[7:]
            else:
                s = s[8:]
                
            s = s.lstrip()
            m = re.search("[a-zA-Z_][a-zA-Z_0-9]*",s)
            if m is not None and m.start() == 0:
                defines.update([m.group()])
            else:
                print("Strange #ifdef/#ifndef: '%s'. Skipping it.",s)
            
        line = fin.readline()
    
    return defines

def filter_code_params():
    """ Find all parameter names that can be read from a parameterfile. """
    with open('src/parameters.c','r') as fin:
        s = fin.read()

    d = re.findall(r'strcpy\(tag\[nt\],\s?"[a-zA-Z_0-9]+"\);',s)

    params = [dd.split('"')[1] for dd in d]

    return set(params)
    
def filter_template_config(fin):
    """ Find all items of Template-Config.sh file. """
    defines = set()
    for line in fin:
        s = line.split()
        if(len(s)>0):
            d = re.findall("^#*([a-zA-Z_][a-zA-Z_0-9]*)",s[0])
            for dd in d:
                defines.update([dd])
    
    return defines

def filter_config(fin):
    """ Find all active items on Config.sh file. """
    defines = set()
    for line in fin:
        s = line.split()
        if(len(s)>0):
            d = re.findall("^([a-zA-Z_][a-zA-Z_0-9]*)",s[0])
            for dd in d:
                defines.update([dd])
    
    return defines
    
def filter_makefile(fin):
    """ Find all macros used in Makefile. """
    defines = set()
    for line in fin:
        s = line.strip()
        if s.startswith("ifeq"):
            d = re.findall("ifeq\s*\(([a-zA-Z_][a-zA-Z_0-9]*)\s*,\s*\$\(findstring",s)
            for dd in d:
                defines.update([dd])
        if s.startswith("ifneq"):
            d = re.findall("ifneq\s*\(([a-zA-Z_][a-zA-Z_0-9]*)\s*,\s*\$\(findstring",s)
            for dd in d:
                defines.update([dd])
    
    return defines

def filter_config_documentation():
    """ Parse for all documented Config.sh options (in documentation/*). These can appear either 
        as base options in core_config_options.md or as specialized options in any individual 
        modules_*.md file. """
    import glob
    defines = set()

    files = ['documentation/core_config_options.md']
    files += glob.glob('documentation/modules_*.md')

    for filename in files:
        with open(filename,'r') as fin:
            s = fin.read()

        d = re.findall(r"\n\n[a-zA-Z_0-9]+\n  ",s)

        defines.update( [dd.strip() for dd in d] )

    return defines
    
def filter_params_documentation():
    """ Parse for all documented Config.sh options (in documentation/*). These can appear either 
        as base options in core_param_options.md or as specialized options in any individual 
        modules_*.md file. """
    import glob
    defines = set()

    # syntax in core_param_options.md file: same as Config.sh options above
    with open('documentation/core_param_options.md','r') as fin:
        s = fin.read()

    d = re.findall(r"\n\n[a-zA-Z_0-9]+\n  ",s)
    defines.update( [dd.strip() for dd in d] )

    # syntax in individual module files: * ``ParamName`` description here.
    files = glob.glob('documentation/modules_*.md')

    for filename in files:
        with open(filename,'r') as fin:
            s = fin.read()

        d = re.findall(r"\n\* ``[a-zA-Z_0-9]+``",s)
        defines.update( [dd.strip()[4:-2] for dd in d] )

    return defines

# ----------------------------

def check_code(fin, fout, template, extra):
    """ Check source files for illegal macros. """
    allowed = filter_template_config(template)
    allowed.update(filter_config(extra))
    
    used = filter_code(fin)    
    diff = sorted(used.difference(allowed))
    
    if len(diff) > 0:
        print("\nIllegal macros/options detected in file %s.\nCheck for potential typos and add them either to the file 'Template-Config.sh' or 'defines_extra'\n('defines_extra' is for macros which are either internally defined or should not appear in Template-Config.sh).\nIn case you want to suppress this check, build with 'make build' instead of 'make'.\n"%fin.name)
        for i in diff:
            print(i)
        print("")
        exit(1)
        
    write(used,fout)
    exit(0)
    
def check_makefile(fin, fout, template, extra):
    """ Check Makefile for illegal options. """
    allowed = filter_template_config(template)

    allowed.update(filter_config(extra))
    
    used = filter_makefile(fin)
    diff = sorted(used.difference(allowed))
    
    if len(diff) > 0:
        print("\nIllegal macros/options detected in file %s.\nCheck for potential typos and add them either to the file 'Template-Config.sh' or 'defines_extra'\n('defines_extra' is for macros which are either internally defined or should not appear in Template-Config.sh).\nIn case you want to suppress this check, build with 'make build' instead of 'make'.\n"%fin.name)
        for i in diff:
            print(i)
        print("")
        exit(1)
        
    write(used,fout)
    exit(0)
    
def check_config(fin, fout, args, extra):
    """ Check Config.sh for illegal options (those which don't appear in any source files). """
    allowed = filter_config(extra)
    
    for file in args:
        allowed.update(load(open(file,'r')))

    used = filter_config(fin)
    diff = sorted(used.difference(allowed))
    
    if len(diff) > 0:
        print("\nThe following options are active in %s, but are not used in any of the\nsource code files being compiled into the final executable.\nPlease check for typos and deactivate the options.\nIn case you want to suppress this check, build with 'make build' instead of 'make'.\n"%fin.name)
        for i in diff:
            print(i)
        print("")
        exit(1)
        
    write(used,fout)
    exit(0)

def check_documentation(fin, fout):
    """ Check whether all Template-Config.sh and parameterfile options are documented. """
    ex = False

    # Template-Config.sh
    documented = filter_config_documentation()
    used = filter_template_config(fin)
    diff = sorted(used.difference(documented))
    
    if len(diff) > 0:
        print("\nWARNING: The following options are undocumented, but appear in %s" % (fin.name))
        print("Please add a proper documentation description:\n")
        for i in diff:
            print(i)
        print("")
        ex = True
        
    diff = sorted(documented.difference(used))
    
    if len(diff) > 0:
        print("\nWARNING: The following options are documented, but are (not used/deleted/misspelled) in %s:\n" % (fin.name))
        for i in diff:
            print(i)
        print("")
        ex = True

    # parameterfile
    documented = filter_params_documentation()
    allowed = filter_code_params()
    diff = sorted(allowed.difference(documented))

    if len(diff) > 0:
        print("\nWARNING: The following parameters are undocumented, but appear in src/parameter.c")
        print("Please add a proper documentation description:\n")
        for i in diff:
            print(i)
        print("")
        ex = True
        
    diff = sorted(documented.difference(allowed))
    
    if len(diff) > 0:
        print("\nWARNING: The following parameters are documented, but are (not used/deleted/misspelled) in src/parameter.c:\n")
        for i in diff:
            print(i)
        print("")
        ex = True
        
    if ex:
        exit(1)
        
    write(used,fout)
    exit(0)    
    
if __name__ == "__main__":
    if len(sys.argv) < 3:
        exit(1)
        
    mode = int(sys.argv[1])
    if mode < 1 or mode > 4:
        print("Unknown mode")
        exit(1)
        
    fin = open(sys.argv[2],'r')
    fout = sys.argv[3]
    
    if mode == 1:
        print("Checking %s for illegal define macros"%sys.argv[2])
        template = open(sys.argv[4],'r')
        extra = open(sys.argv[5],'r')
        
        check_code(fin, fout, template, extra)
        
    if mode == 2:
        print("Checking active options of %s"%sys.argv[2])
        extra = open(sys.argv[4],'r')
        check_config(fin, fout, sys.argv[5:], extra)
        
    if mode == 3:
        print("Checking %s for illegal define macros"%sys.argv[2])
        template = open(sys.argv[4],'r')
        extra = open(sys.argv[5],'r')
        
        check_makefile(fin, fout, template, extra)
        
    if mode == 4:
        print("Checking %s for undocumented options"%sys.argv[2])
        template = open(sys.argv[2],'r')
        
        check_documentation(fin, fout)
        
        