#!/bin/python3

try :
    import os ;
    import configparser ;
    import subprocess ;
    import re ;
    import sys ;
except ImportError as error :
    print(error) ;
    exit(1) ;

def checkInstall() :

    config_path = os.path.abspath(os.path.dirname(__file__))+"/../config/intronSeeker.properties"

    config = configparser.RawConfigParser() ;
    config.read(config_path) ;

    miss_dpdc = [] ;
    for program in config["Commands"] :
        if not checkProgram(program,
                        config["Commands"][program],
                        config["Versions"][program],
                        config["Warnings"][program]) :
            miss_dpdc.append(program) ;

    if not miss_dpdc :
        print("\nAll the dependencies are correctly installed.") ;
        print() ;
    elif len(miss_dpdc) == 1 :
        print("Finally, the following program is missing or is not up to date:") ;
        print("\t"+", ".join(miss_dpdc)) ;
        print() ;
    else : 
        print("Finally, the following programs are missing or are not up to date:") ;
        print("\t"+", ".join(miss_dpdc)) ;
        print() ;

def checkProgram(program : str, command : str, asked_version : str , warning : str) :
    try :
        print("\n"+program + " testing...",end="") ;
        if program == "gffread" :
            sdo = subprocess.run([command, "--version"], stderr=subprocess.PIPE).stderr ;
        else :
            sdo = subprocess.run([command, "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout ;
        current_ver = re.search(r'([\d.]{3,10})',str(sdo)).group(1).rstrip() ;
        if not checkVersion(asked_version,current_ver) :
            raise OldVersion() ;
        else : 
            print("OK ! ") ;
    except OSError as e :
        if e.errno == os.errno.ENOENT : 
            print("\n**WARNING** : {program} is missing.\n{warning} can't be performed.".format(program = program, warning = warning)) ;
            return False ;
    except AttributeError :
        print("\n**WARNING** : {program} is missing.\n{warning} can't be performed.".format(program = program, warning = warning)) ;
        return False ;
    except OldVersion as ov :
        print("\n**WARNING** : {program} version is too old. The {current_ver} is installed, \
            the {asked_version} is necessary.\n{warning} may not be performed.".format(
                program = program,
                current_ver = current_ver,
                asked_version = asked_version,
                warning = warning
                )) ;
        return False ;
    return True ;

class OldVersion(Exception) :
    pass

def checkVersion(asked_version : str, real_version : str) :
    ref = asked_version.split(".") ;
    real = real_version.split(".") ;
    for i in range(min(len(ref),len(real))) :
        if int(real[i]) < int(ref[i]) :
            return False ;
    return True ;


if __name__ == "__main__" :
    checkInstall() ;
