#!/bin/python3

#<IntronSeeker searches introns by splice-realigning reads on contigs.>
#Copyright (C) <2019-2024> INRAE
#<Sarah Maman, Philippe Bardou, Emilien Lasguignes, Faustine Oudin, Floréal Cabanettes, Christophe Klopp>
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <https://www.gnu.org/licenses/>.


try :
    import os ;
    import configparser ;
    import subprocess ;
    import re ;
    import sys ;
    from helpMessages import print_to_stdout ;
except ImportError as error :
    print(error) ;
    exit(1) ;

def checkInstall() :

    print_to_stdout('###  Start to check installation   ###')        
    
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
    print_to_stdout('###  Check installation finished   ###')        

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
    if asked_version == real_version:
        return True;
    for i in range(min(len(ref),len(real))) :
        if int(real[i]) < int(ref[i]) :
            return False ;
        if int(real[i]) > int(ref[i]) :     
            return True;


def print_to_stdout(*a):
 
    # Here a is the array holding the objects
    # passed as the argument of the function
    print(*a, file = sys.stdout)
