#!/usr/bin/env python3

import re
import pickle
import os
import pandas as pd
import pysam
import gzip
import time
import sys
import numpy as np
from pprint import pprint
from Bio import SeqIO
from collections import OrderedDict
import configparser
import concurrent.futures as prl
from itertools import repeat
from intronSeekerPlot import * 

def simulationReport(modifiedfa, gtf) :
    pass

if __name__ == '__main__' :
    
    import argparse 
    
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('-r','--reference', type=argparse.FileType('r'), required=True, dest='modifiedfa')
    parser.add_argument('-g','--gtf', type=argparse.FileType('r'), required=True, dest='gtf') 

    args = vars(parser.parse_args())
    
    simulationReport(**args)
