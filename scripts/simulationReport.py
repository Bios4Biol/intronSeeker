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

def simulationReport(r1,r2) :
    print(r1.name)
    print(r2.name)


if __name__ == '__main__' :
    
    import argparse 
    
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('-1','--R1', type=argparse.FileType('r'), required=True, dest='r1')
    parser.add_argument('-2','--R2', type=argparse.FileType('r'), required=False, dest='r2')
    
    args = vars(parser.parse_args())
    
    simulationReport(**args)
