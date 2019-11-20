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

def simulationReport() :
    print()

