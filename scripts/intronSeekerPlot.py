#!/usr/bin/env python3


#<IntronSeeker searches introns by splice-realigning reads on contigs.>
#Copyright (C) <2019> INRAE
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


import re
import os
import pandas as pd
import pysam
import gzip
from pprint import pprint
from Bio import SeqIO
from collections import OrderedDict
import plotly as py
from plotly import tools
import plotly.graph_objs as go
import plotly.subplots as psp
import plotly.io as pio
from ipywidgets import widgets
import numpy as np


def plot_insertion_in_contig(positions) :
    hist = go.Histogram(
            x=positions,
            xbins=dict(
                start=0,
                end=100,
                size=2),
            marker=dict(
                color='purple'
            )
    )
    layout = go.Layout(title='Distribution of introns insertion position along the contigs',
                       xaxis=dict(
                           title="% of contig length"),
                       yaxis=dict(
                           title="Count"))
    fig = go.Figure(data=[hist],layout=layout)
    print(py.offline.plot(fig, include_plotlyjs=True,output_type='div'))    
    # ~ pio.write_html(fig,'fig.html')
    # ~ return tools.get_embed(fig)
    
    fig.write_html('is_figure.html', auto_open=True)
    fig = go.Figure(data=go.Bar(y=[2, 3, 1]))
    fig.write_html('first_figure.html', auto_open=True)
