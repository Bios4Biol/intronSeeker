#!/usr/bin/env python3

import os
import configparser
import numpy as np
import pandas as pd
import plotly as py
import plotly.graph_objects as go
import plotly.figure_factory as ff
import plotly.subplots as psp
import argparse
import pysam   # To generate a dataframe from a BAM : pysam and pickle
import pickle
from plotly.subplots import make_subplots   #for flagstat pie
import re   #for flagstat pie



import gzip
import time
import sys
import concurrent.futures as prl
import tempfile
from pprint import pprint
from collections import OrderedDict
from itertools import repeat
from Bio import SeqIO

from json import JSONEncoder
import json

# Distribution plot
def plot_dist(len_by_features, feature_names):
    hist_data = len_by_features
    group_labels = feature_names
    colors = ['#333F44', '#37AA9C', '#94F3E4']
    # Create distplot with curve_type set to 'normal'
    fig = ff.create_distplot(hist_data, group_labels, show_hist=False, show_rug=True, colors=colors)
    fig.update_layout(
        margin=go.layout.Margin(
            l=50,
            r=50,
            b=20,
            t=30,
            pad=0
        )
    )
    return py.offline.plot(fig, include_plotlyjs=False, output_type='div')


# Plot introns position on a contigs
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
    layout = go.Layout(xaxis=dict(
                           title="% of contig length"),
                       yaxis=dict(
                           title="Count"))
    fig = go.Figure(data=[hist],layout=layout)
    fig.update_layout(
        margin=go.layout.Margin(
            l=50,
            r=50,
            b=20,
            t=30,
            pad=0
        )
    )
    return py.offline.plot(fig, include_plotlyjs=False, output_type='div')

# Plot ranks file
def abondance_model(rank:dict, real_abund_perc:dict) :
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x = rank,
            y = real_abund_perc,
            mode = 'lines',
            name = 'Simulated abundance model'
        )
    )

    fig.add_trace(
        go.Scatter(
            x = rank,
            y = real_abund_perc,
            mode = 'lines',
            name ='Waited abundance model'
        )
    )

    fig.update_layout(
        xaxis=dict(title="Contigs"),
        yaxis=dict(title="Relative Abundance percentage",
            range=[-0.25,0.5])
    )
    return py.offline.plot(fig, include_plotlyjs=False, output_type='div')

# Barplots of split reads signal detection for HiSAT2
def plot_class_reads(*args,**kwargs) :
    series = [[],[]]
    colors = kwargs['colors']
    
    fig = psp.make_subplots(
        rows = 2, cols=1,
        shared_xaxes=True,
        vertical_spacing=0.03,
        specs=[[{"type":"table"}],
               [{"type":"bar"}]])

    for name,val in args :
        for idx,val in enumerate([val,val.loc[lambda df : df.contig == df.contig.str.rstrip(".ori"),:]]) :
            s = val.loc[:,'classe'].value_counts()
            s.name=name
            
        
            if idx == 0 :
                visibility = True
            else :
                visibility = False
        
            fig.add_trace(
                go.Bar(
                    visible=visibility,
                    x=s.index[1:],
                    y=s.values[1:]/s.sum()*100,
                    name = s.name,
                    marker_color=colors[s.name]
                ),row=2,col=1)
            
            s['Total'] = s.sum()
            series[idx].append(s)

    for idx,serie in enumerate(series) :
        table = pd.concat(serie,axis=1,sort=False).fillna(0)
        
        
        if idx == 0 :
            visibility = True
        else :
            visibility = False
            
        fig.add_trace(
            go.Table(
                visible=visibility,
                columnwidth=[40,80],
                header = dict(
                    values=["",*[c.join(["<b>","</b>"]) for c in table.columns.to_list()]],
                    fill_color='lightgrey'
                    ),
                cells = dict(
                    values=[[v.join(["<b>","</b>"]) for v in list(table.reset_index().T.values)[0]],
                           *list(table.reset_index().T.values)[1:]],
                    fill=dict(color=['lightgrey','lavender'])
                    )
                ),
            row=1,col=1)
    
    fig.update_layout(
        height=700,
        legend=dict(
            x=-0.4,y=0.40
            ),
        updatemenus=[
            go.layout.Updatemenu(
                buttons=[
                    dict(
                        args=[{"visible":[True,False]}],
                        label='All alignements',
                        method='restyle'
                        ),
                    dict(
                        args=[{"visible":[False,True]}],
                        label='Alignements on contigs with introns',
                        method='restyle'
                        ),
                    ]
                )
            ]
        )
    return py.offline.plot(fig, include_plotlyjs=False, output_type='div')

# Plot : Counting table and barplots of mapped covering reads' main characteristics.
#def plot_covering_reads(*args, **kwargs) :
def plot_covering_reads(*args, library:dict, colors:dict) :
    series=[]
    fig = go.Figure()
    for name, val in args :
        s = pd.Series(name=name)
        s['Covering']=len(val)
        s['Unmapped']=len(val.loc[lambda df : df.mapped == False])
        tmp = val.merge(library,left_on='read',right_index=True,suffixes=("","_lib"))
        s['Mismapped'] = len(tmp.loc[lambda df : (df.contig.str.rstrip('.ori') != df.contig_lib.str.rstrip('.ori'))& (df.mapped == True)])
        s['Unsplit'] = len(val.loc[lambda df : (df.mapped==True)&(df.split==False)])
        s['Missplit'] = len(val.loc[lambda df : (df.split==True)&(df.missplit==True)])
        s['Correct splitting'] = len(val.loc[lambda df : df.classe == 'TP'])
        series.append(s)
        
        to_plot = s[['Unmapped','Unsplit','Correct splitting']]/s['Covering']*100
        fig.add_trace(
            go.Bar(
                    x=to_plot.index,
                    y=to_plot.values,
                    name = s.name,
                    marker_color=colors[s.name]
            ))
    table = pd.concat(series,axis=1,sort=False)
        
    fig.update_layout(
        title='Global mapping results on introns-covering reads',
        xaxis=dict(title='Lectures charecteristics'),
        yaxis=dict(title='Percentage of total covering reads alignements')
        )

    return py.offline.plot(fig, include_plotlyjs=False, output_type='div'), table

# Return int from flagstat HISAT mapping in string format (only for the last 3 values : Mapped, Properly paired, Singletons)
def cast_to_value_hisat(str_mapping:str, tot_hisat:int):
    #p = int(re.sub(r'(\([a-zA-Z0-9_]*.[a-zA-Z0-9_]*%\))', r" ", df_flag_all.iloc[2,1]))
    mapping=re.sub(r'(\([a-zA-Z0-9_]*.[a-zA-Z0-9_]*%\))', r" ", str_mapping)
    val=(int(mapping)*100)/tot_hisat
    return val

# Return int from flagstat HISAT mapping in string format
def cast_to_value_star(str_mapping:str, tot_star:int):
    mapping=re.sub(r'(\([a-zA-Z0-9_]*.[a-zA-Z0-9_]*%\))', r" ", str_mapping)
    val=(int(mapping)*100)/tot_star
    return val

# Pie Chart with mapping stats from flagstat files
# https://plot.ly/python/pie-charts/
def plot_flagstat(df_flag_all:dict):
    labels = ["Secondary","Mapped", "Properly paired", "Singletons"]

    # Create subplots: use 'domain' type for Pie subplot
    fig = make_subplots(rows=1, cols=2, specs=[[{'type':'domain'}, {'type':'domain'}]])
    tot_star=int(df_flag_all.iloc[0,1])
    fig.add_trace(go.Pie(labels=labels, values=[round((int(df_flag_all.iloc[1,1])*100)/tot_star, 2), round(cast_to_value_star(df_flag_all.iloc[2,1],tot_star),2), round(cast_to_value_star(df_flag_all.iloc[3,1],tot_star),2), round(cast_to_value_star(df_flag_all.iloc[4,1],tot_star),2)], name="STAR"),
                1, 1)
    tot_hisat=int(df_flag_all.iloc[0,0])
    fig.add_trace(go.Pie(labels=labels, values=[round((int(df_flag_all.iloc[1,0])*100)/tot_hisat, 2), round(cast_to_value_hisat(df_flag_all.iloc[2,0],tot_hisat),2), round(cast_to_value_hisat(df_flag_all.iloc[3,0],tot_hisat),2), round(cast_to_value_hisat(df_flag_all.iloc[4,0],tot_hisat),2)], name="HISAT2"),
                1, 2)

    # Use `hole` to create a donut-like pie chart
    fig.update_traces(hole=.4, hoverinfo="label+percent+name")

    fig.update_layout(
        # Add annotations in the center of the donut pies.
        annotations=[dict(text='STAR', x=0.20, y=0.5, font_size=12, showarrow=False),
                     dict(text='HISAT', x=0.80, y=0.5, font_size=12, showarrow=False)])

    return py.offline.plot(fig, include_plotlyjs=False, output_type='div')