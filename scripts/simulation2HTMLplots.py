#!/usr/bin/env python3

import pandas as pd
import plotly as py
import numpy as np 
import plotly.graph_objects as go
import plotly.figure_factory as ff
import plotly.subplots as psp
import re
import plotly.express as px 


# Histogram   https://plotly.com/python/v3/histograms/
def plot_hist_contigs_len(fastaContigsLen, mFastaContigsLen):
    """
    plot_hist_contigs_len function plot an histogram of contigs length from fasta and modified fasta contigs length data.
    """
    contigs = go.Histogram(
        x=fastaContigsLen,
        name='Contigs',
        opacity=0.75
    )
    modifiedContigs = go.Histogram(
        x=mFastaContigsLen,
        name='Modified contigs',
        opacity=0.75
    )
    data = [contigs, modifiedContigs]
    layout = go.Layout(
        xaxis=dict(
            title='Contigs length'
        ),
        yaxis=dict(
            title='Number of contigs'
        ),
        bargap=0.2,
        bargroupgap=0.1
    )
    fig = go.Figure(data=data, layout=layout)
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

# SARAH 01/2021    
def plot_hist_contigs_len_real_data(fastaContigsLen):
    """
    plot_hist_contigs_len function plot an histogram of contigs length from fasta and modified fasta contigs length data.
    """
    contigs = go.Histogram(
        x=fastaContigsLen,
        name='Contigs',
        opacity=0.75
    )
    data = [contigs]
    layout = go.Layout(
        xaxis=dict(
            title='Contigs length'
        ),
        yaxis=dict(
            title='Number of contigs'
        ),
        bargap=0.2,
        bargroupgap=0.1
    )
    fig = go.Figure(data=data, layout=layout)
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


def plot_hist(values:dict, title:str, xtitle:str, ytitle:str):
    """
    plot_hist function is a reusable function to plot an histogram.
    """
    myhist = go.Histogram(
        x=values,
        name=title,
        xbins=dict( # bins used for histogram
        start=1,
        size=1
        ),
        opacity=0.85
    )
    data = [myhist]
    layout = go.Layout(
         xaxis=dict(
            title=xtitle
        ),
        yaxis=dict(
            title=ytitle
        )
    )
    fig = go.Figure(data=data, layout=layout)
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


# def plot_box(values:dict, title:str, xtitle:str, ytitle:str):
#     layout = go.Layout(
#          xaxis=dict(
#             title=xtitle
#         ),
#         yaxis=dict(
#             title=ytitle
#         )
#     )
#     fig = go.Figure(layout=layout)
#     fig.add_trace(go.Box(x=values))
#     fig.update_layout(
#         margin=go.layout.Margin(
#             l=50,
#             r=50,
#             b=20,
#             t=30,
#             pad=0
#         )
#     )   
#     return py.offline.plot(fig, include_plotlyjs=False, output_type='div')


def plot_dist_features_len(len_by_features, feature_names):
    """
    plot_dist_features_len function plot length by features as a distribution.
    """
    hist_data = len_by_features
    group_labels = feature_names
    colors = ['#333F44', '#37AA9C', '#94F3E4']
    # Create distplot with curve_type set to 'normal'
    fig = ff.create_distplot(hist_data, group_labels, show_hist=False, show_rug=True, colors=colors)
    fig.update_layout(
        xaxis=dict(
            title='Features length'
        ),
        margin=go.layout.Margin(
            l=50,
            r=50,
            b=20,
            t=30,
            pad=0
        )
    )
    return py.offline.plot(fig, include_plotlyjs=False, output_type='div')


def plot_insertion_in_contig(positions) :
    """
    plot_insertion_in_contig function plot introns position on a contig
    """
    hist = go.Histogram(
            x=positions,
            xbins=dict(
                start=0,
                end=100,
                size=2)
    )
    layout = go.Layout(xaxis=dict(title="% of contig length"),
                       yaxis=dict(title="Number of introns"))
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

def plot_abondance_model(df_fasta:dict) :
    """
    plot_abondance_model function plot abundance from df_fasta dataframe
    """
    fig = go.Figure()
    if 'waiting' in df_fasta.columns:
        fig.add_trace(
            go.Scatter(
                x = df_fasta.index,
                y = df_fasta['waiting'],
                mode = 'lines',
                name ='Waited abundance'
            )
        )
    fig.add_trace(
        go.Scatter(
            x = df_fasta.index,
            y = df_fasta['real'],
            mode = 'lines',
            name = 'Abundance'
        )
    )
    fig.add_trace(
        go.Scatter(
            x = df_fasta.index,                              
            y = df_fasta['norm'],
            mode = 'lines',
            name = 'Normalized abundance'
        )
    )
    fig.update_layout(
        xaxis=dict(title="Contig names"),
        yaxis=dict(title="Abundance percentage"#,
        # range=[-0.25,0.5]
        ),
        margin=go.layout.Margin(
            l=50,
            r=50,
            b=20,
            t=30,
            pad=0
        )
    )
    return py.offline.plot(fig, include_plotlyjs=False, output_type='div')


def pourcent(str_mapping:str, tot:int):
    """
    pourcent function return int from flagstat HISAT2/STAR mapping in string format 
    (only for the last 3 values : Mapped, Properly paired, Singletons)
    """
    mapping=re.sub(r'(\([a-zA-Z0-9_]*.[a-zA-Z0-9_]*%\))', r" ", str_mapping)
    val=(int(mapping)*100)/tot
    return val

def nbWithoutPourcent(str_mapping:str):
    """
    nbWithoutPourcent function format a string to remove this substring : (nb %)
    """
    val=re.sub(r'(\([a-zA-Z0-9_]*.[a-zA-Z0-9_]*%\))', r" ", str_mapping)
    return val
