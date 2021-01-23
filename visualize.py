import re
import datetime as dt
import numpy as np
import pandas as pd
from path import Path
from PIL import Image
import base64
from io import BytesIO
import plotly
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from skimage import io
import onion_trees as ot
import visualize as bv
import json
import statsmodels as sm
from statsmodels.formula.api import ols
from data import STATE2ABBREV, COUNTY_CORRECTIONS

import bjorn_support as bs
import mutations as bm



def load_img(img_filepath):
    img = io.imread(img_filepath)
    pil_img = Image.fromarray(img) # PIL image object
    prefix = "data:image/png;base64,"
    with BytesIO() as stream:
        pil_img.save(stream, format="png")
        base64_string = prefix + base64.b64encode(stream.getvalue()).decode("utf-8")
    fig = go.Figure(go.Image(source=base64_string))
    fig.update_layout(margin=dict(l=0, r=0, b=0, t=0),
                    coloraxis_showscale=False, template='plotly_white', autosize=True)
    fig.update_xaxes(showticklabels=False).update_yaxes(showticklabels=False)
    return fig


def world_time_relative(data, feature, values, res, strain='B117', vocs=['B.1.1.7', 'B.1.1.70']):
    if len(values)==1:
        data.loc[:, 'weekday'] = data['date'].dt.weekday
        data.loc[:, 'date'] = data['date'] - data['weekday'] * dt.timedelta(days=1)
        results = (data.loc[(data[feature]==values[0])]
                            .drop_duplicates(subset=['date', 'strain']))
        total_samples = (data[(~data['pangolin_lineage'].isin(vocs))]
                        .groupby('date')
                        .agg(total_samples=('strain', 'nunique')))
    else:
        res = res.copy()
        res.loc[:, 'tmp'] = res['date'].str.split('-')
        res = res[res['tmp'].str.len()>=3]
        res.loc[:, 'date'] = pd.to_datetime(res['date'], errors='coerce')
        res.loc[:, 'weekday'] = res['date'].dt.weekday
        res.loc[:, 'date'] = res['date'] - res['weekday'] * dt.timedelta(days=1)
        total_samples = (res[(~res['pangolin_lineage'].isin(vocs))]
                         .groupby('date')
                         .agg(total_samples=('strain', 'nunique'))
                         .reset_index())
        results = res[(res['is_vui']==True)].drop_duplicates(subset=['date', 'strain'])
        
    b117_world_time = (results.groupby('date')
                              .agg(num_samples=('strain', 'nunique'),
                                   country_counts=('country', 
                                                    lambda x: np.unique(x, 
                                                                        return_counts=True)),
                                   divisions=('division', 'unique'),
                                   locations=('location', 'unique'))
                              .reset_index())
    b117_world_time.loc[:, 'countries'] = b117_world_time['country_counts'].apply(lambda x: list(x[0]))
    b117_world_time.loc[:, 'country_counts'] = b117_world_time['country_counts'].apply(lambda x: list(x[1]))
    b117_world_time = pd.merge(b117_world_time, total_samples, on='date', how='right')
    b117_world_time.loc[:, ['num_samples', 'total_samples']] = b117_world_time[['num_samples', 'total_samples']].fillna(0)
    first_detected = b117_world_time.loc[b117_world_time['num_samples']>0]['date'].min()
    first_countries = b117_world_time.loc[b117_world_time['date']==first_detected, 'countries'].values[0]
    b117_world_time = b117_world_time[b117_world_time['date']>=first_detected]
    b117_world_time['cum_num_samples'] = b117_world_time['num_samples'].cumsum()
    b117_world_time.loc[:, 'cum_total_samples'] = b117_world_time['total_samples'].cumsum()
    b117_world_time.loc[:, 'rel_freq'] = b117_world_time['cum_num_samples'] / b117_world_time['cum_total_samples']
    fig = go.Figure(data=go.Scatter(y=b117_world_time['rel_freq'], 
                                    x=b117_world_time['date'], 
                                    name='B.1.1.7 samples', mode='markers+lines', 
                                    line_color='rgba(220,20,60,.6)',
                                    text=b117_world_time[['num_samples', 'countries', 'country_counts',
                                                          'divisions', 'locations', 
                                                          'date']],
                                    hovertemplate="<b>Number of cases: %{text[0]}</b><br>" +
                                                  "<b>Country(s) Reported: %{text[1]}</b><br>" +
                                                  "<b>Cases Per Country: %{text[2]}</b><br>" +
                                                  "<b>State(s) Reported: %{text[3]}</b><br>" +
                                                  "<b>County(s) Reported: %{text[4]}</b><br>" +
                                                  "<b>Date: %{text[5]}</b><br>"))
    fig.add_annotation(x=first_detected, 
                       y=b117_world_time.loc[b117_world_time['date']==first_detected, 'rel_freq'].values[0],
            text=f"On Earth, {strain} 1st detected in <br> {', '.join(first_countries)} <br> on week of <br> {first_detected.date()}",
            showarrow=True,
            arrowhead=1, yshift=10, arrowsize=2, ay=-250, ax=100)
    fig.update_layout(yaxis_title=f'Relative cumulative frequency of {strain} on Earth', 
                      xaxis_title='Collection Date',
                      template='plotly_white', autosize=True)#, height=850,
    fig.update_yaxes(side = 'right')
    return fig


def world_time(data, feature, values, res, strain='B117'):
    if len(values)==1:
        results = (data.loc[(data[feature]==values[0])]
                            .drop_duplicates(subset=['date', 'strain']))
    else:
        # results = (data.groupby(['date', 'country', 'division', 'purpose_of_sequencing',
        #                          'location', 'pangolin_lineage', 'strain'])
        #                .agg(mutations=('mutation', 'unique')).reset_index())
        # results['is_vui'] = results['mutations'].apply(is_vui, args=(set(values),))
        results = (res[(res['is_vui']==True)]
                        .drop_duplicates(subset=['date', 'strain']))
    b117_world_time = (results.groupby('date')
                              .agg(num_samples=('strain', 'nunique'),
                                   country_counts=('country', 
                                                    lambda x: np.unique(x, 
                                                                        return_counts=True)),
                                   divisions=('division', 'unique'),
                                   locations=('location', 'unique'))
                              .reset_index())
    b117_world_time.loc[:, 'countries'] = b117_world_time['country_counts'].apply(lambda x: list(x[0]))
    b117_world_time.loc[:, 'country_counts'] = b117_world_time['country_counts'].apply(lambda x: list(x[1]))
    b117_world_time.loc[:, 'date'] = pd.to_datetime(b117_world_time['date'], 
                                             errors='coerce')
    b117_world_time['cum_num_samples'] = b117_world_time['num_samples'].cumsum()
    first_detected = b117_world_time['date'].min()
    first_countries = b117_world_time.loc[b117_world_time['date']==first_detected, 'countries']
    fig = go.Figure(data=go.Scatter(y=b117_world_time['cum_num_samples'], 
                                    x=b117_world_time['date'], 
                                    name='B.1.1.7 samples', mode='markers+lines', 
                                    line_color='rgba(220,20,60,.6)',
                                    text=b117_world_time[['num_samples', 'countries', 'country_counts',
                                                          'divisions', 'locations', 
                                                          'date']],
                                    hovertemplate="<b>Number of cases: %{text[0]}</b><br>" +
                                                  "<b>Country(s) Reported: %{text[1]}</b><br>" +
                                                  "<b>Cases Per Country: %{text[2]}</b><br>" +
                                                  "<b>State(s) Reported: %{text[3]}</b><br>" +
                                                  "<b>County(s) Reported: %{text[4]}</b><br>" +
                                                  "<b>Date: %{text[5]}</b><br>"))
    fig.add_annotation(x=first_detected, 
                       y=b117_world_time.loc[b117_world_time['date']==first_detected, 'cum_num_samples'].values[0],
            text=f"On Earth, {strain} 1st detected in <br> {', '.join(first_countries.values[0])} <br> on <br> {first_detected.date()}",
            showarrow=True,
            arrowhead=1, yshift=10, arrowsize=2, ay=-250, ax=100)
    fig.update_layout(yaxis_title='Global cumulative number of cases over time', 
                      xaxis_title='Collection Date',
                      template='plotly_white', autosize=True)#, height=850,
    fig.update_yaxes(side = 'right')
    return fig


def us_time_relative(data, feature, values, res, strain='B117', country='USA', vocs=['B.1.1.7', 'B.1.1.70']):
    if len(values)==1:
        data.loc[:, 'weekday'] = data['date'].dt.weekday
        data.loc[:, 'date'] = data['date'] - data['weekday'] * dt.timedelta(days=1)
        results = (data.loc[(data[feature]==values[0])& 
                             (data['country']=='United States of America')]
                            .drop_duplicates(subset=['date', 'strain']))
        total_samples = (data[(~data['pangolin_lineage'].isin(vocs))& 
                                 (data['country']=='United States of America')]
                        .groupby('date')
                        .agg(total_samples=('strain', 'nunique')))
    else:
        # results = (data.groupby(['date', 'country', 'division', 'purpose_of_sequencing',
        #                          'location', 'pangolin_lineage', 'strain'])
        #                .agg(mutations=('mutation', 'unique')).reset_index())
        # results['is_vui'] = results['mutations'].apply(is_vui, args=(set(values),))
        res = res.copy()
        res.loc[:, 'tmp'] = res['date'].str.split('-')
        res = res[res['tmp'].str.len()>=3]
        res.loc[:, 'date'] = pd.to_datetime(res['date'], errors='coerce')
        res.loc[:, 'weekday'] = res['date'].dt.weekday
        res.loc[:, 'date'] = res['date'] - res['weekday'] * dt.timedelta(days=1)
        total_samples = (res[(~res['pangolin_lineage'].isin(vocs))
                            &(res['country']=='United States of America')]
                         .groupby('date')
                         .agg(total_samples=('strain', 'nunique'))
                         .reset_index())
        results = (res[(res['is_vui']==True)
                        & (res['country']=='United States of America')]
                        .drop_duplicates(subset=['date', 'strain']))
    results['purpose_of_sequencing'] = '?'
    random = results[results['purpose_of_sequencing']=='?']
    biased = results[results['purpose_of_sequencing']!='?']
    b117_us_time = (random.groupby('date')
                           .agg(
                                num_samples=('strain', 'nunique'),
                                state_counts=('division', 
                                              lambda x: np.unique(x, 
                                                                  return_counts=True))
                                )
                           .reset_index())
    b117_us_time.loc[:, 'states'] = b117_us_time['state_counts'].apply(lambda x: list(x[0]))
    b117_us_time.loc[:, 'state_counts'] = b117_us_time['state_counts'].apply(lambda x: list(x[1]))
    b117_us_time = pd.merge(b117_us_time, total_samples, on='date', how='right')
    b117_us_time.loc[:, ['num_samples', 'total_samples']] = b117_us_time[['num_samples', 'total_samples']].fillna(0)

    sdrop_us_time = (biased.groupby('date')
                           .agg(
                                num_samples=('strain', 'nunique'),
                                state_counts=('division', 
                                              lambda x: np.unique(x, 
                                                                  return_counts=True))
                                )
                           .reset_index())
    sdrop_us_time.loc[:, 'states'] = sdrop_us_time['state_counts'].apply(lambda x: list(x[0]))
    sdrop_us_time.loc[:, 'state_counts'] = sdrop_us_time['state_counts'].apply(lambda x: list(x[1]))
    sdrop_us_time = pd.merge(sdrop_us_time, total_samples, on='date', how='right')
    sdrop_us_time.loc[:, ['num_samples', 'total_samples']] = sdrop_us_time[['num_samples', 'total_samples']].fillna(0)
    
    fig = go.Figure()
    if b117_us_time[b117_us_time['num_samples']>0].shape[0] > 0:
        first_detected = b117_us_time.loc[b117_us_time['num_samples']>0]['date'].min()
        first_states = b117_us_time.loc[b117_us_time['date']==first_detected, 'states'].values[0]
        b117_us_time = b117_us_time[b117_us_time['date']>=first_detected]
        b117_us_time.loc[:, 'cum_num_samples'] = b117_us_time['num_samples'].cumsum()
        b117_us_time.loc[:, 'cum_total_samples'] = b117_us_time['total_samples'].cumsum()
        b117_us_time.loc[:, 'rel_freq'] = b117_us_time['cum_num_samples'] / b117_us_time['cum_total_samples']
        fig.add_trace(
            go.Scatter(y=b117_us_time['rel_freq'], 
                                    x=b117_us_time['date'], 
                                    name=f'{strain} samples',
                                    mode='markers+lines', 
                                    line_color='rgba(220,20,60,.6)',
                                    text=b117_us_time[['num_samples', 'states', 
                                                       'state_counts', 'date']],
                                    hovertemplate="<b>Number of cases: %{text[0]}</b><br>" +
                                                  "<b>State(s) Reported: %{text[1]}</b><br>" +
                                                  "<b>Cases per State: %{text[2]}</b><br>" +
                                                  "<b>Date: %{text[3]}</b><br>"))
        fig.add_annotation(x=first_detected, 
                       y=b117_us_time.loc[b117_us_time['date']==first_detected, 'rel_freq'].values[0],
            text=f"In US, {strain} 1st detected in <br> {', '.join(first_states)} <br> on week of <br> {first_detected.date()}",
            showarrow=True,
            arrowhead=1, yshift=10, arrowsize=2, ay=-100)  
    if sdrop_us_time[sdrop_us_time['num_samples']>0].shape[0] > 0:
        first_detected = sdrop_us_time.loc[sdrop_us_time['num_samples']>0]['date'].min()
        first_states = sdrop_us_time.loc[sdrop_us_time['date']==first_detected, 'states'].values[0]
        sdrop_us_time = sdrop_us_time[sdrop_us_time['date']>=first_detected]
        sdrop_us_time.loc[:, 'cum_num_samples'] = sdrop_us_time['num_samples'].cumsum()
        sdrop_us_time.loc[:, 'cum_total_samples'] = sdrop_us_time['total_samples'].cumsum()
        sdrop_us_time.loc[:, 'rel_freq'] = sdrop_us_time['cum_num_samples'] / sdrop_us_time['cum_total_samples']
        fig.add_trace(
            go.Scatter(y=sdrop_us_time['rel_freq'], 
                        x=sdrop_us_time['date'], 
                        name='biased sampling <br> (more info in a later section)',
                        mode='markers+lines', 
                        line_color='rgba(30,144,255,.6)',
                        text=sdrop_us_time[['num_samples', 'states', 
                                            'state_counts', 'date']],
                        hovertemplate="<b>Number of cases: %{text[0]}</b><br>" +
                                        "<b>State(s) Reported: %{text[1]}</b><br>" +
                                        "<b>Cases per State: %{text[2]}</b><br>" +
                                        "<b>Date: %{text[3]}</b><br>"))
        fig.add_annotation(x=first_detected, 
                       y=sdrop_us_time.loc[sdrop_us_time['date']==first_detected, 'rel_freq'].values[0],
            text=f"In US, {strain} 1st detected in <br> {', '.join(first_states)} <br> on week of <br> {first_detected.date()}",
            showarrow=True,
            arrowhead=1, yshift=10, arrowsize=2, ay=-100)        
    
    fig.update_yaxes(side = 'right')
    fig.update_layout(yaxis_title=f'Relative cumulative frequency of {strain} in USA', 
                      xaxis_title='Collection Date',
                      template='plotly_white', autosize=True, showlegend=True,
                      legend=dict(
                                    yanchor="top",
                                    y=0.99,
                                    xanchor="left",
                                    x=0.01
                                ))#, height=850,
    return fig


def us_time(data, feature, values, res, strain='B117', country='USA'):
    if len(values)==1:
        results = (data.loc[(data[feature]==values[0]) & 
                             (data['country']=='United States of America')]
                          .drop_duplicates(subset=['date', 'strain']))
    else:
        # results = (data.groupby(['date', 'country', 'division', 'purpose_of_sequencing',
        #                          'location', 'pangolin_lineage', 'strain'])
        #                .agg(mutations=('mutation', 'unique')).reset_index())
        # results['is_vui'] = results['mutations'].apply(is_vui, args=(set(values),))
        results = (res[(res['is_vui']==True)
                        & (res['country']=='United States of America')]
                        .drop_duplicates(subset=['date', 'strain']))
    results['purpose_of_sequencing'] = '?'
    random = results[results['purpose_of_sequencing']=='?']
    biased = results[results['purpose_of_sequencing']!='?']
    b117_us_time = (random.groupby('date')
                           .agg(
                                num_samples=('strain', 'nunique'),
                                state_counts=('division', 
                                              lambda x: np.unique(x, 
                                                                  return_counts=True))
                                )
                           .reset_index())
    b117_us_time.loc[:, 'states'] = b117_us_time['state_counts'].apply(lambda x: list(x[0]))
    b117_us_time.loc[:, 'state_counts'] = b117_us_time['state_counts'].apply(lambda x: list(x[1]))
    b117_us_time.loc[:, 'date'] = pd.to_datetime(b117_us_time['date'], 
                                             errors='coerce')
    b117_us_time['cum_num_samples'] = b117_us_time['num_samples'].cumsum()

    sdrop_us_time = (biased.groupby('date')
                           .agg(
                                num_samples=('strain', 'nunique'),
                                state_counts=('division', 
                                              lambda x: np.unique(x, 
                                                                  return_counts=True))
                                )
                           .reset_index())
    sdrop_us_time.loc[:, 'states'] = sdrop_us_time['state_counts'].apply(lambda x: list(x[0]))
    sdrop_us_time.loc[:, 'state_counts'] = sdrop_us_time['state_counts'].apply(lambda x: list(x[1]))
    sdrop_us_time.loc[:, 'date'] = pd.to_datetime(sdrop_us_time['date'], 
                                             errors='coerce')
    sdrop_us_time['cum_num_samples'] = sdrop_us_time['num_samples'].cumsum()
    
    
    fig = go.Figure()
    if b117_us_time.shape[0] > 0:

        fig.add_trace(
            go.Scatter(y=b117_us_time['cum_num_samples'], 
                                    x=b117_us_time['date'], 
                                    name=f'{strain} samples',
                                    mode='markers+lines', 
                                    line_color='rgba(220,20,60,.6)',
                                    text=b117_us_time[['num_samples', 'states', 
                                                       'state_counts', 'date']],
                                    hovertemplate="<b>Number of cases: %{text[0]}</b><br>" +
                                                  "<b>State(s) Reported: %{text[1]}</b><br>" +
                                                  "<b>Cases per State: %{text[2]}</b><br>" +
                                                  "<b>Date: %{text[3]}</b><br>"))
        first_detected = b117_us_time['date'].min()
        first_states = b117_us_time.loc[b117_us_time['date']==first_detected, 'states']
        fig.add_annotation(x=first_detected, 
                       y=b117_us_time.loc[b117_us_time['date']==first_detected, 'cum_num_samples'].values[0],
            text=f"In US, {strain} 1st detected in <br> {', '.join(first_states.values[0])} <br> on <br> {first_detected.date()}",
            showarrow=True,
            arrowhead=1, yshift=10, arrowsize=2, ay=-100)  
    if sdrop_us_time.shape[0] > 0:
        fig.add_trace(
            go.Scatter(y=sdrop_us_time['cum_num_samples'], 
                        x=sdrop_us_time['date'], 
                        name='biased sampling <br> (see notes on sampling)',
                        mode='markers+lines', 
                        line_color='rgba(30,144,255,.6)',
                        text=sdrop_us_time[['num_samples', 'states', 
                                            'state_counts', 'date']],
                        hovertemplate="<b>Number of cases: %{text[0]}</b><br>" +
                                        "<b>State(s) Reported: %{text[1]}</b><br>" +
                                        "<b>Cases per State: %{text[2]}</b><br>" +
                                        "<b>Date: %{text[3]}</b><br>"))
        first_detected = sdrop_us_time['date'].min()
        first_states = sdrop_us_time.loc[sdrop_us_time['date']==first_detected, 'states']
        fig.add_annotation(x=first_detected, 
                       y=sdrop_us_time.loc[sdrop_us_time['date']==first_detected, 'cum_num_samples'].values[0],
            text=f"In US, {strain} 1st detected in <br> {', '.join(first_states.values[0])} <br> on <br> {first_detected.date()}",
            showarrow=True,
            arrowhead=1, yshift=10, arrowsize=2, ay=-100)        
    
    fig.update_yaxes(side = 'right')
    fig.update_layout(yaxis_title=f'Cumulative number of cases over time in {country}', 
                      xaxis_title='Collection Date',
                      template='plotly_white', autosize=True, showlegend=True,
                      legend=dict(
                                    yanchor="top",
                                    y=0.99,
                                    xanchor="left",
                                    x=0.01
                                ))#, height=850,
    return fig


def ca_time_relative(data, feature, values, res, 
                     strain='B117', state='California',
                     vocs=['B.1.1.7', 'B.1.1.70']):
    if len(values)==1:
        data.loc[:, 'weekday'] = data['date'].dt.weekday
        data.loc[:, 'date'] = data['date'] - data['weekday'] * dt.timedelta(days=1)
        results = (data.loc[(data[feature]==values[0])& 
                             (data['division']==state)]
                            .drop_duplicates(subset=['date', 'strain']))
        total_samples = (data[(~data['pangolin_lineage'].isin(vocs))& 
                                 (data['division']==state)]
                        .groupby('date')
                        .agg(total_samples=('strain', 'nunique')))
    else:
        res = res.copy()
        res.loc[:, 'tmp'] = res['date'].str.split('-')
        res = res[res['tmp'].str.len()>=3]
        res.loc[:, 'date'] = pd.to_datetime(res['date'], errors='coerce')
        res.loc[:, 'weekday'] = res['date'].dt.weekday
        res.loc[:, 'date'] = res['date'] - res['weekday'] * dt.timedelta(days=1)
        total_samples = (res[(~res['pangolin_lineage'].isin(vocs))
                            &(res['division']==state)]
                         .groupby('date')
                         .agg(total_samples=('strain', 'nunique'))
                         .reset_index())
        results = res[(res['is_vui']==True)
                    & (res['division']==state)].drop_duplicates(subset=['date', 'strain'])
    results.loc[:, 'purpose_of_sequencing'] = '?'
    random = results[results['purpose_of_sequencing']=='?']
    biased = results[results['purpose_of_sequencing']!='?']
    b117_ca_time = (random.groupby('date')
                           .agg(num_samples=('strain', 'nunique'),
                                county_counts=('location', 
                                lambda x: np.unique(x, return_counts=True)))
                           .reset_index())
    b117_ca_time.loc[:, 'counties'] = b117_ca_time['county_counts'].apply(lambda x: list(x[0]))
    b117_ca_time.loc[:, 'county_counts'] = b117_ca_time['county_counts'].apply(lambda x: list(x[1]))
#     b117_ca_time.loc[:, 'date'] = pd.to_datetime(b117_ca_time['date'], 
#                                              errors='coerce')
    b117_ca_time = pd.merge(b117_ca_time, total_samples, on='date', how='right')
    b117_ca_time.loc[:, ['num_samples', 'total_samples']] = b117_ca_time[['num_samples', 'total_samples']].fillna(0)
    sdrop_ca_time = (biased.groupby('date')
                           .agg(
                                num_samples=('strain', 'nunique'),
                                county_counts=('location', 
                                              lambda x: np.unique(x, return_counts=True))
                                )
                           .reset_index())
    sdrop_ca_time.loc[:, 'counties'] = sdrop_ca_time['county_counts'].apply(lambda x: list(x[0]))
    sdrop_ca_time.loc[:, 'county_counts'] = sdrop_ca_time['county_counts'].apply(lambda x: list(x[1]))
#     sdrop_ca_time.loc[:, 'date'] = pd.to_datetime(sdrop_ca_time['date'], errors='coerce')
    sdrop_ca_time = pd.merge(sdrop_ca_time, total_samples, on='date', how='right')
    sdrop_ca_time.loc[:, ['num_samples', 'total_samples']].fillna(0, inplace=True)
    fig = go.Figure()
    if b117_ca_time[b117_ca_time['num_samples']>0].shape[0] > 0:
        first_detected = b117_ca_time.loc[b117_ca_time['num_samples']>0]['date'].min()
        first_counties = b117_ca_time.loc[b117_ca_time['date']==first_detected, 'counties'].values[0]
        b117_ca_time = b117_ca_time[b117_ca_time['date']>=first_detected]
        b117_ca_time.loc[:, 'cum_num_samples'] = b117_ca_time['num_samples'].cumsum()
        b117_ca_time.loc[:, 'cum_total_samples'] = b117_ca_time['total_samples'].cumsum()
        b117_ca_time.loc[:, 'rel_freq'] = b117_ca_time['cum_num_samples'] / b117_ca_time['cum_total_samples']
        fig.add_trace(
            go.Scatter(y=b117_ca_time['rel_freq'], 
                                    x=b117_ca_time['date'], 
                                    name=f'{strain} samples', mode='markers+lines', 
                                    line_color='rgba(220,20,60,.6)',
                                    text=b117_ca_time[['num_samples', 'counties', 
                                                       'county_counts', 'date']],
                                    hovertemplate="<b>Number of cases: %{text[0]}</b><br>" +
                                                  "<b>County(s) Reported: %{text[1]}</b><br>" +
                                                  "<b>Cases per County: %{text[2]}</b><br>" +
                                                  "<b>Date: %{text[3]}</b><br>"))
        fig.add_annotation(x=first_detected, 
                        y=b117_ca_time.loc[b117_ca_time['date']==first_detected, 'rel_freq'].values[0],
                        text=f"In CA, {strain} 1st detected in <br> {', '.join(first_counties)} county(s) <br> on week of <br> {first_detected.date()}",
                        showarrow=True,
                        arrowhead=1, yshift=10, arrowsize=2, ay=-50)
    if sdrop_ca_time[sdrop_ca_time['num_samples']>0].shape[0] > 0:
        first_detected = sdrop_ca_time.loc[sdrop_ca_time['num_samples']>0]['date'].min()
        first_counties = sdrop_ca_time.loc[sdrop_ca_time['date']==first_detected, 'counties'].values[0]
        sdrop_ca_time = sdrop_ca_time[sdrop_ca_time['date']>=first_detected]
        sdrop_ca_time.loc[:, 'cum_num_samples'] = sdrop_ca_time['num_samples'].cumsum()
        sdrop_ca_time.loc[:, 'cum_total_samples'] = sdrop_ca_time['total_samples'].cumsum()
        sdrop_ca_time.loc[:, 'rel_freq'] = sdrop_ca_time['cum_num_samples'] / sdrop_ca_time['cum_total_samples']
        fig.add_trace(
            go.Scatter(y=sdrop_ca_time['rel_freq'], 
                        x=sdrop_ca_time['date'], 
                        name='biased sampling (read next section)', 
                        mode='markers+lines', 
                        line_color='rgba(30,144,255,.6)',
                        text=sdrop_ca_time[['num_samples', 'counties', 
                                            'county_counts', 'date']],
                        hovertemplate="<b>Number of cases: %{text[0]}</b><br>" +
                                        "<b>State(s) Reported: %{text[1]}</b><br>" +
                                        "<b>Cases per State: %{text[2]}</b><br>" +
                                        "<b>Date: %{text[3]}</b><br>"
                        )
                    )
        
        fig.add_annotation(x=first_detected, 
                        y=sdrop_ca_time.loc[sdrop_ca_time['date']==first_detected, 'rel_freq'].values[0],
                text=f"""In CA, {strain} 1st detected in <br> {', '.join(first_counties)} county(s) <br> on week of <br> {first_detected.date()}""",
                showarrow=True,
                arrowhead=1, yshift=10, arrowsize=2, ay=-50)
    fig.update_yaxes(side = 'right')
    fig.update_layout(yaxis_title=f'Relative cumulative frequency of {strain} in CA', 
                      xaxis_title='Collection Date',
                      template='plotly_white', showlegend=True,
                      legend=dict(
                                    yanchor="top",
                                    y=0.99,
                                    xanchor="left",
                                    x=0.01
                                ),
                      autosize=True#, autosize=True
                                )#, height=850,
    return fig


def ca_time(data, feature, values, res, strain='B117', state='California'):
    if len(values)==1:
        results = (data.loc[(data[feature]==values[0]) & 
                                   (data['division']==state)]
                              .drop_duplicates(subset=['date', 'strain']))
    else:
        # results = (data.groupby(['date', 'country', 'division', 
        #                          'location', 'pangolin_lineage', 'strain'])
        #                .agg(mutations=('mutation', 'unique')).reset_index())
        # results['is_vui'] = results['mutations'].apply(is_vui, args=(set(values),))
        results = res[(res['is_vui']==True)
                      &(res['division']==state)].drop_duplicates(subset=['date', 'strain'])
    results['purpose_of_sequencing'] = '?'
    random = results[results['purpose_of_sequencing']=='?']
    biased = results[results['purpose_of_sequencing']!='?']
    b117_ca_time = (random.groupby('date')
                           .agg(num_samples=('strain', 'nunique'),
                                county_counts=('location', 
                                lambda x: np.unique(x, return_counts=True)))
                           .reset_index())
    b117_ca_time.loc[:, 'counties'] = b117_ca_time['county_counts'].apply(lambda x: list(x[0]))
    b117_ca_time.loc[:, 'county_counts'] = b117_ca_time['county_counts'].apply(lambda x: list(x[1]))
    b117_ca_time.loc[:, 'date'] = pd.to_datetime(b117_ca_time['date'], 
                                             errors='coerce')
    b117_ca_time.loc[:, 'cum_num_samples'] = b117_ca_time['num_samples'].cumsum()
    sdrop_ca_time = (biased.groupby('date')
                           .agg(
                                num_samples=('strain', 'nunique'),
                                county_counts=('location', 
                                              lambda x: np.unique(x, return_counts=True))
                                )
                           .reset_index())
    sdrop_ca_time.loc[:, 'counties'] = sdrop_ca_time['county_counts'].apply(lambda x: list(x[0]))
    sdrop_ca_time.loc[:, 'county_counts'] = sdrop_ca_time['county_counts'].apply(lambda x: list(x[1]))
    sdrop_ca_time.loc[:, 'date'] = pd.to_datetime(sdrop_ca_time['date'], errors='coerce')
    sdrop_ca_time['cum_num_samples'] = sdrop_ca_time['num_samples'].cumsum()
    fig = go.Figure()
    if b117_ca_time.shape[0] > 0:
        fig.add_trace(
            go.Scatter(y=b117_ca_time['cum_num_samples'], 
                                    x=b117_ca_time['date'], 
                                    name=f'{strain} samples', mode='markers+lines', 
                                    line_color='rgba(220,20,60,.6)',
                                    text=b117_ca_time[['num_samples', 'counties', 
                                                       'county_counts', 'date']],
                                    hovertemplate="<b>Number of cases: %{text[0]}</b><br>" +
                                                  "<b>County(s) Reported: %{text[1]}</b><br>" +
                                                  "<b>Cases per County: %{text[2]}</b><br>" +
                                                  "<b>Date: %{text[3]}</b><br>"))
        first_detected = b117_ca_time['date'].min()
        first_counties = b117_ca_time.loc[b117_ca_time['date']==first_detected, 'counties']
        fig.add_annotation(x=first_detected, 
                        y=b117_ca_time.loc[b117_ca_time['date']==first_detected, 'cum_num_samples'].values[0],
                text=f"In CA, {strain} 1st detected in <br> {', '.join(first_counties.values[0])} <br> on {first_detected.date()}",
                showarrow=True,
                arrowhead=1, yshift=10, arrowsize=2, ay=-50)
    if sdrop_ca_time.shape[0] > 0:
        fig.add_trace(
            go.Scatter(y=sdrop_ca_time['cum_num_samples'], 
                        x=sdrop_ca_time['date'], 
                        name='biased sampling <br> (see notes on sampling)', 
                        mode='markers+lines', 
                        line_color='rgba(30,144,255,.6)',
                        text=sdrop_ca_time[['num_samples', 'counties', 
                                            'county_counts', 'date']],
                        hovertemplate="<b>Number of cases: %{text[0]}</b><br>" +
                                        "<b>State(s) Reported: %{text[1]}</b><br>" +
                                        "<b>Cases per State: %{text[2]}</b><br>" +
                                        "<b>Date: %{text[3]}</b><br>"
                        )
                    )
        first_detected = sdrop_ca_time['date'].min()
        first_counties = sdrop_ca_time.loc[sdrop_ca_time['date']==first_detected, 'counties']
        fig.add_annotation(x=first_detected, 
                        y=sdrop_ca_time.loc[sdrop_ca_time['date']==first_detected, 'cum_num_samples'].values[0],
                text=f"In CA, {strain} 1st detected in <br> {', '.join(first_counties.values[0])} <br> on {first_detected.date()}",
                showarrow=True,
                arrowhead=1, yshift=10, arrowsize=2, ay=-50)
    fig.update_yaxes(side = 'right')
    fig.update_layout(yaxis_title=f'Cumulative number of {strain} in CA', 
                      xaxis_title='Collection Date',
                      template='plotly_white', autosize=True, showlegend=True,
                      legend=dict(
                                    yanchor="top",
                                    y=0.99,
                                    xanchor="left",
                                    x=0.01
                                ))#, height=850,
    return fig


def strain_nt_distance(data, feature, values, strain='B117', sample_sz=250, vocs=['B.1.1.7', 'B.1.351', 'B.1.1.70']):
    clock_rate = 8e-4
    if feature=='pangolin_lineage':
        dists_df = create_lineage_data(data, feature, values, strain=strain, sample_sz=sample_sz, vocs=vocs)
    elif feature=='mutation':
        dists_df = create_distance_data(data, mutations=set(values), strain=strain, sample_sz=sample_sz, vocs=vocs)
    else:
        raise ValueError(f"Feature of type {feature} is not yet available for analysis. Aborting...")
    dists_df['num_subs'] = dists_df['mutations'].str.len() / 29904
    # ignore seqs with unexpectedly high dists
    dists_df = dists_df[dists_df['num_subs']<=0.0013]
    dists_df = dists_df[~dists_df['date'].isna()]
    dists_df.loc[:, 'date'] = pd.to_datetime(dists_df['date'], errors='coerce')
    dists_df['time'] = dists_df['date'].astype(str).apply(bv.decimal_date)
    b117_model = ols('num_subs ~ time', data=dists_df[dists_df['group']!='outgroup']).fit()
    b117_model.params['time'] = clock_rate
    b117_preds = dists_df[dists_df['group']!='outgroup'].copy()
    b117_model.params['Intercept'] = np.mean(b117_preds['num_subs'] - (clock_rate*b117_preds['time']))
    b117_preds.loc[:, 'predictions'] = b117_model.predict(b117_preds['time'])
    b117_n = int(b117_preds.shape[0] / 2)
    outgrp_model = ols('num_subs ~ time', 
                       data=dists_df[dists_df['group']=='outgroup']).fit()
    outgrp_model.params['time'] = clock_rate
    outgrp_preds = dists_df[dists_df['group']=='outgroup'].copy()
    outgrp_model.params['Intercept'] = np.mean(outgrp_preds['num_subs'] - (clock_rate*outgrp_preds['time']))
    outgrp_preds.loc[:, 'predictions'] = outgrp_model.predict(outgrp_preds['time'])
    outgrp_n = int(outgrp_preds.shape[0] / 3)
    fig = go.Figure(
        data=go.Scatter(y=dists_df[dists_df['group']==f'Lineage {strain} in US']['num_subs'], 
                                    x=dists_df[dists_df['group']==f'Lineage {strain} in US']['date'],
                                    name=f'{strain} (US)', mode='markers',
                                    hovertemplate =
                                    'Sample: %{text}',
                                    marker_color='rgba(220,20,60,.6)'))
    fig.add_trace(
        go.Scatter(y=dists_df[dists_df['group']==f'Lineage {strain}']['num_subs'], 
                   x=dists_df[dists_df['group']==f'Lineage {strain}']['date'],
                   mode='markers', marker_color='rgba(30,144,255,.6)', 
                   name=f'{strain} (non-US)'
                 ))
    # fig.add_trace(go.Scatter(y=b117_preds['predictions'], 
    #                          x=b117_preds['date'], name='OLS (B.1.1.7)', 
    #                          mode='lines', line_color='rgba(0,0,0,1.)'))
    fig.add_annotation(x=b117_preds.iloc[b117_n]['date'], 
                       y=b117_preds.iloc[b117_n]['predictions'],
            text=f"{strain} Lineage",
            showarrow=True,
            arrowhead=1, yshift=10, arrowsize=2, ay=-80)
    fig.add_trace(
        go.Scatter(y=dists_df[dists_df['group']=='outgroup']['num_subs'], 
                   x=dists_df[dists_df['group']=='outgroup']['date'],
                   mode='markers', marker_color='rgb(211,211,211, .6)', 
                   name='outgroup'
                 ))
    # fig.add_trace(go.Scatter(y=outgrp_preds['predictions'], 
    #                          x=outgrp_preds['date'], name='OLS (outgroup)', 
    #                          mode='lines', line_color='rgba(0,0,0,1.)'))
    fig.add_annotation(x=outgrp_preds.iloc[outgrp_n]['date'], 
                       y=outgrp_preds.iloc[outgrp_n]['predictions'],
            text=f"outgroup",
            showarrow=True,
            arrowhead=1, yshift=10, arrowsize=2, ay=-80)
    fig.update_layout(yaxis_title='Genetic Distance (root-to-tip)',
                      xaxis_title='Collection Date',
                      template='plotly_white', autosize=True)#, height=850,
    return fig


def create_lineage_data(data, feature, values, strain, sample_sz=250, vocs=['B.1.1.7', 'B.1.1.70', 'B.1.351']):
    data = (data.groupby(['date', 'country', 'division', 
                          'location', 'pangolin_lineage', 'strain'])
                .agg(mutations=('mutation', 'unique')).reset_index())
    first_detected = data.loc[data[feature].isin(values), 'date'].min()
    mutations = set(data.loc[(data[feature].isin(values)) 
                     &(data['date']==first_detected), 'mutations'].explode().unique())
    data['d_w'] = data['mutations'].apply(compute_similarity, args=(mutations,))
    outgroup = (data[(~data[feature].isin(values))
                    &(~data['pangolin_lineage'].isin(vocs))]
                    .nlargest(sample_sz, 'd_w')['strain']
                    .unique())
    try:
        ingroup = data.loc[(data[feature].isin(values))].sample(sample_sz)['strain'].unique()
    except:
        ingroup = data.loc[(data[feature].isin(values))]['strain'].unique()
    usgroup = data.loc[(data[feature].isin(values)) & (data['country']=='United States of America'), 'strain'].unique()
    data = data.loc[(data['strain'].isin(ingroup)) | (data['strain'].isin(outgroup)) | (data['strain'].isin(usgroup))]
    data.loc[:, 'group'] = 'nan'
    data.loc[data['strain'].isin(outgroup), 'group'] = 'outgroup'
    data.loc[(data['strain'].isin(ingroup)), 'group'] = f'Lineage {strain}'
    data.loc[(data['strain'].isin(usgroup)), 'group'] = f'Lineage {strain} in US'
    return data 


def create_distance_data(data: pd.DataFrame, mutations: set, strain: str, 
                         sample_sz: int=250, vocs: list=['B.1.1.7', 'B.1.351']):
    data = (data.groupby(['date', 'country', 'division', 
                          'location', 'pangolin_lineage', 'strain'])
                .agg(mutations=('mutation', 'unique')).reset_index())
    data['is_vui'] = data['mutations'].apply(is_vui, args=(mutations,))
    ref_muts = extract_mutations(data)
    data['d_w'] = data['mutations'].apply(compute_similarity, args=(ref_muts,))
    outgroup = (data[(data['is_vui']==False)
                &(~data['pangolin_lineage'].isin(vocs))]
                .sample(sample_sz)['strain']
                .unique())
    try:
        ingroup = data.loc[(data['is_vui']==True)].sample(sample_sz)['strain'].unique()
    except:
        ingroup = data.loc[(data['is_vui']==True)]['strain'].unique()
    usgroup = data.loc[(data['is_vui']==True) & (data['country']=='United States of America'), 'strain'].unique()
    data = data.loc[(data['strain'].isin(ingroup)) | (data['strain'].isin(outgroup)) | (data['strain'].isin(usgroup))]
    data['group'] = 'outgroup'
    data.loc[(data['strain'].isin(ingroup)), 'group'] = f'Lineage {strain}'
    data.loc[(data['strain'].isin(usgroup)), 'group'] = f'Lineage {strain} in US'
    return data

def is_vui(x, mutations: set):
    return mutations.issubset(set(x))
    
def extract_mutations(data: pd.DataFrame):
    first_detected = data.loc[data['is_vui']==True, 'date'].min()
    mutations = data.loc[(data['is_vui']==True) 
                         &(data['date']==first_detected), 'mutations'].explode().unique()
    return set(mutations)
        

def compute_similarity(x, reference_mutations: set):
    common_mutations = set(x) & reference_mutations
    return len(common_mutations)


def b117_nt_distance(gisaid_data, tree_fp, b117_meta, sample_sz=250, clock_rate=8e-4):
    # nabla_symbol = u"\u2207"
    croft_meta = pd.read_csv(b117_meta, sep='\t')
    croft_meta = croft_meta[croft_meta['Country']!='USA'].copy()
    # extract B117 samples from Emma Croft's build
    b117_meta = croft_meta[croft_meta['Pangolin Lineage']=='B.1.1.7'].sample(sample_sz)
    # extract outgroup samples from Emma Croft's build
    outgrp_meta = croft_meta[croft_meta['Pangolin Lineage']!='B.1.1.7'].sample(sample_sz)
    # extract B117 US samples from GISAID
    us_b117 = gisaid_data[(gisaid_data['country']=='United States of America')
                      & (gisaid_data['pangolin_lineage']=='B.1.1.7')].copy()
    # consolidate data and analyze
    b117_data = gisaid_data[(gisaid_data['strain'].isin(b117_meta['Strain'].unique()))
                       |(gisaid_data['strain'].isin(outgrp_meta['Strain'].unique()))
                       |(gisaid_data['strain'].isin(us_b117['strain'].unique()))].copy()
    b117_data.drop_duplicates(subset=['strain', 'pos', 'alt_codon'], inplace=True)
#     b117_data = b117_data[b117_data['gene']=='S']
    dists_df = (b117_data.groupby(['strain', 'date'])
                .agg(num_nt_subs=('strain', 'count'))
                .reset_index())
    dists_df['num_nt_subs'] = dists_df['num_nt_subs'] / 29903
    dists_df = dists_df[~dists_df['date'].isna()]           
    dists_df.loc[:, 'group'] = 'outgroup'
    dists_df.loc[dists_df['strain'].isin(b117_meta['Strain'].unique()), 'group'] = 'B.1.1.7 (non-US)'
    dists_df.loc[dists_df['strain'].isin(us_b117['strain'].unique()), 'group'] = 'B.1.1.7 (US)'
    dists_df = dists_df.loc[~((dists_df['group']=='outgroup') & (dists_df['num_nt_subs']>=0.001))]
    dists_df.loc[:, 'date'] = pd.to_datetime(dists_df['date'], errors='coerce')
    dists_df['time'] = dists_df['date'].astype(str).apply(bv.decimal_date)
    b117_model = ols('num_nt_subs ~ time', data=dists_df[dists_df['group']!='outgroup']).fit()
    b117_model.params['time'] = clock_rate
    b117_preds = dists_df[dists_df['group']!='outgroup'].copy()
    b117_model.params['Intercept'] = np.mean(b117_preds['num_nt_subs'] - (clock_rate*b117_preds['time']))
    b117_preds.loc[:, 'predictions'] = b117_model.predict(b117_preds['time'])
    b117_n = int(b117_preds.shape[0] / 2)
    outgrp_model = ols('num_nt_subs ~ time', 
                       data=dists_df[dists_df['group']=='outgroup']).fit()
    outgrp_model.params['time'] = clock_rate
    outgrp_preds = dists_df[dists_df['group']=='outgroup'].copy()
    outgrp_model.params['Intercept'] = np.mean(outgrp_preds['num_nt_subs'] - (clock_rate*outgrp_preds['time']))
    outgrp_preds.loc[:, 'predictions'] = outgrp_model.predict(outgrp_preds['time'])
    outgrp_n = int(outgrp_preds.shape[0] / 3)
    fig = go.Figure(
        data=go.Scatter(y=dists_df[dists_df['group']=='B.1.1.7 (US)']['num_nt_subs'], 
                                    x=dists_df[dists_df['group']=='B.1.1.7 (US)']['date'],
                                    name='B.1.1.7 (US)', mode='markers',
                                    text=dists_df[dists_df['group']=='B.1.1.7 (US)']['strain'],
                                    hovertemplate =
                                    'Sample: %{text}',
                                    marker_color='rgba(220,20,60,.6)'))
    fig.add_trace(
        go.Scatter(y=dists_df[dists_df['group']=='B.1.1.7 (non-US)']['num_nt_subs'], 
                   x=dists_df[dists_df['group']=='B.1.1.7 (non-US)']['date'],
                   mode='markers', marker_color='rgba(30,144,255,.6)', 
                   name='B.1.1.7 (non-US)'
                 ))
    fig.add_trace(go.Scatter(y=b117_preds['predictions'], 
                             x=b117_preds['date'], name='OLS (B.1.1.7)', 
                             mode='lines', line_color='rgba(0,0,0,1.)'))
    fig.add_annotation(x=b117_preds.iloc[b117_n]['date'], 
                       y=b117_preds.iloc[b117_n]['predictions'],
            text=f"B117 Lineage",
            showarrow=True,
            arrowhead=1, yshift=10, arrowsize=2, ay=-80)
    fig.add_trace(
        go.Scatter(y=dists_df[dists_df['group']=='outgroup']['num_nt_subs'], 
                   x=dists_df[dists_df['group']=='outgroup']['date'],
                   mode='markers', marker_color='rgb(211,211,211, .6)', 
                   name='outgroup'
                 ))
    fig.add_trace(go.Scatter(y=outgrp_preds['predictions'], 
                             x=outgrp_preds['date'], name='OLS (outgroup)', 
                             mode='lines', line_color='rgba(0,0,0,1.)'))
    fig.add_annotation(x=outgrp_preds.iloc[outgrp_n]['date'], 
                       y=outgrp_preds.iloc[outgrp_n]['predictions'],
            text=f"outgroup",
            showarrow=True,
            arrowhead=1, yshift=10, arrowsize=2, ay=-80)
    fig.update_layout(yaxis_title='Genetic Distance (root-to-tip)',
                      xaxis_title='Collection Date',
                      template='plotly_white', autosize=True)#, height=850,
    return fig


def b117_aa_distance(gisaid_data, b117_meta, sample_sz=250):
    croft_meta = pd.read_csv(b117_meta, sep='\t')
    croft_meta = croft_meta[croft_meta['Country']!='USA'].copy()
    # extract B117 samples from Emma Croft's build
    b117_meta = croft_meta[croft_meta['Pangolin Lineage']=='B.1.1.7'].sample(sample_sz)
    # extract outgroup samples from Emma Croft's build
    outgrp_meta = croft_meta[croft_meta['Pangolin Lineage']!='B.1.1.7'].sample(sample_sz)
    # extract B117 US samples from GISAID
    us_b117 = gisaid_data[(gisaid_data['country']=='United States of America')
                      & (gisaid_data['pangolin_lineage']=='B.1.1.7')].copy()
    # consolidate data and analyze
    b117_data = gisaid_data[(gisaid_data['strain'].isin(b117_meta['Strain'].unique()))
                       |(gisaid_data['strain'].isin(outgrp_meta['Strain'].unique()))
                       |(gisaid_data['strain'].isin(us_b117['strain'].unique()))]
    b117_data.loc[:, 'nonsyn'] = False
    b117_data.loc[b117_data['ref_aa']!=b117_data['alt_aa'], 
                  'nonsyn'] = True
    b117_data.loc[:, 'S_nonsyn'] = False
    b117_data.loc[(b117_data['gene']=='S') &
                  (b117_data['ref_aa']!=b117_data['alt_aa']), 
                  'S_nonsyn'] = True
    dists_df = (b117_data.groupby(['strain', 'date'])
                .agg(num_nonsyn_muts=('nonsyn', 'sum'), 
                     num_S_nonsyn_muts=('S_nonsyn', 'sum'))
                .reset_index())
    dists_df = dists_df[~dists_df['date'].isna()]           
    dists_df.loc[:, 'group'] = 'outgroup'
    dists_df.loc[dists_df['strain'].isin(b117_meta['Strain'].unique()), 'group'] = 'B.1.1.7 (non-US)'
    dists_df.loc[dists_df['strain'].isin(us_b117['strain'].unique()), 'group'] = 'B.1.1.7 (US)'
    dists_df.loc[:, 'date'] = pd.to_datetime(dists_df['date'], errors='coerce')
    dists_df.loc[:, 'month'] = dists_df['date'].dt.month
    dists_df.loc[:, 'doy'] = dists_df['date'].dt.dayofyear
    dists_df.loc[:, 'time'] = dists_df['date'].astype(int)/1e12
    dists_df = dists_df.loc[~dists_df['doy'].isna()].copy()
    b117_model = ols('num_nonsyn_muts ~ time', data=dists_df[dists_df['group']!='outgroup']).fit()
    b117_preds = dists_df[dists_df['group']!='outgroup'].copy()
    b117_preds.loc[:, 'predictions'] = b117_model.predict(b117_preds['time'])
    
    outgrp_model = ols('num_nonsyn_muts ~ time', 
                       data=dists_df[dists_df['group']=='outgroup']).fit()
    outgrp_preds = dists_df[dists_df['group']=='outgroup'].copy()
    outgrp_preds.loc[:, 'predictions'] = outgrp_model.predict(outgrp_preds['time'])
    fig = go.Figure(
        data=go.Scatter(y=dists_df[dists_df['group']=='B.1.1.7 (US)']['num_nonsyn_muts'], 
                                    x=dists_df[dists_df['group']=='B.1.1.7 (US)']['date'],
                                    name='B.1.1.7 (US)', mode='markers',
                                    marker_color='rgba(220,20,60,.6)'))
    fig.add_trace(
        go.Scatter(y=dists_df[dists_df['group']=='B.1.1.7 (non-US)']['num_nonsyn_muts'], 
                   x=dists_df[dists_df['group']=='B.1.1.7 (non-US)']['date'],
                   mode='markers', marker_color='rgba(30,144,255,.6)', 
                   name='B.1.1.7 (non-US)'
                 ))
    fig.add_trace(go.Scatter(y=b117_preds['predictions'], 
                             x=b117_preds['date'], name='OLS (B.1.1.7)', 
                             mode='lines', line_color='rgba(30,144,255,.6)'))
    fig.add_trace(
        go.Scatter(y=dists_df[dists_df['group']=='outgroup']['num_nonsyn_muts'], 
                   x=dists_df[dists_df['group']=='outgroup']['date'],
                   mode='markers', marker_color='rgba(0,0,0,.6)', 
                   name='outgroup'
                 ))
    fig.add_trace(go.Scatter(y=outgrp_preds['predictions'], 
                             x=outgrp_preds['date'], name='OLS (outgroup)', 
                             mode='lines', line_color='rgba(0,0,0,1.)'))
    fig.update_layout(yaxis_title='Amino Acid Changes (root-to-tip)',
                      xaxis_title='Collection Date',
                      template='plotly_white', autosize=True)#, height=850,
    return fig


def aa_distance(subs_fp, meta_fp, alpha=0.05):
    alab_subs = pd.read_csv(subs_fp)
    alab_subs.loc[:, 'nonsyn'] = False
    alab_subs.loc[alab_subs['ref_aa']!=alab_subs['alt_aa'], 'nonsyn'] = True
    alab_subs.loc[:, 'S_nonsyn'] = False
    alab_subs.loc[(alab_subs['gene']=='S') & (alab_subs['ref_aa']!=alab_subs['alt_aa']), 'S_nonsyn'] = True
    dists_df = (alab_subs.groupby('fasta_hdr')
                .agg(num_nonsyn_muts=('nonsyn', 'sum'), num_S_nonsyn_muts=('S_nonsyn', 'sum'))
                .reset_index())
    meta = pd.read_csv(meta_fp)
    sd_meta = meta[meta['location'].str.contains('San Diego')]
    df = pd.merge(dists_df, sd_meta, on='fasta_hdr')
    df.loc[:, 'date'] = pd.to_datetime(df['collection_date'], errors='coerce')
    df.loc[:, 'month'] = df['date'].dt.month
    df.loc[:, 'doy'] = df['date'].dt.dayofyear
    df = df.loc[~df['doy'].isna()].copy()
    model = ols('num_nonsyn_muts ~ doy', data=df).fit()
    df.loc[:, 'predict'] = model.predict(df['doy'])
    df.loc[:, 'p'] = model.outlier_test(method='fdr_bh')['fdr_bh(p)']
    df.loc[:, 'outlier'] = False
    df.loc[df['p']<alpha, 'outlier'] = True
    fig = go.Figure(
        data=go.Scatter(y=df[df['outlier']==False]['num_nonsyn_muts'], 
                        x=df[df['outlier']==False]['date'], 
                        name='samples', mode='markers', 
                        marker_color='rgba(30,144,255,.6)'))
    fig.add_trace(go.Scatter(y=df[df['outlier']==True]['num_nonsyn_muts'], 
                             x=df[df['outlier']==True]['date'],
                             mode='markers', 
                             marker_color='rgba(220,20,60,.6)', name='SoIs',
                 text=df[df['outlier']==True][['ID', 'date']],
                 hovertemplate = 
                 "<b>%{text[0]}</b><br>" +
                 "<b>%{text[1]}</b><br>"))
    fig.add_trace(go.Scatter(y=df['predict'], x=df['date'], 
                             name='OLS', mode='lines', 
                             line_color='rgba(0,0,0,1.)'))
    fig.update_layout(yaxis_title='Amino Acid Changes (root-to-tip)', xaxis_title='Collection Date',
                      template='plotly_white', autosize=True)#, height=850, width=800)
    return fig


def fetch_s_muts(series):
    muts = [m for m in series.unique() if m[0]=='S']
    return muts


def s_aa_distance(subs_fp, meta_fp, alpha=0.05):
    alab_subs = pd.read_csv(subs_fp)
    alab_subs.loc[:, 'mutation'] = alab_subs['gene']+':'+alab_subs['codon_num'].astype(str)+alab_subs['alt_aa']
    alab_subs.loc[:, 'nonsyn'] = False
    alab_subs.loc[alab_subs['ref_aa']!=alab_subs['alt_aa'], 'nonsyn'] = True
    alab_subs.loc[:, 'S_nonsyn'] = False
    alab_subs.loc[(alab_subs['gene']=='S') & (alab_subs['ref_aa']!=alab_subs['alt_aa']), 'S_nonsyn'] = True
    alab_subs = alab_subs[alab_subs['S_nonsyn']==True]
    dists_df = (alab_subs.groupby('fasta_hdr')
                .agg(num_nonsyn_muts=('nonsyn', 'sum'), 
                     num_S_nonsyn_muts=('S_nonsyn', 'sum'),
                     S_nonsyn_muts=('mutation', fetch_s_muts))
                .reset_index())
    meta = pd.read_csv(meta_fp)
    sd_meta = meta[meta['location'].str.contains('San Diego')]
    df = pd.merge(dists_df, sd_meta, on='fasta_hdr')
    df.loc[:, 'date'] = pd.to_datetime(df['collection_date'], errors='coerce')
    df.loc[:, 'month'] = df['date'].dt.month
    df.loc[:, 'doy'] = df['date'].dt.dayofyear
    df = df.loc[~df['doy'].isna()].copy()
    model = ols('num_S_nonsyn_muts ~ doy', data=df).fit()
    df.loc[:, 'predict'] = model.predict(df['doy'])
    df.loc[:, 'p'] = model.outlier_test(method='fdr_bh')['fdr_bh(p)']
    df.loc[:, 'outlier'] = False
    df.loc[df['p']<alpha, 'outlier'] = True
    fig = go.Figure(data=go.Scatter(y=df[df['outlier']==False]['num_S_nonsyn_muts'],
                                    x=df[df['outlier']==False]['date'],
                                    name='samples', mode='markers', 
                                    marker_color='rgba(30,144,255,.6)'))
    fig.add_trace(go.Scatter(y=df[df['outlier']==True]['num_S_nonsyn_muts'], 
                             x=df[df['outlier']==True]['date'],
                             mode='markers', 
                             marker_color='rgba(220,20,60,.6)', name='SoIs',
                 text=df[df['outlier']==True][['ID', 'date', 'S_nonsyn_muts']],
                 hovertemplate = 
                 "<b>%{text[0]}</b><br>" +
                 "<b>%{text[1]}</b><br>" +
                 "<b>%{text[2]}</b><br>"))
    fig.add_trace(go.Scatter(y=df['predict'], x=df['date'], 
                             name='OLS', mode='lines', 
                             line_color='rgba(0,0,0,1.)'))
    fig.update_layout(yaxis_title='Amino Acid Changes in the S protein(root-to-tip)', 
                      xaxis_title='Collection Date',
                      template='plotly_white', autosize=True)#, height=850, width=800)
    return fig


def genetic_distance(tree_fp, meta_fp, patient_zero, alpha=0.05):
    tree = ot.load_tree(tree_fp, patient_zero)
    dists = {n.name: tree.distance(n.name, patient_zero) for n in tree.get_terminals()}
    dists_df = (pd.DataFrame(index=dists.keys(), data=dists.values(), 
                      columns=['genetic_distance'])
         .reset_index()
         .rename(columns={'index': 'fasta_hdr'}))
    meta = pd.read_csv(meta_fp)
    sd_meta = meta[meta['location'].str.contains('San Diego')]
    df = pd.merge(dists_df, sd_meta, on='fasta_hdr')
    df['date'] = pd.to_datetime(df['collection_date'], errors='coerce')
    df['month'] = df['date'].dt.month
    df['doy'] = df['date'].dt.dayofyear
    df = df.loc[~df['doy'].isna()].copy()
    model = ols('genetic_distance ~ doy', data=df).fit()
    df['predict'] = model.predict(df['doy'])
    df['p'] = model.outlier_test(method='fdr_bh')['fdr_bh(p)']
    df['outlier'] = False
    df.loc[df['p']<alpha, 'outlier'] = True
    fig = go.Figure(data=go.Scatter(y=df[df['outlier']==False]['genetic_distance'], x=df[df['outlier']==False]['doy'], 
                                name='samples', mode='markers', marker_color='rgba(30,144,255,.6)'))
    fig.add_trace(go.Scatter(y=df[df['outlier']==True]['genetic_distance'], x=df[df['outlier']==True]['doy'],
                             mode='markers', marker_color='rgba(220,20,60,.6)', name='SoIs',
                 text=df[df['outlier']==True][['ID', 'date']],
                 hovertemplate = 
                 "<b>%{text[0]}</b><br>" +
                 "<b>%{text[1]}</b><br>"))
    fig.add_trace(go.Scatter(y=df['predict'], x=df['doy'], name='OLS', mode='lines', line_color='rgba(0,0,0,1.)'))
    fig.update_layout(yaxis_title='Genetic Distance (root-to-tip)', xaxis_title='Collection Date',
                      template='plotly_white', autosize=True)#, height=850, width=800)
    return fig


def b117_genetic_distance(gisaid_data, msa_fp, b117_meta, patient_zero, 
                          sample_sz=250, clock_rate=8e-4, vocs: list=['B.1.1.7', 'B.1.1.70', 'B.1.351']):
    # nabla_symbol = u"\u2207"
    croft_meta = pd.read_csv(b117_meta, sep='\t')
    croft_meta = croft_meta[croft_meta['Country']!='USA'].copy()
    # extract B117 samples from Emma Croft's build
    b117_ids = croft_meta[croft_meta['Pangolin Lineage']=='B.1.1.7'].sample(sample_sz)['Strain'].unique().tolist()
    # extract outgroup samples from Emma Croft's build
    outgrp_ids = croft_meta[~croft_meta['Pangolin Lineage'].isin(vocs)].sample(sample_sz)['Strain'].unique().tolist()
    # extract B117 US samples from GISAID
    us_ids = gisaid_data[(gisaid_data['country']=='United States of America')
                      & (gisaid_data['pangolin_lineage']=='B.1.1.7')]['strain'].unique().tolist()
    sois = us_ids+outgrp_ids+b117_ids+[patient_zero]
    tree_fp = Path('tmp/b117_seqs_aligned.fasta' + '.treefile')
    if not Path.isfile(tree_fp):
        bs.fetch_seqs(msa_fp, 'tmp/b117_seqs_aligned.fasta', sois, is_aligned=True, is_gzip=True)
        tree_fp = bs.compute_tree('tmp/b117_seqs_aligned.fasta', num_cpus=20, redo=False)
    tree = ot.load_tree(tree_fp, patient_zero)
    dists = {n.name: tree.distance(n.name, patient_zero) for n in tree.get_terminals()}
    dists_df = (pd.DataFrame(index=dists.keys(), data=dists.values(), 
                      columns=['genetic_distance'])
         .reset_index()
         .rename(columns={'index': 'strain'}))
    b117_meta = gisaid_data[(gisaid_data['strain'].isin(sois))].drop_duplicates(subset=['strain']).copy()
    dists_df = pd.merge(dists_df, b117_meta, on='strain')
    dists_df.loc[:, 'group'] = 'outgroup'
    dists_df.loc[dists_df['strain'].isin(b117_ids), 'group'] = 'B.1.1.7 (non-US)'
    dists_df.loc[dists_df['strain'].isin(us_ids), 'group'] = 'B.1.1.7 (US)'
#     dists_df = dists_df.loc[~((dists_df['group']=='outgroup') & (dists_df['num_nt_subs']>=0.001))]
    dists_df = dists_df[~dists_df['date'].isna()]
    dists_df.loc[:, 'date'] = pd.to_datetime(dists_df['date'], errors='coerce')
    dists_df['time'] = dists_df['date'].astype(str).apply(bv.decimal_date)
    b117_model = ols('genetic_distance ~ time', data=dists_df[dists_df['group']!='outgroup']).fit()
    b117_model.params['time'] = clock_rate
    b117_preds = dists_df[dists_df['group']!='outgroup'].copy()
    b117_model.params['Intercept'] = np.mean(b117_preds['genetic_distance'] - (clock_rate*b117_preds['time']))
    b117_preds.loc[:, 'predictions'] = b117_model.predict(b117_preds['time'])
    b117_n = int(b117_preds.shape[0] / 2)
    outgrp_model = ols('genetic_distance ~ time', 
                       data=dists_df[dists_df['group']=='outgroup']).fit()
    outgrp_model.params['time'] = clock_rate
    outgrp_preds = dists_df[dists_df['group']=='outgroup'].copy()
    outgrp_model.params['Intercept'] = np.mean(outgrp_preds['genetic_distance'] - (clock_rate*outgrp_preds['time']))
    outgrp_preds.loc[:, 'predictions'] = outgrp_model.predict(outgrp_preds['time'])
    outgrp_n = int(outgrp_preds.shape[0] / 2)
    fig = go.Figure(
        data=go.Scatter(y=dists_df[dists_df['group']=='B.1.1.7 (US)']['genetic_distance'], 
                                    x=dists_df[dists_df['group']=='B.1.1.7 (US)']['date'],
                                    name='B.1.1.7 (US)', mode='markers',
                                    text=dists_df[dists_df['group']=='B.1.1.7 (US)']['strain'],
                                    hovertemplate =
                                    'Sample: %{text}',
                                    marker_color='rgba(220,20,60,.6)'))
    fig.add_trace(
        go.Scatter(y=dists_df[dists_df['group']=='B.1.1.7 (non-US)']['genetic_distance'], 
                   x=dists_df[dists_df['group']=='B.1.1.7 (non-US)']['date'],
                   mode='markers', marker_color='rgba(30,144,255,.6)', 
                   name='B.1.1.7 (non-US)'
                 ))
    fig.add_trace(go.Scatter(y=b117_preds['predictions'], 
                             x=b117_preds['date'], name='OLS (B.1.1.7)', 
                             mode='lines', line_color='rgba(0,0,0,1.)'))
    fig.add_annotation(x=b117_preds.iloc[b117_n]['date'], 
                       y=b117_preds.iloc[b117_n]['predictions'],
            text=f"B117 Lineage",
            showarrow=True,
            arrowhead=1, yshift=10, arrowsize=2, ay=-80)
    fig.add_trace(
        go.Scatter(y=dists_df[dists_df['group']=='outgroup']['genetic_distance'], 
                   x=dists_df[dists_df['group']=='outgroup']['date'],
                   mode='markers', marker_color='rgb(211,211,211, .6)', 
                   name='outgroup'
                 ))
    fig.add_trace(go.Scatter(y=outgrp_preds['predictions'], 
                             x=outgrp_preds['date'], name='OLS (outgroup)', 
                             mode='lines', line_color='rgba(0,0,0,1.)'))
    fig.add_annotation(x=outgrp_preds.iloc[outgrp_n]['date'], 
                       y=outgrp_preds.iloc[outgrp_n]['predictions'],
            text=f"outgroup",
            showarrow=True,
            arrowhead=1, yshift=10, arrowsize=2, ay=-80)
    fig.update_layout(yaxis_title='Genetic Distance (root-to-tip)',
                      xaxis_title='Collection Date',
                      template='plotly_white', autosize=True,
                      legend=dict(
                                    yanchor="top",
                                    y=0.99,
                                    xanchor="left",
                                    x=0.01
                                )
                                )
    fig.update_yaxes(side = 'right')
    return fig, dists_df


def mutation_diversity(data, mutation, strain='S:L452R'):
    ref_codon = data.loc[data['mutation']==mutation, 'ref_codon'].unique()[0]
    ref_pos = ','.join(str(i) for i in data.loc[data['mutation']==mutation, 'pos'].unique())
    title1 = 'A'
#     f"A: Nucleotide mutations that led to {strain}. Reference codon: {ref_codon}. Positions: {ref_pos}"
    title2 = 'B'
#     f"B: Distribution of {strain} across lineages"
    fig = make_subplots(rows=2, cols=1, row_heights=[0.1, 0.9], 
                        subplot_titles=(title1, title2), vertical_spacing=0.05)
    codon_counts = (data.loc[data['mutation']==mutation, 'alt_codon']
                      .value_counts()
                      .to_frame()
                      .reset_index()
                      .rename(columns={'index': 'alt_codon', 'alt_codon': 'num_samples'}))
    codon_counts['pct_samples'] = codon_counts['num_samples'] / codon_counts['num_samples'].sum()
    fig.add_trace(go.Bar(
                y=codon_counts['alt_codon'], x=codon_counts['num_samples'], orientation='h',
                text=codon_counts['pct_samples'],
                textposition='inside'
            ), 
                row=1, col=1)
    lineage_counts = (data.loc[data['mutation']==mutation, 'pangolin_lineage']
                  .value_counts()
                  .to_frame()
                  .reset_index()
                  .rename(columns={'index': 'lineage', 'pangolin_lineage': 'num_samples'}))
    lineage_counts['pct_samples'] = lineage_counts['num_samples'] / lineage_counts['num_samples'].sum()
    fig.add_trace(go.Bar(
            y=lineage_counts['lineage'], x=lineage_counts['num_samples'], orientation='h',
            text=lineage_counts['pct_samples'],
            textposition='outside'
        ), 
                row=2, col=1)
#     fig.for_each_xaxis(lambda axis: axis.title.update(font=dict(color = 'blue', size=8)))
    fig.update_traces(texttemplate='%{text:.2p}')
    fig.update_yaxes(title_text="Alternative <br>codons", row=1, col=1)
    fig.update_yaxes(title_text="Lineage (pangolin)", row=2, col=1)
    fig.update_xaxes(title_text="Number of Sequences", row=2, col=1)
    fig.update_layout(title=f"Reference codon: {ref_codon}, Position: {ref_pos}", 
                      template='plotly_white', width=500, showlegend=False,
                     margin={"r":0})
    return fig


def mutation_diversity_multi(data, mutations, strain='CAVUI1'):
    ref_codons, ref_positions = {}, {}
    # for each mutation type e.g. S:N501Y
    for m in mutations:
        # fetch reference codon
        ref_codons[m] = data.loc[data['mutation']==m, 'ref_codon'].unique()[0]
        # fetch all alternative codons that lead to this amino acid change
        ref_positions[m] = ','.join(str(i) for i in data.loc[data['mutation']==m, 'pos'].unique())
    titles = [m+f' ({ref_codons[m]})' for m in mutations]
    fig = make_subplots(rows=2, cols=len(mutations), horizontal_spacing=0.1,
                        subplot_titles=titles,
                        specs=[[{} for m in range(len(mutations))],
                              [{"colspan": len(mutations)}]+[None for m in range(len(mutations)-1)]],
                        )
    res = (data.groupby(['date', 'strain', 'pangolin_lineage'])
               .agg(mutations=('mutation', 'unique'))
               .reset_index())
    res['is_vui'] = res['mutations'].apply(bv.is_vui, args=(set(mutations),))
    res = res[res['is_vui']==True]
    sois = res['strain'].unique()
    data = data.loc[data['strain'].isin(sois)]
    for i, m in enumerate(mutations):
        codon_counts = (data.loc[data['mutation']==m, 'alt_codon']
                            .value_counts()
                            .to_frame()
                            .reset_index()
                            .rename(columns={'index': 'alt_codon', 'alt_codon': 'num_samples'}))
        codon_counts['pct_samples'] = codon_counts['num_samples'] / codon_counts['num_samples'].sum()
        fig.add_trace(go.Bar(
                        y=codon_counts['alt_codon'], x=codon_counts['num_samples'], orientation='h',
                        text=codon_counts['pct_samples'],
                        textposition='inside', width=0.6
                    ), 
                        row=1, col=i+1)
    lineage_counts = (res['pangolin_lineage']
                  .value_counts()
                  .to_frame()
                  .reset_index()
                  .rename(columns={'index': 'lineage', 'pangolin_lineage': 'num_samples'}))
    lineage_counts['pct_samples'] = lineage_counts['num_samples'] / lineage_counts['num_samples'].sum()
    fig.add_trace(go.Bar(
            y=lineage_counts['lineage'], x=lineage_counts['num_samples'], orientation='h',
            text=lineage_counts['pct_samples'],
            textposition='outside'
        ), 
                row=2, col=1)
    fig.update_traces(texttemplate='%{text:.2p}')
    fig.update_xaxes(showticklabels=False, row=1, col=1)
    fig.update_yaxes(title_text="Alternative <br>codons", row=1, col=1)
    fig.update_yaxes(title_text="Lineage (pangolin)", row=2, col=1)
    fig.update_xaxes(title_text="Number of Sequences", row=2, col=1)
    fig.update_layout(template='plotly_white', showlegend=False, 
                      width=500, margin={"r":0})
    for i in range(len(mutations)):
        fig.update_xaxes(showticklabels=False, row=1, col=i+1)
    for i in fig['layout']['annotations']:
        i['font'] = dict(size=7,color='#ff0000')
    return fig


def map_by_state(data: pd.DataFrame, feature: str, values: list, states_fp: str, res, strain: str='B117'):
    with open(states_fp) as f:
        states = json.load(f)
    state_map = {x['properties']['name']: x['id'] for x in states['features']}
    total_samples_by_state = data.groupby('division').agg(total_samples=('strain', 'nunique')).reset_index()
    if len(values)==1:
        results = data.loc[(data[feature].isin(values)) 
                         & (data['country']=='United States of America')].copy()
    else:    
        # results = (data.groupby(['date', 'country', 'division', 'location', 
        #                 'pangolin_lineage', 'strain'])
        #        .agg(mutations=('mutation', 'unique'))
        #        .reset_index())
        # results['is_vui'] = results['mutations'].apply(is_vui, args=(set(values),))
        results = res[(res['is_vui']==True) & (res['country']=='United States of America')]
    results_by_state = results.groupby('division').agg(num_samples=('strain', 'nunique')).reset_index()
    results_by_state = pd.merge(total_samples_by_state, results_by_state, on='division', how='left')
    results_by_state['num_samples'].fillna(0, inplace=True)
    results_by_state['id'] = results_by_state['division'].apply(lambda x: state_map.get(x, 'unk'))
    results_by_state = results_by_state[results_by_state['num_samples']>0]
    results_by_state['log_num_samples'] = results_by_state['num_samples'].apply(lambda x: np.log(x))
    fig = px.choropleth_mapbox(results_by_state, geojson=states, 
                               locations='id', color='log_num_samples',
                               color_continuous_scale='Bluyl', center={"lat": 37.0902, "lon": -100.7129},
                               range_color=(0, results_by_state['log_num_samples'].max()),
                               mapbox_style="carto-positron", zoom=2.75,
                               opacity=0.5,
                               hover_data=['division', 'num_samples', 'total_samples'],
                               labels={'num_samples':f'Sequences with {strain}', 'total_samples': 'Total Sequences'}
                              )
    tickvals = np.linspace(0, results_by_state['log_num_samples'].max(), num=5)
    fig.update_coloraxes(colorbar=dict(showticklabels=True, tickvals=tickvals, ticks='inside',
                                       ticktext=[int(np.exp(i)) for i in tickvals],
                                       title=f"{strain} Cases <br>(logarithmic)",
                                       y=0.5, x=0))
    fig.update_layout(margin={"t":0,"l":0,"b":0}, autosize=True)
    return fig, state_map, results_by_state


def map_by_country(data: pd.DataFrame, feature: str, values: list, countries_fp, res, strain: str='B117'):
    with open(countries_fp) as f:
        countries = json.load(f)
    country_map = {x['properties']['name']: x['id'] for x in countries['features']}
    total_samples_by_country = data.groupby('country').agg(total_samples=('strain', 'nunique')).reset_index()
    if len(values)==1:
        results = data.loc[data[feature]==values[0]].copy()
    else:
        # results = (data.groupby(['date', 'country', 'division', 'location', 
        #                 'pangolin_lineage', 'strain'])
        #        .agg(mutations=('mutation', 'unique'))
        #        .reset_index())
        # results['is_vui'] = results['mutations'].apply(is_vui, args=(set(values),))
        results = res[res['is_vui']==True] 
    results_by_cntry = results.groupby('country').agg(num_samples=('strain', 'nunique')).reset_index()
    results_by_cntry = pd.merge(total_samples_by_country, results_by_cntry, on='country', how='left')
    results_by_cntry['num_samples'].fillna(0, inplace=True)
    results_by_cntry['id'] = results_by_cntry['country'].apply(lambda x: country_map.get(x, 'unk'))
    results_by_cntry = results_by_cntry[results_by_cntry['num_samples']>0]
    results_by_cntry['log_num_samples'] = results_by_cntry['num_samples'].apply(lambda x: np.log(x))
    fig = px.choropleth_mapbox(results_by_cntry, geojson=countries, 
                               locations='id', color='log_num_samples',
                               color_continuous_scale="Bluyl",
                               range_color=(0, results_by_cntry['log_num_samples'].max()),
                               mapbox_style="carto-positron", zoom=0.2,
                               opacity=0.5, center={"lat": 40.0902, "lon": 15.0},
                               hover_data=['country', 'num_samples', 'total_samples'],
                               labels={'num_samples':f'Sequences with {strain}', 'total_samples': 'Total Sequences'}
                              )
    tickvals = np.linspace(0, results_by_cntry['log_num_samples'].max(), num=5)
    fig.update_coloraxes(colorbar=dict(showticklabels=True, tickvals=tickvals, 
                                       ticks='inside',
                                       ticktext=[int(np.exp(i)) for i in tickvals],
                                       title=f"{strain} Cases <br>(logarithmic)",
                                       y=0.5, x=0))
    fig.update_layout(margin={"r": 0, "t":0,"l":0,"b":0}, autosize=True)
    return fig, country_map, results_by_cntry


def map_by_county(data: pd.DataFrame, feature: str, values: list, 
                  counties_fp: str, states_fp: str, strain: str='B117',
                  corrections: dict=COUNTY_CORRECTIONS, state2abbrev: dict=STATE2ABBREV):
    data = data[data['country']=='United States of America'].copy()
    data = process_counties(data, corrections, state2abbrev)
    with open(states_fp) as f:
        states = json.load(f)
    state_map = {x['id']: state2abbrev[x['properties']['name']] for x in states['features']}
    with open(counties_fp) as f:
        counties = json.load(f)
    counties_map = {x['properties']['NAME']+'-'+state_map[x['properties']['STATE']]: x['id'] for x in counties['features']}
    total_samples_by_county = (data.groupby('county')
                                   .agg(total_samples=('strain', 'nunique'))
                                   .reset_index())
    if len(values)==1:
        results = data.loc[(data[feature]==values[0]) 
                         & (data['country']=='United States of America')]
    else:
        results = (data.groupby(['date', 'country', 'division', 'county', 
                        'pangolin_lineage', 'strain'])
               .agg(mutations=('mutation', 'unique'))
               .reset_index())
        results['is_vui'] = results['mutations'].apply(is_vui, args=(set(values),))
        results = results[(results['is_vui']==True)
                        & (results['country']=='United States of America')]
    results_by_county = results.groupby('county').agg(num_samples=('strain', 'nunique')).reset_index()
    results_by_county = results_by_county[results_by_county['county']!='unk']
    results_by_county = pd.merge(total_samples_by_county, results_by_county, on='county', how='left')
    results_by_county['num_samples'].fillna(0, inplace=True)
    results_by_county.loc[:, 'id'] = results_by_county['county'].apply(lambda x: counties_map.get(x, 'unk'))
    results_by_county = results_by_county[results_by_county['num_samples']>0]
    results_by_county['log_num_samples'] = results_by_county['num_samples'].apply(lambda x: np.log(x))
    results_by_county = results_by_county[results_by_county['county']!='unk']
    fig = px.choropleth_mapbox(results_by_county, geojson=counties, 
                               locations='id', color='log_num_samples',
                               color_continuous_scale="Bluyl", center={"lat": 40, "lon": -120},
                               range_color=(0, results_by_county['log_num_samples'].max()),
                               mapbox_style="carto-positron", zoom=4,
                               opacity=0.5,
                               hover_data=['county', 'num_samples', 'total_samples'],
                               labels={'num_samples':f'Cases with {strain}', 'total_samples': 'Total Cases'}
                              )
    # get tick values of color bar 
    tickvals = np.linspace(0, results_by_county['log_num_samples'].max(), num=5)
    # show absolute values as tick labels, other colorbar configs
    fig.update_coloraxes(colorbar=dict(showticklabels=True, tickvals=tickvals, 
                                       ticktext=[int(np.exp(i)) for i in tickvals],
                                       title=f"{strain} Cases <br>(logarithmic)",
                                       y=0.5, x=0))
    fig.update_layout(margin={"t":0,"l":0,"b":0}, autosize=True)
    return fig, counties_map, results_by_county


def process_counties(us: pd.DataFrame, corrections: dict, state2abbrev: dict):
    for key, val in corrections.items():
        us.loc[:, 'location'] = us['location'].str.replace(key, val)
    us.loc[:, 'location'] = us['location'].str.replace(' County', '')
    us.loc[:, 'location'] = us['location'].apply(check_state)
    us.loc[:, 'county'] = us['location'] + '-' + us['division'].apply(lambda x: state2abbrev.get(x, 'unk'))
    return us


def check_state(x):
    if x[-2:].isupper():
        x = x[:-3]
    return x


def decimal_date(date,fmt="%Y-%m-%d",variable=False):
#     date = str(date)
    """ Converts calendar dates in specified format to decimal date. """
    if fmt == "":
        return date
    delimiter=re.search('[^0-9A-Za-z%]',fmt) ## search for non-alphanumeric symbols in fmt (should be field delimiter)
    delimit=None
    if delimiter is not None:
        delimit=delimiter.group()

    if variable==True: ## if date is variable - extract what is available
        if delimit is not None:
            dateL=len(date.split(delimit)) ## split date based on symbol
        else:
            dateL=1 ## no non-alphanumeric characters in date, assume dealing with an imprecise date (something like just year)

        if dateL==2:
            fmt=delimit.join(fmt.split(delimit)[:-1]) ## reduce fmt down to what's available
        elif dateL==1:
            fmt=delimit.join(fmt.split(delimit)[:-2])

    adatetime=dt.datetime.strptime(date,fmt) ## convert to datetime object
    year = adatetime.year ## get year
    boy = dt.datetime(year, 1, 1) ## get beginning of the year
    eoy = dt.datetime(year + 1, 1, 1) ## get beginning of next year
    return year + ((adatetime - boy).total_seconds() / ((eoy - boy).total_seconds())) ## return fractional year