import pandas as pd 
import plotly
from path import Path
from jinja2 import Environment, FileSystemLoader  # html template engine
from flask import url_for
import visualize as bv



def generate_voc_html(feature: str, values: list, results: dict, template_name: str='voc.html'):
    # express plots in html and JS
    mutation_diversity = ''
    # config = dict({'displaylogo': False})
    config = {'displaylogo': False,
              'scrollZoom': False,
              'modeBarButtonsToAdd':['drawline',
                                    'drawopenpath',
                                    'drawrect',
                                    'eraseshape'
                                       ],
              'modeBarButtonsToRemove': ['toggleSpikelines','hoverCompareCartesian','lasso2d']}
    # config = {'displayModeBar': False}
    if results.get('mutation_diversity', None):
        mutation_diversity = plotly.offline.plot(results['mutation_diversity'], include_plotlyjs=False, output_type='div', config=config)
    sampling_img = plotly.offline.plot(results['sampling_fig'], include_plotlyjs=False, output_type='div', config=config)
    world_time = plotly.offline.plot(results['world_time'], include_plotlyjs=False, output_type='div', config=config)
    us_time = plotly.offline.plot(results['us_time'], include_plotlyjs=False, output_type='div', config=config)
    ca_time = plotly.offline.plot(results['ca_time'], include_plotlyjs=False, output_type='div', config=config)
    world_rtime = plotly.offline.plot(results['world_rtime'], include_plotlyjs=False, output_type='div', config=config)
    us_rtime = plotly.offline.plot(results['us_rtime'], include_plotlyjs=False, output_type='div', config=config)
    ca_rtime = plotly.offline.plot(results['ca_rtime'], include_plotlyjs=False, output_type='div', config=config)
    world_map = plotly.offline.plot(results['world_map'],
                                    include_plotlyjs=False, output_type='div', config=config)
    state_map = plotly.offline.plot(results['state_map'], include_plotlyjs=False, output_type='div', config=config)
    county_map = plotly.offline.plot(results['county_map'], include_plotlyjs=False, output_type='div', config=config)
    # genetic_distance_plot = plotly.offline.plot(results['genetic_distance_plot'], include_plotlyjs=False, output_type='div')
    strain_distance_plot = plotly.offline.plot(results['strain_distance_plot'], include_plotlyjs=False, output_type='div', config=config)
    # aa_distance_plot = plotly.offline.plot(results['aa_distance_plot'], include_plotlyjs=False, output_type='div')
    # s_aa_distance_plot = plotly.offline.plot(results['s_aa_distance_plot'], include_plotlyjs=False, output_type='div')
    # generate output messages
    #TODO: expt_name, first_detected
    date = results['date']
    strain = results['strain']
    total_num = results['total_num']
    num_countries = results['num_countries']
    us_num = results['us_num']
    num_states = results['num_states']
    ca_num = results['ca_num']
    num_lineages = results.get('num_lineages', '')
    mutations = results.get('mutations', '')
    # dir containing our template
    file_loader = FileSystemLoader('templates')
    # load the environment
    env = Environment(loader=file_loader)
    # load the template
    template = env.get_template(template_name)
    # render data in our template format
    html_output = template.render(feature=feature, values=values,
                                  total_num=total_num, num_countries=num_countries,
                                  us_num=us_num, num_states=num_states, ca_num=ca_num,
                                  num_lineages=num_lineages, strain=strain, 
                                  mutations=mutations,
                                  date=date, world_time=world_time, us_time=us_time,
                                  ca_time=ca_time, world_rtime=world_rtime, 
                                  ca_rtime=ca_rtime, us_rtime=us_rtime,
                                  world_map=world_map, 
                                  state_map=state_map, county_map=county_map,
                                #   genetic_distance_plot=genetic_distance_plot, 
                                  strain_distance_plot=strain_distance_plot,
                                #   aa_distance_plot=aa_distance_plot, 
                                #   s_aa_distance_plot=s_aa_distance_plot,
                                  first_detected=results['first_detected'],
                                  sampling_img=sampling_img,
                                  mutation_diversity=mutation_diversity)
    print(f"Results for {values} embedded in HTML report")
    return html_output


def generate_voc_data(feature, values, input_params):
    results = pd.DataFrame()
    res = pd.DataFrame()
    if feature == 'mutation':
        print(f"Loading variant data...")
        gisaid_data = pd.read_csv(input_params['gisaid_data_fp'], compression='gzip')
        if len(values) > 1:
            res = (gisaid_data.groupby(['date', 'country', 'division', 
                                        'purpose_of_sequencing',
                                        'location', 'pangolin_lineage', 'strain'])
                       .agg(mutations=('mutation', 'unique')).reset_index())
            res['is_vui'] = res['mutations'].apply(bv.is_vui, args=(set(values),))
    else:
        print(f"Loading metadata...")
        gisaid_data = pd.read_csv(input_params['gisaid_meta_fp'], sep='\t', compression='gzip')
        gisaid_data.loc[gisaid_data['location'].isna(), 'location'] = 'unk'
    
    gisaid_data.loc[gisaid_data['country']=='USA', 'country'] = 'United States of America'
    print(f"Collecting input parameters...")
    date = input_params['date']
    sampling_type = input_params['sampling_type']
    sampling_img_fp = input_params['sampling_img_fp']
    msa_fp = input_params['msa_fp']
    # tree_fp = input_params['tree_fp']
    b117_meta = input_params['b117_meta']
    sample_sz = input_params['sample_sz']
    # subs_fp = input_params['subs_fp']
    # meta_fp = input_params['meta_fp']
    countries_fp = input_params['countries_fp']
    states_fp = input_params['states_fp']
    counties_fp = input_params['counties_fp']
    patient_zero = input_params['patient_zero']
    print(f"Fetching strain data...")
    strain_data = get_strain_data(gisaid_data, feature, values)
    # TEXT results
    print(f"Generating text-based results")
    results = get_text_results(strain_data, feature, values)
    results['strain'] = input_params['strain']
    results['date'] = input_params['date']
    results['sampling_fig'] = bv.load_img(sampling_img_fp)
    print(f"Generating geo-based results")
    results['state_map'], _, _ = bv.map_by_state(gisaid_data, feature, values, states_fp, res, strain=results['strain'])
    results['world_map'], _, _ = bv.map_by_country(gisaid_data, feature, values, countries_fp, res, strain=results['strain'])
    results['county_map'], _, _ = bv.map_by_county(gisaid_data, feature, values, counties_fp, states_fp, strain=results['strain'])
    # filter out records with bad dates
    gisaid_data['tmp'] = gisaid_data['date'].str.split('-')
    gisaid_data = gisaid_data[gisaid_data['tmp'].str.len()>=3]
    gisaid_data['date'] = pd.to_datetime(gisaid_data['date'], errors='coerce')
    gisaid_data = gisaid_data[(gisaid_data['date']<date)&(gisaid_data['date']>'2020-03-01')]
    if res.shape[0]!=0:
        res['tmp'] = res['date'].astype(str).str.split('-')
        res = res[res['tmp'].str.len()>=3]
        res['date'] = pd.to_datetime(res['date'], errors='coerce')
        res = res.loc[(res['date']<date)&(res['date']>'2020-03-01')]
#     gisaid_data = gisaid_data[~((gisaid_data['pangolin_lineage']=='B.1.1.7')
#                             &(gisaid_data['date'].dt.month==1)) & 
#                           (gisaid_data['date'].dt.year>=2020) &
#                           ~(gisaid_data['date']=='2020-01-01 00:00:00')]
    print(f"Generating time-based results...")
    results['world_time'] = bv.world_time(gisaid_data, feature, values, res, strain=results['strain'], sampling_type=sampling_type)
    results['us_time'] = bv.us_time(gisaid_data, feature, values, res, strain=results['strain'], sampling_type=sampling_type)
    results['ca_time'] = bv.ca_time(gisaid_data, feature, values, res, strain=results['strain'], sampling_type=sampling_type)
    results['world_rtime'] = bv.world_time_relative(gisaid_data, feature, values, res, strain=results['strain'])
    results['us_rtime'] = bv.us_time_relative(gisaid_data, feature, values, res, strain=results['strain'])
    results['ca_rtime'] = bv.ca_time_relative(gisaid_data, feature, values, res, strain=results['strain']) 
    # results['genetic_distance_plot'] = bv.genetic_distance(tree_fp, meta_fp, patient_zero)
    print(f"Generating genomic results...")
    if 'B.1.1.7' in values:
        results['strain_distance_plot'], _ = bv.b117_genetic_distance(gisaid_data, msa_fp, b117_meta, 
                                                                    patient_zero=patient_zero, sample_sz=sample_sz)
    else:
        results['strain_distance_plot'] = bv.strain_nt_distance(gisaid_data, feature, values, strain=results['strain'], sample_sz=sample_sz)
    if feature=='mutation' and len(values)==1:
        results['mutation_diversity'] = bv.mutation_diversity(gisaid_data, values[0], strain=results['strain'])
    elif feature=='mutation':
        results['mutation_diversity'] = bv.mutation_diversity_multi(gisaid_data, values, res, strain=results['strain'])
    # results['aa_distance_plot'] = bv.aa_distance(subs_fp, meta_fp)
    # results['s_aa_distance_plot'] = bv.s_aa_distance(subs_fp, meta_fp)
    print(f"Results generated on {values}...")
    return results


def get_text_results(strain_data: pd.DataFrame, feature, values):
    results = {}
    num_lineages = strain_data['pangolin_lineage'].unique().shape[0]
    date = strain_data['date'].min()
    state = strain_data[strain_data['date']==date]['division'].unique()
    cntry = strain_data[strain_data['date']==date]['country'].unique()
    results['first_detected'] = f"The {values} {feature} was first detected on {date} in {state}, {cntry}"
    results['total_num'] = strain_data['strain'].unique().shape[0]
    results['num_countries'] = strain_data['country'].unique().shape[0]
    results['us_num'] = strain_data.loc[(strain_data['country']=='United States of America'), 'strain'].unique().shape[0]
    results['num_states'] = strain_data.loc[(strain_data['country']=='United States of America'), 'division'].unique().shape[0]
    results['ca_num'] = strain_data.loc[(strain_data['division']=='California'), 'strain'].unique().shape[0]
    if feature=='mutation':
        results['num_lineages'] = f"""The {', '.join(values)} mutation(s) has been detected in {num_lineages} lineage(s) (see fig 1.1B)."""
        results['mutations'] = ', '.join(values)
    return results


def get_strain_data(data, feature, values):
    if len(values)==1:
        strain_data = data.loc[data[feature]==values[0]]
    elif feature=='mutation':
        strain_data = (data.groupby(['date', 'pangolin_lineage', 
                                     'country', 'division', 'strain'])
                           .agg(mutations=('mutation', 'unique'))
                           .reset_index())
        strain_data['is_vui'] = strain_data['mutations'].apply(bv.is_vui, args=(set(values),)) 
        strain_data = strain_data.loc[strain_data['is_vui']==True]
    else:
        strain_data = (data.groupby(['date', 'country', 
                                     'division', 'strain'])
                           .agg(lineages=('pangolin_lineage', 'unique'))
                           .reset_index())
        strain_data['is_vui'] = strain_data['lineages'].apply(bv.is_vui, args=(set(values),)) 
        strain_data = strain_data.loc[strain_data['is_vui']==True]
    return strain_data


def save_html(html_output: str, filename: str):
    with open(filename, 'w') as f:
        f.write(html_output)
    print(f"Results saved in {filename}")
    return 0