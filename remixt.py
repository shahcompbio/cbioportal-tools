import pandas as pd


def load_cnv(sample, filename, min_ploidy=None, max_ploidy=None):
    cn_data = {}
    stats_data = []
    
    with pd.HDFStore(filename) as store:
        stats = store['stats']
        stats = stats[stats['proportion_divergent'] < 0.5]
        
        if max_ploidy:
            stats = stats[stats['ploidy'] < max_ploidy]
        
        if min_ploidy:
            stats = stats[stats['ploidy'] > min_ploidy]
        
        stats = stats.sort_values('elbo').iloc[-1]
        stats['sample'] = sample

        init_id = stats['init_id']

        cn = store[f'/solutions/solution_{init_id}/cn']
        cn['segment_length'] = cn['end'] - cn['start'] + 1
        cn['length_ratio'] = cn['length'] / cn['segment_length']

        mix = store[f'/solutions/solution_{init_id}/mix']

        stats['normal_proportion'] = mix[0]

        cn_data[sample] = cn
        stats_data.append(stats)

    stats_data = pd.DataFrame(stats_data)
    stats_data['tumour_proportion'] = 1. - stats_data['normal_proportion']

    stats_data[[
        'sample', 'ploidy', 'proportion_divergent',
        'tumour_proportion', 'proportion_divergent', 'elbo']].sort_values('sample').to_csv('remixt_stats.csv', index=False)

    stats_data[[
        'sample', 'ploidy', 'proportion_divergent',
        'tumour_proportion', 'proportion_divergent', 'elbo']].sort_values('sample')

    print(cn_data)
    print(stats_data)
    
    return cn_data, stats_data
