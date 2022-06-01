import yaml

def get_config(cluster_name):
    '''
    Get the conficuration of perzonalized variable definitions
    '''
    # gearshift
    if cluster_name == 'gearshift':
        with open('/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/config_gearshift.yaml', 'r') as stream:
            config = yaml.safe_load(stream)
    else:
        # Calculon
        with open('/groups/umcg-wijmenga/tmp04/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/config_calculon.yaml', 'r') as stream:
            config = yaml.safe_load(stream)
    return config