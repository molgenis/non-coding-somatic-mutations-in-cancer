import yaml

def get_config():
    '''
    Get the conficuration of perzonalized variable definitions
    '''
    with open('/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/config.yaml', 'r') as stream:
        config = yaml.safe_load(stream)
    return config