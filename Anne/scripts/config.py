import yaml

def get_config():
    '''
    Get the conficuration of perzonalized variable definitions
    '''
    with open('config_full.yaml', 'r') as stream:
        config = yaml.safe_load(stream)
    return config