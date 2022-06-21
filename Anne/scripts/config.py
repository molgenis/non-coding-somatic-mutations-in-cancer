#!/usr/bin/env python3

# Imports
import yaml

def get_config(cluster_name):
    '''
    Get the conficuration of perzonalized variable definitions
    :param cluster_name: Name of the cluster
    :return: config:     Dictionary with as keys the name of the paths and as value the paths
    '''
    # gearshift
    if cluster_name == 'gearshift':
        with open('/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/config_gearshift.yaml', 'r') as stream:
            config = yaml.safe_load(stream)
    elif cluster_name == 'Anne':
        # Anne
        with open("D:/Hanze_Groningen/STAGE/00git/non-coding-somatic-mutations-in-cancer/Anne/scripts/config_Anne.yaml", 'r') as stream:
            config = yaml.safe_load(stream)
    else:
        # Calculon
        with open('/groups/umcg-wijmenga/tmp04/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/config_calculon.yaml', 'r') as stream:
            config = yaml.safe_load(stream)
    return config