#!/usr/bin/env python3
import sys
import yaml

from outbreak_tools import outbreak_tools
from outbreak_tools import outbreak_clustering

def get_alias_key(lineages_yml=None):
    tree = outbreak_tools.get_tree()
    lineage_key = outbreak_clustering.get_lineage_key(tree)
    return lineage_key

def crumbs(lin, alias_key, depth=100):
    lineage_key = alias_key
    if lin in lineage_key:
        lin = lineage_key[lin]
    else:
        lin = {'name': lin, 'alias': lin, 'parent': lin.split('.')[0]}
    return [lin['alias'], lin['name']] + \
           (crumbs(lin['parent'], lineage_key, depth-1) \
               if lin['parent'] not in ['*', lin] and depth>0 else [])
