#!/usr/bin/env python3
import sys
import yaml

def get_alias_key(lineages_yml='lineages.yml'):
    with open(lineages_yml, 'r') as alias_key:
        lineage_key = yaml.load(alias_key, Loader=yaml.Loader)
    alias_key = dict([(lin['name'], lin['parent']) for lin in lineage_key if 'parent' in lin])
    alias_key.update([(lin['name'], lin['alias']) for lin in lineage_key if lin['name'] != lin['alias']])
    alias_key.update([r for lin in lineage_key for r in \
        [(lin['name'], lin['name'].split('.')[0]), (lin['name'].split('.')[0], lin['alias'])] \
        if (lin['name'] != lin['alias']) and len(lin['name'].split('.')) == 2 ])
    for n in range(4):
        alias_key.update([(alias, '.'.join(alias.split('.')[:-1])) for name, alias in alias_key.items() if not alias in alias_key and len(alias.split('.')) > 1])
    alias_key.update({'A.1': 'A', 'B.1': 'B'})
    return alias_key

def _crumbs(lin, alias_key):
    return [lin] + ( _crumbs(alias_key[lin], alias_key) if lin in alias_key else [] )
def crumbs(lin, alias_key):
    lin = lin.upper()
    return _crumbs(lin, alias_key) if lin in alias_key else crumbs(lin[:-1], alias_key) if len(lin.split('.')) > 1 else []
