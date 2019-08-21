#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 2019
@author: jruzicka
"""

# Libraries
import analysis as an
import script as sc


def get_unpr_list(d):
    unpr = []
    for l in [*d.values()]:
        for uniprotid in l:
            if uniprotid not in unpr:
                unpr.append(uniprotid)
    return unpr


def find_pores(l):
    pores_db = sc.check_db()
    pores_cdb = sc.get_pores_from_cdb()
    print(str(len(pores_db[0])) + " pores from MemProtDB to add in ChannelsDB.")
    print("Pores")
    mapping_pores = sc.get_all_uniprotid(pores_db[0]+pores_cdb)
    print("Release")
    mapping_rel = sc.get_all_uniprotid(l)
    unpr_pores = get_unpr_list(mapping_pores)
    unpr_rel = get_unpr_list(mapping_rel)
    intrsc = an.intersection_list(unpr_pores, unpr_rel)
    print(intrsc)
    ptn_pores = []
    for pdbid in mapping_rel:
        for uniprotid in mapping_rel[pdbid]:
            if uniprotid in intrsc:
                if pdbid not in ptn_pores:
                    ptn_pores.append(pdbid)
    print(ptn_pores)
    sc.download_xml(intrsc)
    if ptn_pores is not []:
        new_pores = []
        list_data = sc.load_all(intrsc)
        print('Data loaded')
        up_l = sc.list_upper(ptn_pores)
        tm_comparison = sc.compare_all(list_data, up_l)
        print(tm_comparison)
        print(len(tm_comparison))
        for pdbid in up_l:
            if pdbid not in tm_comparison:
                print(pdbid + ': Not in UniProt page')
            elif tm_comparison[pdbid] == 'No TM':
                print(pdbid + ': No TM')
            else:
                if tm_comparison[pdbid]:
                    new_pores.append(pdbid)
                elif not tm_comparison[pdbid]:
                    print(pdbid + ': Filter')
                    lines = sc.search_uniprot_line(list_data, pdbid)
                    for line in lines:
                        if sc.filter_tm(line, pdbid):
                            new_pores.append(pdbid)
                            break
                        else:
                            if sc.compare_family_tm(line, pdbid):
                                new_pores.append(pdbid)
                                break
                            else:
                                print(pdbid + ': Filter not passed.')
        print('New ' + str(len(new_pores)) + ' pores')
        return new_pores
    else:
        return 'No new pores'


if __name__ == "__main__":
    with open('pores.txt') as f:  # get the list of pores
        pores = f.readlines()
        for i in range(len(pores)):
            pores[i] = pores[i].rstrip('\n')
    print(len(pores))
    print(find_pores(pores))

    with open('rel2019.txt') as f:  # get the list of pores
        rel2019 = f.readlines()
        for i in range(len(rel2019)):
            rel2019[i] = rel2019[i].rstrip('\n')
    print(len(rel2019))
    print(find_pores(sc.list_lower(rel2019)))

