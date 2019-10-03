#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 2019
@author: jruzicka
"""

# Libraries
import analysis as an
import script as sc
import numpy as np


def get_unpr_list(d):
    """
    gets a list of uniprot ids from a dictionary pdb:uniprot
    :param d: dictionary containing uniprot ids as values
    :return: list of uniprot ids
    """
    unpr = []
    for l in [*d.values()]:
        for uniprotid in l:
            if uniprotid not in unpr:
                unpr.append(uniprotid)
    return unpr


def find_pores(l, per, abs):
    """
    detection methode of new pores
    :param l: list, pdb release
    :param per: percentage of the filter
    :param abs: absolute deviation
    :return: list of potential pores - pdb ids
    """
    pores_db = sc.check_db()  # membrane databases control
    pores_cdb = sc.get_pores_from_cdb()  # channelsdb entries
    print(str(len(pores_db[0])) + " pores from MemProtDB to add in ChannelsDB.")
    print("Pores")
    mapping_pores = sc.get_all_uniprotid(pores_db[0] + pores_cdb)  # UniProt ids from sifts, pores
    print("Release")
    mapping_rel = sc.get_all_uniprotid(l)  # UniProt ids from sifts, release
    unpr_pores = get_unpr_list(mapping_pores)  # list from dict
    unpr_rel = get_unpr_list(mapping_rel)  # list from dict
    intrsc = an.intersection_list(unpr_pores, unpr_rel)  # intersection of UniProts
    # print(intrsc)
    ptn_pores = []  # potential new pores
    for pdbid in mapping_rel:
        for uniprotid in mapping_rel[pdbid]:
            if uniprotid in intrsc:
                if pdbid not in ptn_pores:
                    ptn_pores.append(pdbid)
    print(len(ptn_pores))
    sc.download_xml(intrsc)  # download UniProt files for potential pores
    if ptn_pores is not []:
        new_pores = []
        list_data = sc.load_all(intrsc)  # loads information about TM regions
        print('Data loaded')
        up_l = sc.list_upper(ptn_pores)
        tm_comparison = sc.compare_all(list_data, up_l)  # compares TM regions
        print(tm_comparison)
        # print(len(tm_comparison))
        for pdbid in up_l:
            if pdbid not in tm_comparison:
                pass
                # print(pdbid + ': Not in UniProt page')
            elif tm_comparison[pdbid] == 'No TM':
                pass
                # print(pdbid + ': No TM')
            else:
                if tm_comparison[pdbid]:
                    new_pores.append(pdbid)
                elif not tm_comparison[pdbid]:
                    # print(pdbid + ': Filter')
                    lines = sc.search_uniprot_line(list_data, pdbid)
                    for line in lines:
                        if sc.filter_tm(line, pdbid, per, abs):  # Filter comparison
                            new_pores.append(pdbid)
                            print(pdbid + ': Filter passed directly.')
                            break
                        else:
                            if sc.compare_family_tm(line, pdbid):
                                print(pdbid + ': Filter passed by family.')
                                new_pores.append(pdbid)
                                break
                            else:
                                print(pdbid + ': Filter not passed.')
        print('New ' + str(len(new_pores)) + ' pores, of which ' +
              str(len(new_pores) - an.count_presence_lists(sc.list_lower(new_pores), pores_cdb)) +
              ' are not in ChannelsDB and ' +
              str(len(new_pores) - an.count_presence_lists(sc.list_lower(new_pores), pores_db[1])) +
              ' are not in any database.')
        return new_pores
    else:
        return 'No new pores.'


if __name__ == "__main__":
    with open('pores.txt') as f:  # get the list of pores
        pores = f.readlines()
        for i in range(len(pores)):
            pores[i] = pores[i].rstrip('\n')
    print(len(pores))
    # print(find_pores(pores))

    with open('rel2019.txt') as f:  # get the list of pores
        rel2019 = f.readlines()
        for i in range(len(rel2019)):
            rel2019[i] = rel2019[i].rstrip('\n')
    print(len(rel2019))
    """
    with open('result_filt.txt', 'w') as f:
        for j in np.arange(0.0, 0.5, 0.1):
            for k in range(0, 10, 1):
                rel_pores = find_pores(sc.list_lower(rel2019), j, k)
                print(rel_pores)
                print(str(len(rel_pores)) + ' pores with percentage ' + str(j) + ' and abs ' + str(k))
                f.write("%s\n" % (str(len(rel_pores)) + ' pores with percentage ' + str(j) + ' and abs ' + str(k)))
    """
    rel_pores = find_pores(sc.list_lower(rel2019), 0.1, 5)
    print(rel_pores)

    # print(sc.get_all_uniprotid(sc.list_lower(rel_pores)))


