#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 2019
@author: jruzicka
"""

import analysis as an
import script as sc
import batch as ba

if __name__ == "__main__":

    # Upload files
    with open('pdbid_to_mole.txt') as f:  # get the list of pores
        pdbid_from_db = f.readlines()
        for i in range(len(pdbid_from_db)):
            pdbid_from_db[i] = pdbid_from_db[i].rstrip('\n')
    print(len(pdbid_from_db))
    with open('pdbid_to_mole_fromTM.txt') as f:  # get the list of pores
        pdbid_from_tm = f.readlines()
        for i in range(len(pdbid_from_tm)):
            pdbid_from_tm[i] = pdbid_from_tm[i].rstrip('\n')
    print(len(pdbid_from_tm))
    pdbid_from_cdb = sc.get_pores_from_cdb()
    print(len(pdbid_from_cdb))

    # Change
    all_pdbid = pdbid_from_db + pdbid_from_tm + pdbid_from_cdb
    print(len(all_pdbid))
    all_pdbid.remove('1gmk')
    pdbid_remove = ba.get_list_exist("C:\Computations\Calc\FromMole\\", pdbid_from_db + pdbid_from_tm)[1]
    print(len(pdbid_remove))
    for pdbid in all_pdbid:
        if pdbid in pdbid_remove:
            all_pdbid.remove(pdbid)
    for pdbid in all_pdbid:
        if pdbid in pdbid_remove:
            all_pdbid.remove(pdbid)
    print(len(all_pdbid))
    with open('noPath.txt') as f:  # get the list of pores
        no_path = f.readlines()
        for i in range(len(no_path)):
            no_path[i] = no_path[i].rstrip('\n')
    for pdbid in all_pdbid:
        if pdbid in no_path:
            all_pdbid.remove(pdbid)
    for pdbid in all_pdbid:
        if pdbid in no_path:
            all_pdbid.remove(pdbid)
    print(len(all_pdbid))

    # Properties
    """
    an.histogram_property(all_pdbid, 'charge')
    an.histogram_property(all_pdbid, 'hydropathy')
    an.histogram_property(all_pdbid, 'hydrophobicity')
    an.histogram_property(all_pdbid, 'polarity')
    an.histogram_property(all_pdbid, 'mutability')
    an.histogram_property(all_pdbid, 'length')
    """
    charge = an.analyze_property(all_pdbid, 'charge')
    for i in range(len(charge)):
        if charge[i] == -31:
            print(all_pdbid[i])

    # Residues
    """"
    print(an.show_residues_ascending(an.analyze_residues(all_pdbid, 'general')))
    print(an.show_residues_ascending(an.analyze_residues(all_pdbid, 'bottleneck')))
    print(an.show_residues_ascending(an.get_percentage(an.analyze_residues(all_pdbid, 'general'))))
    print(an.show_residues_ascending(an.get_percentage(an.analyze_residues(all_pdbid, 'bottleneck'))))
    print(an.get_stat_res_number(all_pdbid, 'general'))
    print(an.get_stat_res_number(all_pdbid, 'bottleneck'))
    """
