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
    print(all_pdbid)

    # Properties
    """
    an.histogram_property(all_pdbid, 'charge')
    an.histogram_property(all_pdbid, 'hydropathy')
    an.histogram_property(all_pdbid, 'hydrophobicity')
    an.histogram_property(all_pdbid, 'polarity')
    an.histogram_property(all_pdbid, 'mutability')
    an.histogram_property(all_pdbid, 'length')
    
    print(an.get_stat_property(all_pdbid, 'charge'))
    print(an.get_stat_property(all_pdbid, 'hydropathy'))
    print(an.get_stat_property(all_pdbid, 'hydrophobicity'))
    print(an.get_stat_property(all_pdbid, 'polarity'))
    print(an.get_stat_property(all_pdbid, 'mutability'))
    print(an.get_stat_property(all_pdbid, 'length'))
    """

    mpstruc = sc.get_dict_classes('mpstruc')
    chanpot = [*mpstruc['channels-potassium-sodium-proton-ion-selective']]
    chancal = [*mpstruc['channels-calcium-ion-selective']]
    chanoth = [*mpstruc['channels-other-ion-channels']]
    porins = [*mpstruc['channels-aquaporins-and-glyceroporins']]
    print(len(chanpot), len(chancal), len(chanoth), len(porins))
    list_classes = [chanpot, chancal, chanoth, porins]
    """
    pdbid_remove_classes = ba.get_list_exist("C:\Computations\Calc\FromMole\\", chanpot + chancal + chanoth + porins)[1]\
                           + ['1k4d', '3ldc', '3lde', '3ous','3stl','3zjz','4p2z','4p30','5ec2','5wuc','6c1e','6c1k', '4ycr']
    print(pdbid_remove_classes)
    
    for c in list_classes:
        for pdbid in c:
            if pdbid in pdbid_remove_classes:
                c.remove(pdbid)
        for pdbid in c:
            if pdbid in pdbid_remove_classes:
                c.remove(pdbid)
    print(len(chanpot), len(chancal), len(chanoth), len(porins))
    

    for c in list_classes:
        print(len(c))
        print(len(sc.get_data_from_cdb(c)))
        an.download_jsons(sc.get_data_from_cdb(c))

        print('charge: ' + str(an.get_stat_property(c, 'charge')))
        print('hydropathy: ' + str(an.get_stat_property(c, 'hydropathy')))
        print('hydrophobicity: ' + str(an.get_stat_property(c, 'hydrophobicity')))
        print('polarity: ' + str(an.get_stat_property(c, 'polarity')))
        print('mutability: ' + str(an.get_stat_property(c, 'mutability')))
        print('length: ' + str(an.get_stat_property(c, 'length')))
        print(an.show_residues_ascending(an.get_percentage(an.analyze_residues(c, 'general'))))
        print(an.show_residues_ascending(an.get_percentage(an.analyze_residues(c, 'bottleneck'))))
        print('')
    """
    print(len(porins))
    print(len(sc.get_data_from_cdb(porins)))
    for p in porins:
         if p not in sc.get_data_from_cdb(porins):
             print(p)
    porins.remove('2o9f')
    porins.remove('2w1p')

    print('charge: ' + str(an.get_stat_property(porins, 'charge')))
    print('hydropathy: ' + str(an.get_stat_property(porins, 'hydropathy')))
    print('hydrophobicity: ' + str(an.get_stat_property(porins, 'hydrophobicity')))
    print('polarity: ' + str(an.get_stat_property(porins, 'polarity')))
    print('mutability: ' + str(an.get_stat_property(porins, 'mutability')))
    print('length: ' + str(an.get_stat_property(porins, 'length')))
    print(an.show_residues_ascending(an.get_percentage(an.analyze_residues(porins, 'general'))))
    print(an.show_residues_ascending(an.get_percentage(an.analyze_residues(porins, 'bottleneck'))))
    print('')

    print(an.analyze_property(porins,'charge'))
    uni_porins = sc.get_all_uniprotid(porins)
    print(len(an.list_from_dict(uni_porins)))
    an.histogram_property(porins, 'charge')

    # an.download_jsons(sc.get_pores_from_cdb())

    # Residues
    """
    print(an.show_residues_ascending(an.analyze_residues(all_pdbid, 'general')))
    print(an.show_residues_ascending(an.analyze_residues(all_pdbid, 'bottleneck')))
    print(an.show_residues_ascending(an.get_percentage(an.analyze_residues(all_pdbid, 'general'))))
    print(an.show_residues_ascending(an.get_percentage(an.analyze_residues(all_pdbid, 'bottleneck'))))
    """

