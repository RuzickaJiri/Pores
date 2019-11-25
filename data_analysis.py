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

    # an.download_jsons(an.get_pores_from_channelsdb('Content.txt'))

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
    # Mem ids
    all_pdbid = pdbid_from_db + pdbid_from_tm # + pdbid_from_cdb
    print(len(all_pdbid))
    print(len(all_pdbid))
    # all_pdbid.remove('1gmk')
    pdbid_remove = ba.get_list_exist("C:\Computations\Calc\FromMoleMem\\", pdbid_from_db + pdbid_from_tm)[1] \
                   + ['2bob','2r9r','3m75','3mra','3k04','4pdv','3f5w','3or6','4pa9','4uuj', '3m76','3ous', '3m76',
                      '3ous','3k06','4r6z','5e1j','3m77','3k08','3m78','3k0d']
    print(len(pdbid_remove))
    for pdbid in all_pdbid:
        if pdbid in pdbid_remove:
            all_pdbid.remove(pdbid)
    for pdbid in all_pdbid:
        if pdbid in pdbid_remove:
            all_pdbid.remove(pdbid)
    print(len(all_pdbid))

    with open('noPathMem.txt') as f:  # get the list of pores
        no_path = f.readlines()
        for i in range(len(no_path)):
            no_path[i] = no_path[i].rstrip('\n')
    for pdbid in all_pdbid:
        if pdbid in no_path:
            all_pdbid.remove(pdbid)
    for pdbid in all_pdbid:
        if pdbid in no_path:
            all_pdbid.remove(pdbid)
    for pdbid in all_pdbid:
        if pdbid in no_path:
            all_pdbid.remove(pdbid)

    # Full structure ids
    # Change
    all_pdbid_full = pdbid_from_db + pdbid_from_tm  # + pdbid_from_cdb
    print(len(all_pdbid_full))
    # all_pdbid.remove('1gmk')
    pdbid_remove = ba.get_list_exist("C:\Computations\Calc\FromMole\\", pdbid_from_db + pdbid_from_tm)[1]
    print(len(pdbid_remove))
    for pdbid in all_pdbid_full:
        if pdbid in pdbid_remove:
            all_pdbid_full.remove(pdbid)
    for pdbid in all_pdbid_full:
        if pdbid in pdbid_remove:
            all_pdbid_full.remove(pdbid)
    print(len(all_pdbid_full))

    with open('noPath.txt') as f:  # get the list of pores
        no_path = f.readlines()
        for i in range(len(no_path)):
            no_path[i] = no_path[i].rstrip('\n')
    for pdbid in all_pdbid_full:
        if pdbid in no_path:
            all_pdbid_full.remove(pdbid)
    for pdbid in all_pdbid_full:
        if pdbid in no_path:
            all_pdbid_full.remove(pdbid)
    for pdbid in all_pdbid_full:
        if pdbid in no_path:
            all_pdbid_full.remove(pdbid)

    print(len(all_pdbid_full))
    print(all_pdbid_full)

    print(len(all_pdbid))
    print(all_pdbid)

    pdbid_intr = an.intersection_list(all_pdbid, all_pdbid_full)
    print(len(pdbid_intr))
    print('Intersection: ' + str(pdbid_intr))

    # Residues #

    print(an.get_property(an.load_json('1a0s', True), 'residues#'))
    # print(an.get_stat_property(all_pdbid_full, 'residues#', False))
    # print(an.histogram_property(all_pdbid_full, 'residues#', False))

    # Bottleneck

    print(an.get_property(an.load_json('1a0s', True), 'bottleneck'))
    print(an.get_stat_property(all_pdbid_full, 'bottleneck', False))
    print(an.histogram_property(all_pdbid_full, 'bottleneck', False))

    btn = an.analyze_property(all_pdbid, 'bottleneck', True)
    for i in range(len(btn)):
        if btn[i] < 0:
            print(all_pdbid[i])

    # Properties
    """
    an.histogram_property(all_pdbid, 'charge', False)
    an.histogram_property(all_pdbid, 'hydropathy', False)
    an.histogram_property(all_pdbid, 'hydrophobicity', False)
    an.histogram_property(all_pdbid, 'polarity', False)
    an.histogram_property(all_pdbid, 'mutability', False)
    an.histogram_property(all_pdbid, 'length', False)
    """
    print(an.histogram_property(all_pdbid, 'charge', True))
    print(an.histogram_property(all_pdbid, 'hydropathy', True))
    print(an.histogram_property(all_pdbid, 'hydrophobicity', True))
    print(an.histogram_property(all_pdbid, 'polarity', True))
    print(an.histogram_property(all_pdbid, 'mutability', True))
    print(an.histogram_property(all_pdbid, 'length', True))


    mpstruc = sc.get_dict_classes('mpstruc')
    chanpot = [*mpstruc['channels-potassium-sodium-proton-ion-selective']]
    chancal = [*mpstruc['channels-calcium-ion-selective']]
    chanoth = [*mpstruc['channels-other-ion-channels']]
    porins = [*mpstruc['channels-aquaporins-and-glyceroporins']]
    print(len(chanpot), len(chancal), len(chanoth), len(porins))
    list_classes = [chanpot, chancal, chanoth, porins]

    pdbid_remove_classes = ba.get_list_exist("C:\Computations\Calc\FromMole\\", chanpot + chancal + chanoth + porins)[1]\
                           + ['1k4d', '3ldc', '3lde', '3ous','3stl','3zjz','4p2z','4p30','5ec2','5wuc','6c1e','6c1k', '4ycr', '4ltr']
    print(pdbid_remove_classes)
    
    for c in list_classes:
        for pdbid in c:
            if pdbid in pdbid_remove_classes:
                c.remove(pdbid)
        for pdbid in c:
            if pdbid in pdbid_remove_classes:
                c.remove(pdbid)
    print(len(chanpot), len(chancal), len(chanoth), len(porins))
    
    """
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
    # porins.remove('6poj')

    """
    print('charge: ' + str(an.get_stat_property(porins, 'charge')))
    print('hydropathy: ' + str(an.get_stat_property(porins, 'hydropathy')))
    print('hydrophobicity: ' + str(an.get_stat_property(porins, 'hydrophobicity')))
    print('polarity: ' + str(an.get_stat_property(porins, 'polarity')))
    print('mutability: ' + str(an.get_stat_property(porins, 'mutability')))
    print('length: ' + str(an.get_stat_property(porins, 'length')))
    print(an.show_residues_ascending(an.get_percentage(an.analyze_residues(porins, 'general'))))
    print(an.show_residues_ascending(an.get_percentage(an.analyze_residues(porins, 'bottleneck'))))
    print('')
    """

    print(an.analyze_property(porins,'charge', False))
    uni_porins = sc.get_all_uniprotid(porins)
    print(len(an.list_from_dict(uni_porins)))
    # an.histogram_property(porins, 'charge')

    length_all = an.analyze_property(all_pdbid, 'length', False)
    for i in range(len(length_all)):
        if length_all[i] < 10 or length_all[i] > 300:
            print(all_pdbid[i], str(length_all[i]))

    # an.download_jsons(sc.get_pores_from_cdb())

    # Residues
    """
    print(an.show_residues_ascending(an.analyze_residues(all_pdbid, 'general')))
    print(an.show_residues_ascending(an.analyze_residues(all_pdbid, 'bottleneck')))
    print(an.show_residues_ascending(an.get_percentage(an.analyze_residues(all_pdbid, 'general'))))
    print(an.show_residues_ascending(an.get_percentage(an.analyze_residues(all_pdbid, 'bottleneck'))))
    """

