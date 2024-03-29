#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 2019
@author: jruzicka
"""

# Librairies
import os
import shutil
import analysis as an
import script as sc


def mole_basic(l, start, end):
    """
    mole 2.5 computation of given pdb entries between the two positions, basic mode
    :param l: list od pdb ids
    :param start: starting position
    :param end: ending position
    """
    pores_trunc = l[start:end]

    for pdbid in pores_trunc:
        try:
            os.mkdir("D:\calculation\jiri\mole_pores\\" + pdbid)

            os.popen("D:\\calculation\\jiri\\Pores\\Pores.exe ""D:\calculation\jiri\mole_pores\\" + pdbid + ".json")

        except :
            print(pdbid)


def mole_tm(l, start, end):
    """
    mole 2.5 computation of given pdb entries between the two positions, TM mode
    :param l: list od pdb ids
    :param start: starting position
    :param end: ending position
    """
    pores_trunc = l[start:end]

    for pdbid in pores_trunc:
        try:
            os.mkdir("D:\calculation\jiri\mole_pores\\" + pdbid + "_TM")

            os.popen("D:\\calculation\\jiri\\Pores\\Pores.exe ""D:\calculation\jiri\mole_pores\\" + pdbid + "_TM.json")

        except FileExistsError:
            print(pdbid + "_TM")


def get_exist(path):
    """
    verifies if given path exists
    :param path: path to be verified
    :return: binary
    """
    return os.path.exists(path + "\\json")


def get_list_exist(path, l):
    """
    finds the entries for which the computation was not successful, basic mode
    :param path: given path
    :param l: list of pdb ids
    :return: tuple (number of non existing computations, its list)
    """
    does_not = 0
    pdb_list = []
    for pdbid in l:
        if not get_exist(path + pdbid):
            does_not += 1
            pdb_list.append(pdbid)
    return does_not, pdb_list


def get_list_exist_tm(path, l):
    """
    finds the entries for which the computation was not successful, TM mode
    :param path: given path
    :param l: list of pdb ids
    :return: tuple (number of non existing computations, its list)
    """
    does_not = 0
    pdb_list = []
    for pdbid in l:
        if not get_exist(path + pdbid + "_TM"):
            does_not += 1
            pdb_list.append(pdbid)
    return does_not, pdb_list


if __name__ == "__main__":
    with open('pdbid_to_mole_fromTM.txt') as f:  # get the list of pores
        my_pores_to_moleTM = f.readlines()
        for i in range(len(my_pores_to_moleTM)):
            my_pores_to_moleTM[i] = my_pores_to_moleTM[i].rstrip('\n')
    print(my_pores_to_moleTM)
    # an.generate_mole_jsons("D:\calculation\jiri\mole_pores\\", False, my_pores_to_mole)
    # an.download_all_cif(my_pores_to_mole)
    # mole_basic(my_pores_to_moleTM, 0, 38)
    with open('pdbid_to_mole.txt') as f:  # get the list of pores
        my_pores_to_mole = f.readlines()
        for i in range(len(my_pores_to_mole)):
            my_pores_to_mole[i] = my_pores_to_mole[i].rstrip('\n')
    print(my_pores_to_mole)
    # an.generate_mole_jsons("D:\calculation\jiri\mole_pores\\", False, my_pores_to_mole)
    # an.download_all_cif(my_pores_to_mole)
    # mole_basic(my_pores_to_mole, 0, 976)
    print(len(my_pores_to_mole+my_pores_to_moleTM))
    print(get_list_exist("D:\calculation\jiri\mole_pores\\", my_pores_to_mole+my_pores_to_moleTM))
    redo = get_list_exist("D:\calculation\jiri\mole_pores\\", my_pores_to_mole+my_pores_to_moleTM)[1]
