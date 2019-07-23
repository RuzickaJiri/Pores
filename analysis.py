#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 2019
@author: jruzicka
"""

# Libraries
import os
import ast
import urllib.request
import json
import statistics
import operator
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

# download JSON


def get_list_from_dir(path):
    """"
    returns a list of all files in a directory
    the path specified as a parameter
    """
    new_list = os.listdir(path)
    return new_list


def get_pores_from_channelsdb(file):
    """"
    return a list of pores from channelsdb
    the pores satisfy the condition array[2] != 0
    """
    with open(file) as f:
        lines = f.readlines()
        str = lines[0]
    content = ast.literal_eval(str)
    pores = []
    for i in range(len(content)):
        if list(list(content.values())[i].values())[0][2] != 0:
            pores.append(list(content.keys())[i])
    return pores


def download_file(pdbid):
    """"
    downloads a json file from channelsdb
    the pdbid as the parameter
    """
    urllib.request.urlretrieve("https://webchem.ncbr.muni.cz/API/ChannelsDB/Download/" + pdbid + "?type=json",
                               pdbid + ".json")


def download_jsons(l):
    """"
    downloads all json files from a list
    (in this case a list of pores we want to get)
    list contains the pdbids of the pores
    """
    for name in l:
        download_file(name)


def load_json(file):
    """"
    returns a dictionary issued from the json file specified as a parameter
    """
    with open(file) as f:
        py_json = json.load(f)
    return py_json


def get_list_json_ending(l):
    """"
    changes the instances of a list, adds .json ending
    """
    str_list_json = []
    for i in range(len(l)):
        str_list_json.append(l[i] + ".json")
    return str_list_json


def load_all(l):
    """"
    loads all jsons from a list and returns a list of dictionaries (previosly json files)
    l is a list that contains the names of json files
    """
    list_json = []
    for i in range(len(l)):
        list_json.append(load_json(l[i]))
    return list_json


def get_charge(d):
    """"
    returns the charge of the pore
    parameter - dictionarz from json file
    """
    charge = d['channels']['transmembranePores'][0]['properties']['charge']
    return charge


def get_stat_charge(l):
    """"
    returns the mean of charges of all pores in the list
    parameter - list of dictionaries from json files
    """
    charges = []
    for d in l:
        charges.append(get_charge(d))
    return statistics.mean(charges), statistics.stdev(charges), min(charges), max(charges)


def get_hydrophobicity(d):
    """"
    returns the hydrophobicity of the pore
    parameter - dictionary from json file
    """
    hydrophobicity = d['channels']['transmembranePores'][0]['properties']['hydrophobicity']
    return hydrophobicity


def get_stat_hydrophobicity(l):
    """"
    returns the mean of hydrophobicities of all pores in the list
    parameter - list of dictionaries from json files
    """
    hydrophobicities = []
    for d in l:
        hydrophobicities.append(get_hydrophobicity(d))
    return statistics.mean(hydrophobicities), statistics.stdev(hydrophobicities), min(hydrophobicities), max(hydrophobicities)


def get_hydropathy(d):
    """"
    returns the hydropathy of the pore
    parameter - dictionary from json file
    """
    hydropathy = d['channels']['transmembranePores'][0]['properties']['hydropathy']
    return hydropathy


def get_stat_hydropathy(l):
    """"
    returns the mean of hydropathies of all pores in the list
    parameter - list of dictionaries from json files
    """
    hydropathies = []
    for d in l:
        hydropathies.append(get_hydropathy(d))
    return statistics.mean(hydropathies), statistics.stdev(hydropathies), min(hydropathies), max(hydropathies)


def get_polarity(d):
    """"
    returns the polarity of the pore
    parameter - dictionary from json file
    """
    polarity = d['channels']['transmembranePores'][0]['properties']['polarity']
    return polarity


def get_stat_polarity(l):
    """"
    returns the mean of polarities of all pores in the list
    parameter - list of dictionaries from json files
    """
    polarities = []
    for d in l:
        polarities.append(get_polarity(d))
    return statistics.mean(polarities), statistics.stdev(polarities), min(polarities), max(polarities)


def get_mutability(d):
    """"
    returns the mutability of the pore
    parameter - dictionary from json file
    """
    mutability = d['channels']['transmembranePores'][0]['properties']['mutability']
    return mutability


def get_stat_mutability(l):
    """"
    returns the mean of mutabilities of all pores in the list
    parameter - list of dictionaries from json files
    """
    mutabilities = []
    for d in l:
        mutabilities.append(get_mutability(d))
    return statistics.mean(mutabilities), statistics.stdev(mutabilities), min(mutabilities), max(mutabilities)


def hist_charge(l):
    """"
    returns the mean of mutabilities of all pores in the list
    parameter - list of dictionaries from json files
    """
    charges = []
    for d in l:
        charges.append(get_charge(d))
    num_bins = 50
    n, bins, patches = plt.hist(charges, num_bins, facecolor='black', alpha=0.5)
    plt.xlabel('Charge')
    plt.ylabel('Quantity')
    plt.title('Histogram of Charge')
    return plt.show()


def hist_hydropathy(l):
    """"
    returns the mean of mutabilities of all pores in the list
    parameter - list of dictionaries from json files
    """
    hydropathies = []
    for d in l:
        hydropathies.append(get_hydropathy(d))
    num_bins = 50
    n, bins, patches = plt.hist(hydropathies, num_bins, facecolor='black', alpha=0.5)
    plt.xlabel('Hydropathy')
    plt.ylabel('Quantity')
    plt.title('Histogram of Hydropathy')
    return plt.show()


def hist_hydrophobicity(l):
    """"
    returns the mean of mutabilities of all pores in the list
    parameter - list of dictionaries from json files
    """
    hydrophobicities = []
    for d in l:
        hydrophobicities.append(get_hydrophobicity(d))
    num_bins = 50
    n, bins, patches = plt.hist(hydrophobicities, num_bins, facecolor='black', alpha=0.5)
    plt.xlabel('Hydrophobicity')
    plt.ylabel('Quantity')
    plt.title('Histogram of Hydrophobicity')
    return plt.show()



def hist_polarity(l):
    """"
    returns the mean of mutabilities of all pores in the list
    parameter - list of dictionaries from json files
    """
    polarities = []
    for d in l:
        polarities.append(get_polarity(d))
    num_bins = 50
    n, bins, patches = plt.hist(polarities, num_bins, facecolor='black', alpha=0.5)
    plt.xlabel('Polarity')
    plt.ylabel('Quantity')
    plt.title('Histogram of Polarity')
    return plt.show()



def hist_mutability(l):
    """"
    returns the mean of mutabilities of all pores in the list
    parameter - list of dictionaries from json files
    """
    mutabilities = []
    for d in l:
        mutabilities.append(get_mutability(d))
    num_bins = 50
    n, bins, patches = plt.hist(mutabilities, num_bins, facecolor='black', alpha=0.5)
    plt.xlabel('Mutability')
    plt.ylabel('Quantity')
    plt.title('Histogram of Mutability')
    return plt.show()


def get_residues(d):
    residues = {}
    residues_json = d['channels']['transmembranePores'][0]['layers']['residueFlow']
    for aa in residues_json:
        aa = aa[:3]
        if aa not in residues:
            residues[aa] = 1
        else:
            residues[aa] += 1
    return residues


def sort_residues(d):
    residues = ['ala', 'arg', 'asn', 'asp', 'cys', 'glu', 'gln', 'gly', 'his', 'ile', 'leu', 'lys', 'met', 'phe', 'pro', 'ser', 'thr', 'trp', 'tyr', 'val']
    Residues = []
    for i in range(len(residues)):
        Residues.append(residues[i].upper())
    new_d = {}
    for molecule in d:
        if molecule in Residues:
            new_d[molecule] = d[molecule]
    return new_d


def get_stat_residues(l):
    all_residues = {}
    residues = []
    for d in l:
        residues.append(get_residues(d))
    for d in residues:
        for aa in d:
            aa = aa[:3]
            if aa not in all_residues:
                all_residues[aa] = d[aa]
            else:
                all_residues[aa] += d[aa]
    all_residues = sort_residues(all_residues)
    return all_residues


def get_stat_res_number(l):
    d = get_stat_residues(l)
    res_nb = 0
    for aa in d:
        res_nb += d[aa]
    return res_nb, res_nb/len(d), res_nb/len(l)


# my_pores = get_pores_from_channelsdb("Content.txt")
# download_jsons(my_pores) #works


with open('pores_no_1gmk.txt') as f:  # get the list of pores
        my_pores = f.readlines()
        for i in range(len(my_pores)):
            my_pores[i] = my_pores[i].rstrip('\n')
print(my_pores)
my_list = load_all(get_list_json_ending(my_pores))  # load a list containing all jsons in form of python dicts
print(len(my_list))
print(my_list[0])
print(my_list[1])
print(str(get_charge(my_list[0])) + '\n')
# compute properties
print(get_stat_charge(my_list))
print(get_stat_hydrophobicity(my_list))
print(get_stat_hydropathy(my_list))
print(get_stat_polarity(my_list))
print(get_stat_mutability(my_list))
print(get_residues(my_list[0]))
print(get_residues(my_list[1]))
print(get_residues(my_list[2]))
print(get_residues(my_list[3]))
print(get_residues(my_list[4]))
resid = get_stat_residues(my_list)
print(resid)
print(len(resid))
print(sorted(resid.items(), key=operator.itemgetter(1)))
print(get_stat_res_number(my_list))
hist_charge(my_list)
hist_hydropathy(my_list)
hist_hydrophobicity(my_list)
hist_polarity(my_list)
hist_mutability(my_list)



