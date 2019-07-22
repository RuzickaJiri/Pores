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


def get_mean_charge(l):
    """"
    returns the mean of charges of all pores in the list
    parameter - list of dictionaries from json files
    """
    charges = []
    for d in l:
        charges.append(get_charge(d))
    return statistics.mean(charges)


def get_hydrophobicity(d):
    """"
    returns the hydrophobicity of the pore
    parameter - dictionary from json file
    """
    hydrophobicity = d['channels']['transmembranePores'][0]['properties']['hydrophobicity']
    return hydrophobicity


def get_mean_hydrophobicity(l):
    """"
    returns the mean of hydrophobicities of all pores in the list
    parameter - list of dictionaries from json files
    """
    hydrophobicities = []
    for d in l:
        hydrophobicities.append(get_hydrophobicity(d))
    return statistics.mean(hydrophobicities)


def get_hydropathy(d):
    """"
    returns the hydropathy of the pore
    parameter - dictionary from json file
    """
    hydropathy = d['channels']['transmembranePores'][0]['properties']['hydropathy']
    return hydropathy


def get_mean_hydropathy(l):
    """"
    returns the mean of hydropathies of all pores in the list
    parameter - list of dictionaries from json files
    """
    hydropathies = []
    for d in l:
        hydropathies.append(get_hydropathy(d))
    return statistics.mean(hydropathies)


def get_polarity(d):
    """"
    returns the polarity of the pore
    parameter - dictionary from json file
    """
    polarity = d['channels']['transmembranePores'][0]['properties']['polarity']
    return polarity


def get_mean_polarity(l):
    """"
    returns the mean of polarities of all pores in the list
    parameter - list of dictionaries from json files
    """
    polarities = []
    for d in l:
        polarities.append(get_polarity(d))
    return statistics.mean(polarities)


def get_mutability(d):
    """"
    returns the mutability of the pore
    parameter - dictionary from json file
    """
    mutability = d['channels']['transmembranePores'][0]['properties']['mutability']
    return mutability


def get_mean_mutability(l):
    """"
    returns the mean of mutabilities of all pores in the list
    parameter - list of dictionaries from json files
    """
    mutabilities = []
    for d in l:
        mutabilities.append(get_mutability(d))
    return statistics.mean(mutabilities)


# my_pores = get_pores_from_channelsdb("Content.txt")
# download_jsons(my_pores) #works


with open('pores.txt') as f:  # get the list of pores
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
print(get_mean_charge(my_list))
print(get_mean_hydrophobicity(my_list))
print(get_mean_hydropathy(my_list))
print(get_mean_polarity(my_list))
print(get_mean_mutability(my_list))
