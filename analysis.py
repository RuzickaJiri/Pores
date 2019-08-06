#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 2019
@author: jruzicka
"""

# Libraries
import os
import ast
import requests
import urllib.request
import json
import statistics
import operator
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt


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


def get_property(d, prop):
    """"
    returns the specified property of the pore
    parameter - dictionary from json file
    """
    if prop == 'length':
        return d['channels']['transmembranePores'][0]['layers']['layersInfo'][-1]['layerGeometry']['endDistance']
    else:
        return d['channels']['transmembranePores'][0]['properties'][prop]


def get_stat_property(l, prop):
    """"
    returns the given property statistics (mean, stdev, min, max) of all pores in the list
    parameter - list of dictionaries from json files, string property
    """
    properties = []
    for d in l:
        properties.append(get_property(d, prop))
    return statistics.mean(properties), statistics.stdev(properties), min(properties), max(properties)


def hist_random_mole(l):
    """"
    returns the histogram of mutabilities of all pores in the list
    parameter - list of dictionaries from json files
    """
    mutabilities = []
    for d in l:
        try:
            mutabilities.append(d['Channels']['Paths'][0]['Layers']['LayersInfo'][-1]['LayerGeometry']['EndDistance'])
        except:
            pass
    num_bins = 50
    n, bins, patches = plt.hist(mutabilities, num_bins, facecolor='black', alpha=0.5, range=(0, 200))
    plt.xlabel('Length')
    plt.ylabel('Quantity')
    plt.title('Histogram of Length')
    return plt.show()


def histogram_property(l, prop):
    """"
    returns the histogram of the given property of all pores in the list
    parameter - list of dictionaries from json files, the string property
    """
    properties = []
    for d in l:
        properties.append(get_property(d, prop))
    num_bins = 50
    n, bins, patches = plt.hist(properties, num_bins, facecolor='black', alpha=0.5)
    plt.xlabel(prop)
    plt.ylabel('Quantity')
    plt.title('Histogram of ' + prop)
    return plt.show()


def get_residues(d, place):
    """"
     returns the dictionary containing the residues as keys and its quantity as values
     parameter - dictionary from json file
    """
    residues = {}
    residues_json = []
    if place == 'general':
        residues_json = d['channels']['transmembranePores'][0]['layers']['residueFlow']  # list
    elif place == 'bottleneck':
        for i in range(len(d['channels']['transmembranePores'][0]['layers']['layersInfo'])):
            if d['channels']['transmembranePores'][0]['layers']['layersInfo'][i]['layerGeometry']['bottleneck']:
                residues_json = d['channels']['transmembranePores'][0]['layers']['layersInfo'][i]['residues']
    elif place == 'mole':
        try:
            residues_json = d['Channels']['Paths'][0]['Layers']['ResidueFlow']  # list
        except IndexError:
            pass

    backbones = []
    backbones_d = {}
    backbones_spare = []
    backbones_d_spare = {}
    for aa in residues_json:
        if aa[-8:] == 'Backbone' and aa[:-9] not in residues_json:
            backbones.append(aa)
        elif aa[-8:] == 'Backbone' and aa[:-9] in residues_json:
            backbones_spare.append(aa)
    for aa in backbones:
        aa = aa[:3]
        if aa not in backbones_d:
            backbones_d[aa] = 1
        else:
            backbones_d[aa] += 1
    for aa in backbones_spare:
        aa = aa[:3]
        if aa not in backbones_d_spare:
            backbones_d_spare[aa] = 1
        else:
            backbones_d_spare[aa] += 1

    for aa in residues_json:
        aa = aa[:3]
        if aa not in residues:
            residues[aa] = 1
        else:
            residues[aa] += 1
    residues['MNC'] = 0
    for aa in residues:
        if aa in backbones_d and aa != 'GLY':
            residues[aa] -= backbones_d[aa]
            if residues['MNC'] == 0:
                residues['MNC'] = backbones_d[aa]
            else:
                residues['MNC'] += backbones_d[aa]
        if aa in backbones_d_spare and aa != 'GLY':
            residues[aa] -= backbones_d_spare[aa]

    return residues


def get_stat_residues(l, place):
    """"
    returns the dictionary containing the residue quantity of all pores in the list
    parameter - list of dictionaries from json files
    """
    all_residues = {}
    residues = []
    for d in l:
        residues.append(get_residues(d, place))
    for d in residues:
        for aa in d:
            aa = aa[:3]
            if aa not in all_residues:
                all_residues[aa] = d[aa]
            else:
                all_residues[aa] += d[aa]
    all_residues = sort_residues(all_residues)
    return all_residues


def get_stat_res_number(l, place):
    """"
    returns the residue statistics (total_residues, mean_by_residue, mean_by_pore) of all pores in the list
    parameter - list of dictionaries from json files
    """
    d = get_stat_residues(l, place)
    res_nb = 0
    for aa in d:
        res_nb += d[aa]
    return res_nb, res_nb/len(d), res_nb/len(l)


def sort_residues(d):
    """"
    returns the dictionary containing the sorted residues as keys and its quantity as values
    sorts the residues by eliminating the het molecules
    parameter - dictionary containing the unsorted residues
    """
    residues = ['ala', 'arg', 'asn', 'asp', 'cys', 'glu', 'gln', 'gly', 'his', 'ile', 'leu', 'lys', 'met', 'phe', 'pro', 'ser', 'thr', 'trp', 'tyr', 'val', 'mnc']
    Residues = []
    for i in range(len(residues)):
        Residues.append(residues[i].upper())
    new_d = {}
    for molecule in d:
        if molecule in Residues:
            new_d[molecule] = d[molecule]
    return new_d


def show_residues_ascending(d):
    """"
    sorts a dictionary by ascending value
    parameter - dictionary with residues
    """
    return sorted(d.items(), key=operator.itemgetter(1))


def average_d(d):
    """"
    returns the residue propensity i.e. divides its quantity by the total number of residues
    parameter - dictionary from json files
    """
    ave_dict = {}
    nb = 0
    for aa in d:
        nb += d[aa]
    for key in d:
        ave_dict[key] = round(d[key]/(nb/len(d)), 3)
    return ave_dict


def get_percentage(d):
    new_d = {}
    ave_d = average_d(d)
    for aa in ave_d:
        new_d[aa] = round((100/len(ave_d))*ave_d[aa], 2)
    return new_d


def get_all_memprotmd_references(database):
    """
    commands all the references from the mpm database from the MemProtMD site
    :return: json file (list)
    """
    MEMPROTMD_ROOT_URI = "http://memprotmd.bioch.ox.ac.uk/"
    return requests.post(MEMPROTMD_ROOT_URI + "api/references/all/" + database).json()


def get_dict_classes(database):

    new_d = {}
    l = get_all_memprotmd_references(database)
    for d in l:
        new_d[d['accession']] = []  # 'accession/title'
        for item in d['simulations']:
            new_d[d['accession']].append(item[:4])  # 'accession/title'
    return new_d


def get_pores_tcdb(d):
    new_d = {}
    for key in d:
        if key[0] == '1':
            new_d[key] = d[key]
    return new_d


def get_database_classes(l, database):
    """
    creates a text file with all the references from mpm
    :param l: list of references
    :return: None
    """
    with open(database + '.txt', 'w') as f:
        for d in l:
            f.write("%s\n" % d['accession'])  # 'accession/title'
            for item in d['simulations']:
                f.write("%s\n" % item[:4])


def find_uniprot_from_sifts(pdbid):
    """
    finds the uniprot id from pdb id thanks to sifts API
    :param pdbid: pdb id
    :return: uniprot id
    """
    urllib.request.urlretrieve("http://www.ebi.ac.uk/pdbe/api/mappings/uniprot/" + pdbid,  pdbid + "up.json")
    d = load_json(pdbid + "up.json")
    id = list(d[pdbid]['UniProt'].keys())[0]
    return id


def get_all_uniprotid(l):
    """
    gets all uniprot ids from a list
    :param l: list of pdb ids
    :return: list of uniprot ids
    """
    uniprotid = []
    i = 0
    for pdbid in l:
        try:
            uniprotid.append(find_uniprot_from_sifts(pdbid))
        except:
            i += 1
    print(i)
    return uniprotid


def download_uniprot_fasta(uniprotid, path):
    """"
    downloads the Uniprot fasta file from uniprot id
    the uniprot id as the parameter
    """
    urllib.request.urlretrieve("https://www.uniprot.org/uniprot/" + uniprotid + ".fasta", path + uniprotid + ".fasta")


def download_fasta(l, path):
    """"
    downloads all xml files from a list
    (in this case a list of pores we want to get)
    list contains the pdbids of the pores
    """
    for name in l:
        download_uniprot_fasta(name, path)


def count_presence(file, l):
    """
    returns the number of equalities between file and list
    :param file: file of structures (pdbid)
    :param l: list of pores (pdbid)
    :return: integer k
    """
    with open(file) as f:
        c = f.readlines()
        for i in range(len(c)):
            c[i] = c[i].rstrip('\n')
    k = 0
    for item in c:
        if item in l:
            k += 1
    return k


def get_res_fasta(file):
    """
    returns a string containing residues from a fasta file
    """
    with open(file) as f:
        c = f.readlines()
        for i in range(len(c)):
            c[i] = c[i].rstrip('\n')
        del c[0]
    residues = ""
    for line in c:
        residues += line
    return residues


def get_res_all_tm(l):
    """
    gets the total of all TM residues
    :param l: list of fasta files
    :return: dictionary of all residues in all membrane proteins
    """
    all_residues = {}
    residues = [] # list of strings
    for fasta in l:
        residues.append(get_res_fasta('C:/Users/Jirkův NB/Documents/CEITEC/fasta/membrane_proteins/' + fasta))
    for s in residues:
        for char in s:
            if char not in all_residues:
                all_residues[char] = 1
            else:
                all_residues[char] += 1
    return all_residues


def get_stat_res_number_all_tm(l):
    """"
    returns the residue statistics (total_residues, mean_by_residue, mean_by_pore) of all pores in the list,
    in this case all TM proteins
    parameter - list of dictionaries from json files
    """
    d = get_res_all_tm(l)
    res_nb = 0
    for aa in d:
        res_nb += d[aa]
    return res_nb, res_nb/len(d), res_nb/len(l)


def download_cif(pdbid):
    """"
    downloads the PDB cif file from pdb id
    the pdb id as the parameter
    """
    urllib.request.urlretrieve("https://files.rcsb.org/download/" + pdbid + ".cif",
                               pdbid + ".cif")


def download_all_cif(l):
    """"
    downloads all cif files from a list
    (in this case a list of pores we want to get)
    list contains the pdbids of the pores
    """
    for name in l:
        download_cif(name)


def change_template(tm, pdbid):
    if tm:
        data = load_json('template_TM.json')
    else:
        data = load_json('template.json')
    data["PdbId"] = pdbid
    if tm:
        data["WorkingDirectory"] = "C:\\Computations\\Calc\\" + pdbid + '_TM'
    else:
        data["WorkingDirectory"] = "C:\\Computations\\Calc\\" + pdbid

    if tm:
        with open("C:/Computations/Calc/" + pdbid + "_TM.json", "w") as f:
            json.dump(data, f)
    else:
        with open("C:/Computations/Calc/" + pdbid + ".json", "w") as f:
            json.dump(data, f)


def generate_mole_jsons(template, l):
    for pdbid in l:
        change_template(template, pdbid)


def load_all_basic(l):
    my_list = []
    for pdbid in l:
        if os.path.exists("C:/Computations/Calc/" + pdbid + "\\json"):
            my_path = "C:/Computations/Calc/" + pdbid + "/json/"
            with open(my_path + 'data.json') as f:
                py_json = json.load(f)
            my_list.append(py_json)
    return my_list


def load_all_tm(l):
    my_list = []
    for pdbid in l:
        if os.path.exists("C:/Computations/Calc/" + pdbid + "_TM/json"):
            my_path = "C:/Computations/Calc/" + pdbid + "_TM/json/"
            with open(my_path + 'data.json') as f:
                py_json = json.load(f)
            my_list.append(py_json)
    return my_list


# my_pores = get_pores_from_channelsdb("Content.txt")
# download_jsons(my_pores) #works


with open('pores_no_1gmk.txt') as f:  # get the list of pores
        my_pores = f.readlines()
        for i in range(len(my_pores)):
            my_pores[i] = my_pores[i].rstrip('\n')
print(len(my_pores))
with open('pores.txt') as f:  # get the list of pores
    my_pores_all = f.readlines()
    for i in range(len(my_pores_all)):
        my_pores_all[i] = my_pores_all[i].rstrip('\n')
print(len(my_pores_all))

my_list = load_all(get_list_json_ending(my_pores))  # load a list containing all jsons in form of python dicts
print(len(my_list))
print(my_list[1])
print(get_stat_property(my_list, 'charge'))
print(get_stat_property(my_list, 'hydrophobicity'))
print(get_stat_property(my_list, 'hydropathy'))
print(get_stat_property(my_list, 'polarity'))
print(get_stat_property(my_list, 'mutability'))
print(get_stat_property(my_list, 'length'))
print('')
# histogram_property(my_list, 'charge')
# histogram_property(my_list, 'length')

list_fasta = get_list_from_dir('C:/Users/Jirkův NB/Documents/CEITEC/fasta/membrane_proteins/')
my_tm_proteins = get_res_all_tm(list_fasta)
my_tm_proteins.pop('X')
my_basic = load_all_basic(my_pores_all)
basic_stat = get_stat_residues(my_basic, 'mole')
del my_basic
my_tm = load_all_tm(my_pores_all)

print(get_res_fasta('C:/Users/Jirkův NB/Documents/CEITEC/fasta/membrane_proteins/A0A0B4ZYM1.fasta'))
print(get_residues(load_json('1z98.json'), 'general'))
print(get_residues(load_json('1z98.json'), 'bottleneck'))
print(show_residues_ascending(get_percentage(my_tm_proteins)))
print(show_residues_ascending(get_percentage(get_stat_residues(my_list, 'general'))))
print(show_residues_ascending(get_percentage(get_stat_residues(my_list, 'bottleneck'))))
print(show_residues_ascending(get_percentage(basic_stat)))
print(show_residues_ascending(get_percentage(get_stat_residues(my_tm, 'mole'))))

print('')

"""
with open('residues.txt', 'w') as f:
    for item in show_residues_ascending(get_percentage(my_tm_proteins)):
        f.write(str(item))
    f.write('\n')
    for item in show_residues_ascending(get_percentage(get_stat_residues(my_list, 'general'))):
            f.write(str(item))
    f.write('\n')
    for item in show_residues_ascending(get_percentage(get_stat_residues(my_list, 'bottleneck'))):
        f.write(str(item))
    f.write('\n')
    for item in show_residues_ascending(get_percentage(basic_stat)):
        f.write(str(item))
    f.write('\n')
    for item in show_residues_ascending(get_percentage(get_stat_residues(my_tm, 'mole'))):
        f.write(str(item))
"""

"""
no_tm_pores = []
for item in my_pores:
    if item not in membrane_p:
        no_tm_pores.append(item)
print(no_tm_pores)
print(len(no_tm_pores))

"""

# download_cif('6g8z')

# download_all_cif(my_pores_all)

# change_template('template.json', '3s3w')
# generate_mole_jsons(False, my_pores_all)
# generate_mole_jsons(True, my_pores_all)

"""
my_basic = load_all_basic(my_pores_all)
print(len(my_basic))
print(get_residues_mole(my_basic[0]))
basic_stat = get_stat_residues_mole(my_basic)
# hist_random_mole(my_basic)
del my_basic
my_tm = load_all_tm(my_pores_all)
print(len(my_tm))
print(get_residues_mole(my_tm[0]))
"""

print(type(get_all_memprotmd_references('TCDB')))
# get_database_classes(get_all_memprotmd_references('TCDB'), 'TCDB')
print(len(get_all_memprotmd_references('TCDB')))

"""
# mem_uniprot = get_all_uniprotid(membrane_p)
with open('mem_prot_uniprot.txt', 'r') as f:
    mem_prot_uniprot = f.readlines()
    for i in range(len(mem_prot_uniprot)):
        mem_prot_uniprot[i] = mem_prot_uniprot[i].rstrip('\n')
# download_fasta(mem_prot_uniprot, 'C:/Users/Jirkův NB/Documents/CEITEC/fasta/membrane_proteins/')
"""

tcdb_class = get_dict_classes('TCDB')
list_class_names = []
for key in tcdb_class:
    list_class_names.append(key)
print(list_class_names)
print(len(list_class_names))

"""
with open('classes_tcdb.txt', 'w') as f:
    with open('classes_tcdb_filled.txt', 'w') as fi:
        for key in tcdb_class:
            i = 0
            for pdbid in my_pores_all:
                if pdbid in tcdb_class[key]:
                    i += 1
            # print(key + ': ' + str(i) + '/' + str(len(tcdb_class[key])))
            try:
                f.write(key + ': ' + str(i) + '/' + str(len(tcdb_class[key])) + '\n')
                if i > 0:
                    fi.write(key + ': ' + str(i) + '/' + str(len(tcdb_class[key])) + '\n')
                if i > 10:
                    print(key + ': ' + str(i) + '/' + str(len(tcdb_class[key])))
            except UnicodeEncodeError:
                pass
"""
"""
accesion_tcbd = get_dict_classes('TCDB')
print(accesion_tcbd)
with open('accession_tcdb_filled.txt', 'w') as f:
    for key in accesion_tcbd:
        i = 0
        for pdbid in my_pores_all:
            if pdbid in accesion_tcbd[key]:
                i += 1
        # print(key + ': ' + str(i) + '/' + str(len(tcdb_class[key])))
        if i > 0:
            f.write(key + ': ' + str(i) + '/' + str(len(accesion_tcbd[key])) + '\n')
pores_tcbd = get_pores_tcdb(get_dict_classes('TCDB'))
with open('pores_tcdb_list.txt', 'w') as f:
    for key in pores_tcbd:
        f.write("%s\n" % key)  # 'accession/title'
        for item in pores_tcbd[key]:
            f.write("%s\n" % item)
"""

