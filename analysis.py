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


def get_charge(d):
    """"
    returns the charge of the pore
    parameter - dictionary from json file
    """
    charge = d['channels']['transmembranePores'][0]['properties']['charge']
    return charge


def get_stat_charge(l):
    """"
    returns the charge statistics (mean, stdev, min, max) of all pores in the list
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
    returns the hydrophobicity statistics (mean, stdev, min, max) of all pores in the list
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
    returns the hydropathy statistics (mean, stdev, min, max) of all pores in the list
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
    returns the polarity statistics (mean, stdev, min, max) of all pores in the list
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
    returns the mutability statistics (mean, stdev, min, max) of all pores in the list
    parameter - list of dictionaries from json files
    """
    mutabilities = []
    for d in l:
        mutabilities.append(get_mutability(d))
    return statistics.mean(mutabilities), statistics.stdev(mutabilities), min(mutabilities), max(mutabilities)


def hist_charge(l):
    """"
    returns the histogram of charges of all pores in the list
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
    returns the histogram of hydropathies of all pores in the list
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
    returns the histogram of hydrophobicities of all pores in the list
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
    returns the histogram of polarities of all pores in the list
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
    returns the histogram of mutabilities of all pores in the list
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
    """"
     returns the dictionary containing the residues as keys and its quantity as values
     parameter - dictionary from json file
    """
    residues = {}
    residues_json = d['channels']['transmembranePores'][0]['layers']['residueFlow']  # list

    for aa in residues_json:
        aa = aa[:3]
        if aa not in residues:
            residues[aa] = 1
        else:
            residues[aa] += 1
    for aa in residues:
        if aa in compare_residues(residues_json):
            residues[aa] -= compare_residues(residues_json)[aa]
    return residues


def sort_residues(d):
    """"
    returns the dictionary containing the sorted residues as keys and its quantity as values
    sorts the residues by eliminating the het molecules
    parameter - dictionary containing the unsorted residues
    """
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
    """"
    returns the dictionary containing the residue quantity of all pores in the list
    parameter - list of dictionaries from json files
    """
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
    """"
    returns the residue statistics (total_residues, mean_by_residue, mean_by_pore) of all pores in the list
    parameter - list of dictionaries from json files
    """
    d = get_stat_residues(l)
    res_nb = 0
    for aa in d:
        res_nb += d[aa]
    return res_nb, res_nb/len(d), res_nb/len(l)


def show_residues_ascending(d):
    """"
    sorts a dictionary by ascending value
    parameter - dictionary with residues
    """
    return sorted(d.items(), key=operator.itemgetter(1))


def get_residues_from_bottleneck(d):
    """"
    returns the dictionary containing the residue quantity in the bottleneck of this json file
    parameter - dictionary from json files
    """
    residues = {}
    residues_json = []
    for i in range(len(d['channels']['transmembranePores'][0]['layers']['layersInfo'])):
        if d['channels']['transmembranePores'][0]['layers']['layersInfo'][i]['layerGeometry']['bottleneck']:
            residues_json = d['channels']['transmembranePores'][0]['layers']['layersInfo'][i]['residues']
    for aa in residues_json:
        aa = aa[:3]
        if aa not in residues:
            residues[aa] = 1
        else:
            residues[aa] += 1
    for aa in residues:
        if aa in compare_residues(residues_json):
            residues[aa] -= compare_residues(residues_json)[aa]
    return residues


def get_stat_bottleneck(l):
    """"
    returns the dictionary containing the residue quantity in the bottleneck of all pores in the list
    parameter - list of dictionaries from json files
    """
    btn_residues = {}
    residues = []
    for d in l:
        residues.append(get_residues_from_bottleneck(d))
    for d in residues:
        for aa in d:
            aa = aa[:3]
            if aa not in btn_residues:
                btn_residues[aa] = d[aa]
            else:
                btn_residues[aa] += d[aa]
        btn_residues = sort_residues(btn_residues)
    return btn_residues


def get_stat_btn_number(l):
    """"
    returns the residue statistics (total_residues, mean_by_residue, mean_by_pore) of all pores in the list
    parameter - list of dictionaries from json files
    """
    d = get_stat_bottleneck(l)
    btn_nb = 0
    for aa in d:
        btn_nb += d[aa]
    return btn_nb, btn_nb/len(d), btn_nb/len(l)


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
        ave_dict[key] = round(d[key]/(nb/len(d)), 2)
    return ave_dict


def get_length(d):
    """"
    returns the length of the pore from json file
    parameter - dictionary from json files
    """
    return d['channels']['transmembranePores'][0]['layers']['layersInfo'][-1]['layerGeometry']['endDistance']


def get_stat_length(l):
    """"
    returns the length statistics (mean, stdev, min, max) of all pores in the list
    parameter - list of dictionaries from json files
    """
    lengths = []
    for d in l:
        lengths.append(get_length(d))
    return statistics.mean(lengths), statistics.stdev(lengths), min(lengths), max(lengths)


def hist_length(l):
    """"
    returns the histogram of lengths of all pores in the list
    parameter - list of dictionaries from json files
    """
    lengths = []
    for d in l:
        lengths.append(get_length(d))
    num_bins = 100
    n, bins, patches = plt.hist(lengths, num_bins, facecolor='black', alpha=0.5)
    plt.xlabel('Length')
    plt.ylabel('Quantity')
    plt.title('Histogram of Length')
    return plt.show()


def compare_residues(l):
    """
    finds residues aligning two times in the ResidueFlow (backbone)
    :param l: list of residues
    :return: residues present two times
    """
    backbones = []
    backbones_d = {}
    for aa in l:
        if aa + ' Backbone' in l:
            backbones.append(aa)
    for aa in backbones:
        aa = aa[:3]
        if aa not in backbones_d:
            backbones_d[aa] = 1
        else:
            backbones_d[aa] += 1
    return backbones_d


def get_all_memprotmd_references():
    """
    commands all the references from the mpm database from the MemProtMD site
    :return: json file
    """
    MEMPROTMD_ROOT_URI = "http://memprotmd.bioch.ox.ac.uk/"
    return requests.post(MEMPROTMD_ROOT_URI + "api/references/all/mpm").json()


def get_mpm_classes(l):
    """
    creates a text file with all the references from mpm
    :param l: list of references
    :return: None
    """
    with open('mpm.txt', 'w') as f:
        for d in l:
            f.write("%s\n" % d['accession'])
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


def get_percentage(d):
    new_d = {}
    ave_d = average_d(d)
    for aa in ave_d:
        new_d[aa] = 5*ave_d[aa]
    return new_d


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

# hist_charge(my_list)
# hist_hydropathy(my_list)
# hist_hydrophobicity(my_list)
# hist_polarity(my_list)
# hist_mutability(my_list)
print(get_residues_from_bottleneck(my_list[0]))
btnres = get_stat_bottleneck(my_list)
print(btnres)
print('')


"""
with open('residues.txt', 'w') as f:
    for item in show_residues_ascending(average_d(resid)):
        f.write(str(item))
    f.write('\n')
    for item in show_residues_ascending(average_d(btnres)):
            f.write(str(item))
"""

print(get_length(my_list[1]))
print(get_stat_length(my_list))
# hist_length(my_list)

# ResidueFlow correction
print(get_residues(load_json('3jaf.json')))

get_mpm_classes(get_all_memprotmd_references())
print(len(get_all_memprotmd_references()))

membrane_p = []
for item in get_all_memprotmd_references()[3]['simulations']:
    membrane_p.append(item[:4])
print(len(membrane_p))

"""
# mem_uniprot = get_all_uniprotid(membrane_p)
with open('mem_prot_uniprot.txt', 'r') as f:
    mem_prot_uniprot = f.readlines()
    for i in range(len(mem_prot_uniprot)):
        mem_prot_uniprot[i] = mem_prot_uniprot[i].rstrip('\n')
# download_fasta(mem_prot_uniprot, 'C:/Users/Jirkův NB/Documents/CEITEC/fasta/membrane_proteins/')
"""

print(len(my_pores))
l_classes = get_list_from_dir('C:/Users/Jirkův NB/Documents/CEITEC/fasta/')
print(l_classes)
l_classes.remove('membrane_proteins')
for f in l_classes:
    print(f + ": " + str(count_presence(f, my_pores)))

no_tm_pores = []
for item in my_pores:
    if item not in membrane_p:
        no_tm_pores.append(item)
print(no_tm_pores)
print(len(no_tm_pores))

print(get_res_fasta('C:/Users/Jirkův NB/Documents/CEITEC/fasta/membrane_proteins/A0A0B4ZYM1.fasta'))
list_fasta = get_list_from_dir('C:/Users/Jirkův NB/Documents/CEITEC/fasta/membrane_proteins/')
my_d_all_tm = get_res_all_tm(list_fasta)
print(show_residues_ascending(average_d(my_d_all_tm)))
print(len(my_d_all_tm))
print(get_stat_res_number_all_tm(list_fasta))
my_d_all_tm.pop('X')
print(show_residues_ascending(get_percentage(my_d_all_tm)))
print(show_residues_ascending(get_percentage(resid)))
print(show_residues_ascending(get_percentage(btnres)))

# download_cif('6g8z')
with open('pores.txt') as f:  # get the list of pores
    my_pores_all = f.readlines()
    for i in range(len(my_pores_all)):
        my_pores_all[i] = my_pores_all[i].rstrip('\n')
print(len(my_pores_all))
# download_all_cif(my_pores_all)

change_template('template.json', '3s3w')
generate_mole_jsons(False, my_pores_all)
generate_mole_jsons(True, my_pores_all)
