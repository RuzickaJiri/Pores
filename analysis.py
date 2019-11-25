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
import script as sc

# Downloading and loading data


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


def load_json(pdbid, mem):
    """"
    returns a dictionary issued from the json file specified as a parameter
    """

    # if os.path.exists(pdbid + ".json"):
      #  with open(pdbid + ".json") as f:
       #     py_json = json.load(f)
        #return py_json
    if mem:
        with open("C:\Computations\Calc\FromMoleMem\\" + pdbid + "\\json\data.json") as f:
            py_json = json.load(f)
        return py_json
    else:
        with open("C:\Computations\Calc\FromMole\\" + pdbid + "\\json\data.json") as f:
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
    loads all jsons from a list and returns a list of dictionaries (previously json files)
    l is a list that contains the names of json files
    """
    list_json = []
    for i in range(len(l)):
        list_json.append(load_json(l[i]))
    return list_json


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
    """
    changes the template of json file destined to be sent to MOLE
    creates a new file with the given data
    :param tm: binary - TM-only parameter
    :param pdbid: pdb id of the structure
    """
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
    """
    generates json files for all the structures in the list
    :param template: template of json file
    :param l: list of pdbids
    """
    for pdbid in l:
        change_template(template, pdbid)


def load_all_mole(path, tm, l):
    """
    loads all json files from the list in a new list
    :param path: string, path to the files to upload
    :param tm: binary - TM-only parameter
    :param l: list of pdbids
    :return: list of json data
    """
    my_list = []
    if tm:
        for pdbid in l:
            if os.path.exists(path + pdbid + "_TM/json"):
                my_path = path + pdbid + "_TM/json/"
                with open(my_path + 'data.json') as f:
                    py_json = json.load(f)
                my_list.append(py_json)
    else:
        for pdbid in l:
            if os.path.exists(path + pdbid + "\\json"):
                my_path = path + pdbid + "/json/"
                with open(my_path + 'data.json') as f:
                    py_json = json.load(f)
                my_list.append(py_json)
    return my_list

# Properties analysis


def analyze_property(l, prop, mem):
    """
    rewritten method for faster loading of the given property from a json file
    :param l: list of pdb ids
    :param prop: string, property to analyze
    :return: array of all property values from the list
    """
    properties = []
    for pdbid in l:
        temp = get_property(load_json(pdbid, mem), prop)
        if temp is not None:
            properties.append(temp)
    return properties


def get_property(d, prop):
    """"
    returns the specified property of the pore
    parameter - dictionary from json file
    """
    try:
        if prop == 'length':
            if 'channels' in d:
                try:
                    return d['channels']['transmembranePores'][0]['layers']['layersInfo'][-1]['layerGeometry']['endDistance']
                except IndexError:
                    return d['channels']['reviewedChannels'][0]['layers']['layersInfo'][-1]['layerGeometry']['endDistance']
            else:
                return d['Channels']['Paths'][0]['Layers']['LayersInfo'][-1]['LayerGeometry']['EndDistance']
        elif prop == 'bottleneck':
            if 'channels' in d:
                try:
                    for i in range(len(d['channels']['transmembranePores'][0]['layers']['layersInfo'])):
                        if d['channels']['transmembranePores'][0]['layers']['layersInfo'][i]['layerGeometry']['bottleneck']:
                            return d['channels']['transmembranePores'][0]['layers']['layersInfo'][i]['minRadius']
                except IndexError:
                    for i in range(len(d['channels']['transmembranePores'][0]['layers']['layersInfo'])):
                        if d['channels']['reviewedChannels'][0]['layers']['layersInfo'][i]['layerGeometry']['bottleneck']:
                            return d['channels']['reviewedChannels'][0]['layers']['layersInfo'][i]['minRadius']
            else:
                for i in range(len(d['Channels']['Paths'][0]['Layers']['LayersInfo'])):
                    if 'Bottleneck' in d['Channels']['Paths'][0]['Layers']['LayersInfo'][i]['LayerGeometry']:
                        return d['Channels']['Paths'][0]['Layers']['LayersInfo'][i]['LayerGeometry']['MinRadius']
        elif prop == 'residues#':
            if 'channels' in d:
                return len(d['channels']['transmembranePores'][0]['layers']['residueFlow'])  # list
            else:
                return len(d['Channels']['Paths'][0]['Layers']['ResidueFlow'])
        else:
            if 'channels' in d:
                try:
                    return d['channels']['transmembranePores'][0]['properties'][prop]
                except IndexError:
                    return d['channels']['reviewedChannels'][0]['properties'][prop]
            else:
                propy = prop[0].upper() + prop[1:]
                return d['Channels']['Paths'][0]['Properties'][propy]
    except IndexError:
        try:
            print(d['Config']['PdbId'])
        except KeyError:
            print(d)# ['Config']['PdbId'])


def get_stat_property(l, prop, mem):
    """"
    returns the given property statistics (mean, stdev, min, max) of all pores in the list
    parameter - list of dictionaries from json files, string property
    """
    properties = analyze_property(l, prop, mem)
    return round(statistics.mean(properties),2), round(statistics.stdev(properties),2), min(properties), max(properties)


def histogram_property(l, prop, mem):
    """"
    returns the histogram of the given property of all pores in the list
    parameter - list of dictionaries from json files, the string property
    """
    num_bins = 75
    arr = analyze_property(l, prop, mem)
    n, bins, patches = plt.hist(arr, num_bins, facecolor='black', alpha=0.5)
    plt.xlabel(prop)
    plt.ylabel('Quantity')
    plt.title('Histogram of ' + prop)
    return plt.show()


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

# Residue analysis


def analyze_residues(l, place):
    """
    rewritten method for faster loading of residues from a json file
    :param l: list of pdb ids
    :param place: string, place specification to analyze (general, bottleneck)
    :return: array of all residues from the json files of the given list
    """
    all_residues = {}
    residues = []
    for pdbid in l:
        residues.append(get_residues(load_json(pdbid), place))
    for d in residues:
        for aa in d:
            aa = aa[:3]
            if aa not in all_residues:
                all_residues[aa] = d[aa]
            else:
                all_residues[aa] += d[aa]
    all_residues = sort_residues(all_residues)
    return all_residues


def get_residues(d, place):
    """"
     returns the dictionary containing the residues as keys and its quantity as values
     parameter - dictionary from json file
    """
    residues = {}
    residues_json = []
    if place == 'general':
        if 'channels' in d:
            residues_json = d['channels']['transmembranePores'][0]['layers']['residueFlow']  # list
        else:
            residues_json = d['Channels']['Paths'][0]['Layers']['ResidueFlow']
    elif place == 'bottleneck':
        if 'channels' in d:
            for i in range(len(d['channels']['transmembranePores'][0]['layers']['layersInfo'])):
                if d['channels']['transmembranePores'][0]['layers']['layersInfo'][i]['layerGeometry']['bottleneck']:
                    residues_json = d['channels']['transmembranePores'][0]['layers']['layersInfo'][i]['residues']
        else:
            for i in range(len(d['Channels']['Paths'][0]['Layers']['LayersInfo'])):
                if 'Bottleneck' in d['Channels']['Paths'][0]['Layers']['LayersInfo'][i]['LayerGeometry']:
                    residues_json = d['Channels']['Paths'][0]['Layers']['LayersInfo'][i]['Residues']

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
    # d = get_stat_residues(l, place)
    d = analyze_residues(l, place)
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

# Additional


def count_presence_lists(ll, lr):
    """
    computes the number of intersections between two lists
    :param ll: left list
    :param lr: right list
    :return: number of common items
    """
    i = 0
    for item in ll:
        if item in lr:
            i += 1
    return i


def count_presences_lists(ll, lr, lt):
    """
    computes the number of intersections between three lists
    :param ll: left list
    :param lr: right list
    :param lt: third list
    :return: number of common items
    """
    i = 0
    for item in ll:
        if item in lr:
            if item in lt:
                i += 1
    return i


def intersection_list(ll, lr):
    """
    computes the intersection of two lists
    :param ll: left list
    :param lr: right list
    :return: new list of intersection
    """
    new_l = []
    for item in ll:
        if item in lr:
            new_l.append(item)
    return new_l


def no_intersection_list(ll, lr):
    """
    finds if the items in left list are also in the right list, returns the unique items
    :param ll: left list
    :param lr: right list
    :return: list with items only in the left list
    """
    new_l = []
    for item in ll:
        if item not in lr:
            new_l.append(item)
    return new_l


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


def list_from_dict(d):
    """
    creates a list from a dictionary keys
    :param d: dictionary
    :return: list
    """
    new_l = []
    for key in d:
        for pdbid in d[key]:
            if pdbid not in new_l:
                new_l.append(pdbid)
    return new_l


def memprotmd_text_search(search_term, in_databases=["mpm", "TCDB", "mpstruc"], size_bias=0.1, n_results=10):
    """
    text research method from MemProt API GitHub
    :param search_term: search term
    :param in_databases: database
    :param size_bias: size bias
    :param n_results: number of results
    :return: research result
    """
    return requests.post(

        MEMPROTMD_ROOT_URI
        + "api/search/simple",

        json={
            "searchTerm": search_term,
            "inDatabases": in_databases,
            "numResults": n_results,
            "sizeBias": size_bias
        }

    ).json()


def memprotmd_advanced_search(collection_name, query, projection, options):
    """
    advanced research method from MemProt API GitHub
    :param collection_name: name of the collection
    :param query: query
    :param projection: projection
    :param options: options
    :return: advanced research
    """
    return requests.post(

        MEMPROTMD_ROOT_URI
        + "api/search/advanced",

        json={
            "collectionName": collection_name,
            "query": query,
            "projection": projection,
            "options": options
        }

    ).json()


def get_database_classes(l, database):
    """
    creates a text file with all the references from mpm
    :param l: list of references
    :param database: string, the given database
    :return: None
    """
    with open(database + '.txt', 'w') as f:
        for d in l:
            f.write("%s\n" % d['accession'])  # 'accession/title'
            for item in d['simulations']:
                f.write("%s\n" % item[:4])


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
    """
    computes the percentage of every residue in the dictinary, rounds the percentage with two decimals
    :param d: dictionary of residues
    :return: dictionary with the percentage of residues
    """
    new_d = {}
    ave_d = average_d(d)
    for aa in ave_d:
        new_d[aa] = round((100/len(ave_d))*ave_d[aa], 2)
    return new_d


if __name__ == "__main__":

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

    # my_list = load_all(my_pores)  # load a list containing all jsons in form of python dicts
    # print(len(my_list))
    # print(my_list[1])
    print(get_stat_property(my_pores, 'charge'))
    print(get_stat_property(my_pores, 'hydrophobicity'))
    print(get_stat_property(my_pores, 'hydropathy'))
    print(get_stat_property(my_pores, 'polarity'))
    print(get_stat_property(my_pores, 'mutability'))
    print(get_stat_property(my_pores, 'length'))
    print('')
    # histogram_property(my_pores, 'charge')
    # histogram_property(my_pores, 'length')

    # Residues analysis
    """
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
    """

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

    print(type(sc.get_all_memprotmd_references('TCDB')))
    # get_database_classes(get_all_memprotmd_references('TCDB'), 'TCDB')
    print(len(sc.get_all_memprotmd_references('TCDB')))

    """
    # mem_uniprot = get_all_uniprotid(membrane_p)
    with open('mem_prot_uniprot.txt', 'r') as f:
        mem_prot_uniprot = f.readlines()
        for i in range(len(mem_prot_uniprot)):
            mem_prot_uniprot[i] = mem_prot_uniprot[i].rstrip('\n')
    # download_fasta(mem_prot_uniprot, 'C:/Users/Jirkův NB/Documents/CEITEC/fasta/membrane_proteins/')
    """

    mpm_class = sc.get_dict_classes('mpm')
    list_class_names = []
    for key in mpm_class:
        list_class_names.append(key)
    print(list_class_names)
    print(len(list_class_names))

    """
    with open('classes_mpm.txt', 'w') as f:
        with open('classes_mpm_filled.txt', 'w') as fi:
            for key in mpm_class:
                i = 0
                for pdbid in my_pores_all:
                    if pdbid in mpm_class[key]:
                        i += 1
                # print(key + ': ' + str(i) + '/' + str(len(tcdb_class[key])))
                try:
                    f.write(key + ': ' + str(i) + '/' + str(len(mpm_class[key])) + '\n')
                    if i > 0:
                        fi.write(key + ': ' + str(i) + '/' + str(len(mpm_class[key])) + '\n')
                    if i > 10:
                        print(key + ': ' + str(i) + '/' + str(len(mpm_class[key])))
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

    MEMPROTMD_ROOT_URI = "http://memprotmd.bioch.ox.ac.uk/"
    print(memprotmd_text_search("channel", size_bias=10))

    """
    print(memprotmd_advanced_search(
        # Use the simulation collection
        "refs",
        # Look in the chains array of each simulation and see if
        # any element in the array has the field tm_alpha equal to 7
        {
            "accession": {
                "$in": ["outer-membrane-carboxylate-channels-occ", "beta-barrel-membrane-proteins-porins-and-relatives"]
            }
        },
        # Projection - choose the fields to return. ID is returned
        # by default.
        {
        },
        # Options - sort and then limit
        {
        }
    ))
    
    print(memprotmd_advanced_search("refs",
        {"accession": {"$in": ["channels-mechanosensitive", "channels-potassium-sodium-proton-ion-selective",
                               "channels-calcium-ion-selective", "channels-transient-receptor-potential-trp",
                               "channels-other-ion-channels", "channels-fluc-family", "channels-urea-transporters",
                               "channels-aquaporins-and-glyceroporins", "channels-formate-nitrite-transporter-fnt-family",
                               "channels-gap-junctions", "channels-amt-mep-rh-proteins"]}},
                                    {}, {}))
    """

    pores_tcdb_list = list_from_dict(sc.get_pores_from_db(sc.get_dict_classes('TCDB'), 'TCDB'))
    alpha_mpstruc_list = list_from_dict(sc.get_pores_from_db(sc.get_dict_classes('mpstruc'), 'mpstruc-alpha'))
    beta_mpstruc_list = list_from_dict(sc.get_pores_from_db(sc.get_dict_classes('mpstruc'), 'mpstruc-beta'))
    pores_mpm_list = list_from_dict(sc.get_pores_from_db(sc.get_dict_classes('mpm'), 'mpm'))

    mpstruc_list = alpha_mpstruc_list + beta_mpstruc_list
    all_pdbid_pores = []
    list_db = []
    list_db.append(pores_tcdb_list), list_db.append(mpstruc_list), list_db.append(pores_mpm_list)
    for db in list_db:
        for pdbid in db:
            if pdbid not in all_pdbid_pores:
                all_pdbid_pores.append(pdbid)

    print(len(pores_mpm_list))
    print(len(pores_tcdb_list))
    print(len(mpstruc_list))
    print('')
    print(count_presence_lists(my_pores_all, pores_mpm_list))
    print(count_presence_lists(my_pores_all, pores_tcdb_list))
    print(count_presence_lists(my_pores_all, mpstruc_list))
    print(count_presence_lists(my_pores_all, all_pdbid_pores))
    print(count_presence_lists(mpstruc_list, pores_tcdb_list))
    print(count_presence_lists(mpstruc_list, pores_mpm_list))
    print(count_presence_lists(pores_mpm_list, pores_tcdb_list))
    print(count_presences_lists(pores_mpm_list, pores_tcdb_list, mpstruc_list))
    print(len(all_pdbid_pores))
    # download_all_cif(check_db()[0])
    print(analyze_property(my_pores, 'charge'))
