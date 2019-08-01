#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 2019

@author: jruzicka
"""

# Libraries
import os
import urllib.request
import xml.etree.ElementTree as ET
import statistics
import json
import re


def load_json(file):
    """"
    returns a dictionary issued from the json file specified as a parameter
    """
    with open(file) as f:
        py_json = json.load(f)
    return py_json


def find_uniprot_from_sifts(pdbid):
    urllib.request.urlretrieve("http://www.ebi.ac.uk/pdbe/api/mappings/uniprot/" + pdbid,  pdbid + "up.json")
    d = load_json(pdbid + "up.json")
    id = list(d[pdbid]['UniProt'].keys())[0]
    return id


def get_all_uniprotid(l):
    uniprotid = []
    for pdbid in l:
        try:
            uniprotid.append(find_uniprot_from_sifts(pdbid))
        except:
            print(pdbid)
    return uniprotid


def download_uniprot_xml(uniprotid):
    """"
    downloads the Uniprot xml file from uniprot id
    the uniprot id as the parameter
    """
    urllib.request.urlretrieve("https://www.uniprot.org/uniprot/" + uniprotid + ".xml",
                               uniprotid + ".xml")


def download_xml(l):
    """"
    downloads all xml files from a list
    (in this case a list of pores we want to get)
    list contains the pdbids of the pores
    """
    for name in l:
        download_uniprot_xml(name)


def load_xml(file):
    """"
    opens an xml file, returns a tree
    """
    with open(file) as f:
        tree = ET.parse(f)
    return tree


def find_uniprot_id(tree):
    """"
    returns the uniprot id from a tree
    """
    return tree.getroot()[0][0].text


def find_tm_regions(tree):
    """"
    returns a list of transmembrane regions
    """
    l = []
    for elements in tree.getroot()[0].findall("./{http://uniprot.org/uniprot}feature/[@type='transmembrane region']"):
        begin = elements[0][0].attrib['position']
        end = elements[0][1].attrib['position']
        l.append(begin + '-' + end)
    return l


def find_pdb(tree):
    """"
    returns a dictionary in form {pdbid : chains}
    """
    d = {}
    for elements in tree.getroot()[0].findall("./{http://uniprot.org/uniprot}dbReference/[@type='PDB']"):
        id = elements.attrib['id']
        try:
            chains = elements[2].attrib['value']
        except IndexError:
            chains = elements[1].attrib['value']
        d.update({id: chains})
    return d


def find_all(file):
    """"
    returns a tuple containing uniprot id, tm regions and pdb
    """
    tree = load_xml(file)
    return find_uniprot_id(tree), find_tm_regions(tree), find_pdb(tree)


def load_all(l):
    """"
    returns a list of tuple from all pores in the list specified as a parameter
    """
    list_data = []
    for i in range(len(l)):
        list_data.append(find_all(l[i] + '.xml'))
        # print(find_all(l[i] + '.xml'))
    return list_data


def get_first_st(t, pdb):
    return int(re.findall('\d+', t[2][pdb])[0])


def get_last_st(t, pdb):
    return int(re.findall('\d+', t[2][pdb])[-1])


def get_first_tm(t):
    return int(re.findall('\d+', t[1][0])[0])


def get_last_tm(t):
    return int(re.findall('\d+', t[1][-1])[-1])


def compare_tm_structure(t, pdb):
    """"
    compares if the pdb structure satisfies the transmembrane region
    """
    tm1 = get_first_tm(t)
    tm2 = get_last_tm(t)
    st1 = get_first_st(t, pdb)
    st2 = get_last_st(t, pdb)
    if tm1 >= st1 and tm2 <= st2:
        return True
    else:
        return False

def compare_tm_structure_rupt(t, pdb): # not finished
    tm1 = get_first_tm(t)
    tm2 = get_last_tm(t)
    st1 = get_first_st(t, pdb)
    st2 = get_last_st(t, pdb)
    rup1 = int(re.findall('\d+', t[2][pdb])[1])
    rup2 = int(re.findall('\d+', t[2][pdb])[2])
    if tm1 >= st1 and tm2 <= st2:
        for i in list(map(int, re.findall('\d+', t[1][0]))):
            if i > rup1 and i < rup2:
                return False
            else:
                return True
    else:
        return False


def compare_all(list_tmr, list_por):
    result = {}
    pupper = []
    for i in range(len(list_tmr)):
        pupper.append(list_por[i].upper())
        if pupper[i] in list_tmr[i][2] and len(re.findall('\d+', list_tmr[i][2][pupper[i]])) < 3:
            result.update({pupper[i]: compare_tm_structure(list_tmr[i], pupper[i])})
        # elif len(re.findall('\d+', list_tmr[i][2][pupper[i]])) > 2:
            # result.update({pupper[i]: compare_tm_structure_rupt(list_tmr[i], pupper[i])})

    return result


with open('pores.txt') as f:  # get the list of pores
    my_pores = f.readlines()
    for i in range(len(my_pores)):
        my_pores[i] = my_pores[i].rstrip('\n')
print(len(my_pores))
# download_xml(my_pores) # works
# my_uni  = get_all_uniprotid(my_pores)

with open('uniprotid.txt') as f:  # get the list of proteins
    my_uni = f.readlines()
    for i in range(len(my_uni)):
        my_uni[i] = my_uni[i].rstrip('\n')
print(my_uni)
# download_xml(my_uni) #  works

# my_tree = load_xml('3s3w.xml')
# print(find_uniprot_id(my_tree))
# print(find_tm_regions(my_tree))
# print(find_pdb(my_tree))
# print(find_all('3s3w.xml'))
print(str(find_all(my_uni[0] + '.xml')) + '\n')
# print(load_all(my_pores)) #  works
list_pores = load_all(my_uni)

with open('no_uniprot.txt') as f:  # get the list of proteins
    my_no_uni = f.readlines()
    for i in range(len(my_no_uni)):
        my_no_uni[i] = my_no_uni[i].rstrip('\n')
print(my_no_uni)
for id in my_no_uni:
    if id in my_pores:
        my_pores.remove(id)
print(len(my_pores))

print(len(list_pores))
list_tm = []  # list of pores with tm and pdb id
to_delete = []
for i in range(len(list_pores)):
    if list_pores[i][1] != [] and list_pores[i][2] != {}:
        list_tm.append(list_pores[i])
    else:
        to_delete.append(my_pores[i])

for id in to_delete:
    my_pores.remove(id)
"""
for i in range(len(list_tm)):
    print(str(list_tm[i]) + '/n')
print(len(list_tm))

PORES_UP = []
for i in range(len(list_tm)):
    PORES_UP.append(my_pores[i].upper())
    if PORES_UP[i] in list_tm[i][2]:
        print(re.findall('\d+', list_tm[i][2][PORES_UP[i]]))
"""
PORES_UP = []
for i in range(len(list_tm)):
    PORES_UP.append(my_pores[i].upper())
print(PORES_UP)
print(get_first_st(list_tm[0], '3S3W'))
print(get_last_st(list_tm[0], '3S3W'))
print(get_first_tm(list_tm[0]))
print(get_last_tm(list_tm[0]))
print(compare_tm_structure(list_tm[0], '3S3W'))
my_di_str = compare_all(list_tm, PORES_UP)
print(my_di_str)
print(len(my_di_str))
c = 0
for key in my_di_str:
    if my_di_str[key]:
        c += 1
print(c)

pdb_keys = []
for key in my_di_str:
    for k in list_pores:
        if key in k[2]:
            for pdbid in k[2]:
                if pdbid not in pdb_keys:
                    pdb_keys.append(pdbid)
print(pdb_keys)
print(len(pdb_keys))
print(len(pdb_keys)-c)

