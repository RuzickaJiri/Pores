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


def download_uniprot_xml(pdbid):
    """"
    downloads the Uniprot xml file from pdb id
    the pdbid as the parameter
    """
    urllib.request.urlretrieve("https://www.uniprot.org/uniprot/?query=" + pdbid + "&sort=score&format=xml",
                               pdbid + ".xml")


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
        d.update({id:chains})
    return d


def find_all(file):
    """"
    returns a tuple containing uniprot id, tm regions and pdb
    """
    tree = load_xml(file)
    return find_uniprot_id(tree),find_tm_regions(tree), find_pdb(tree)


def load_all(l):
    """"
    returns a list of tuple from all pores in the list specified as a parameter
    """
    list_data = []
    for i in range(len(l)):
        list_data.append(find_all(l[i] + '.xml'))
        # print(find_all(l[i] + '.xml'))
    return list_data


with open('pores.txt') as f:  # get the list of pores
    my_pores = f.readlines()
    for i in range(len(my_pores)):
        my_pores[i] = my_pores[i].rstrip('\n')
print(my_pores)
# download_xml(my_pores) # works

my_tree = load_xml('3s3w.xml')
print(find_uniprot_id(my_tree))
print(find_tm_regions(my_tree))
print(find_pdb(my_tree))
print(find_all('3s3w.xml'))
print(str(find_all(my_pores[0] + '.xml')) + '\n')
# print(load_all(my_pores)) #  works
list_pores = load_all(my_pores)
for i in range(len(list_pores)):
    print(str(list_pores[i]) + '/n')
print(len(list_pores))
list_tm = []
for i in range(len(list_pores)):
    if list_pores[i][1] != [] and list_pores[i][2] != {}:
        list_tm.append(list_pores[i])
print(len(list_tm))

k = 0
PORES_UP = []
for i in range(len(list_tm)):
    PORES_UP.append(my_pores[i].upper())
    if PORES_UP[i] in list_tm[3]:
        k += 1
print(k)
print(PORES_UP)
