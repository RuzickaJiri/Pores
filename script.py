#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 2019

@author: jruzicka
"""

# Libraries
import os
import urllib.request
import requests
import xml.etree.ElementTree as ET
import statistics
import json
import re
import analysis as an
import selection as se

# Downloading and loading data


def load_json(file):
    """"
    returns a dictionary issued from the json file specified as a parameter
    """
    with open(file) as f:
        py_json = json.load(f)
    return py_json


def get_all_uniprotid(l):
    """
    gets all uniprot ids from a list
    :param l: list of pdb ids
    :return: list of uniprot ids
    """
    uniprotid = {}
    i = 0
    for pdbid in l:
        try:
            uniprotid[pdbid] = find_uniprot_from_sifts(pdbid)
        except urllib.error.HTTPError:
            i += 1
            # print('No UniProt: ' + pdbid)
            # pass
    print(str(i) + " PDBids have no UniProt mapping.")
    return uniprotid


def download_uniprot_xml(uniprotid):
    """"
    downloads the Uniprot xml file from uniprot id
    the uniprot id as the parameter
    """
    if not os.path.isfile(uniprotid + ".xml"):
        urllib.request.urlretrieve("https://www.uniprot.org/uniprot/" + uniprotid + ".xml",
                               uniprotid + ".xml")


def download_xml(l):
    """"
    downloads all xml files from a list
    (in this case a list of pores we want to get)
    list contains the uniprot of the pores
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


def find_best_uniprot(l, d, pdbid):
    """
    finds the uniprot containing the longest sequence of the structure and returns its position from the list
    :param l: list
    :param d: dictionary
    :param pdbid: pdb id
    :return: integer
    """
    p = 0
    if len(l) > 1:
        max_length = 0
        for i in range(len(l)):
            length = d[pdbid]['UniProt'][l[i]]['mappings'][0]['unp_end'] - d[pdbid]['UniProt'][l[i]]['mappings'][0]['unp_start']
            j = 1
            while d[pdbid]['UniProt'][l[i]]['mappings'][j]['chain_id'] == d[pdbid]['UniProt'][l[i]]['mappings'][0]['chain_id']:
                length += d[pdbid]['UniProt'][l[i]]['mappings'][j]['unp_end'] - d[pdbid]['UniProt'][l[i]]['mappings'][j]['unp_start']
                j += 1
            if length > max_length:
                max_length = length
                p = i
    return p


def find_uniprot_from_sifts(pdbid):
    """
    downloads a json file from sitfs api, finds the uniprot ids and returns the best one. deletes the json file
    :param pdbid: pdbid
    :return: uniprot id
    """
    if not os.path.isfile(pdbid + 'up.json'):
        urllib.request.urlretrieve("http://www.ebi.ac.uk/pdbe/api/mappings/uniprot/" + pdbid,  pdbid + "up.json")
    d = load_json(pdbid + "up.json")
    id = []
    for key in d[pdbid]['UniProt']:
        id.append(key)
    # i = find_best_uniprot(id, d, pdbid)
    # os.remove(pdbid + "up.json")
    return id

# UniProt (TM) selection


def get_first_st(t, pdb):
    """
    finds the first position of the structure in the protein
    :param t: tuple
    :param pdb: pdb id
    :return: integer, first position of the structure
    """
    return int(re.findall('\d+', t[2][pdb])[0])


def get_last_st(t, pdb):
    """
    finds the last position of the structure in the protein
    :param t: tuple
    :param pdb: pdbid
    :return: integer, last position of the structure
    """
    return int(re.findall('\d+', t[2][pdb])[-1])


def get_first_tm(t):
    """
    finds the first position of the TM part in the given protein
    :param t: tuple
    :return: integer, first position of the TM part
    """
    return int(re.findall('\d+', t[1][0])[0])


def get_last_tm(t):
    """
    finds the last position of the TM part in the given protein
    :param t: tuple
    :return: integer, last position of the TM part
    """
    return int(re.findall('\d+', t[1][-1])[-1])


def compare_tm_structure(t, pdb):
    """
    compares if the pdb structure satisfies the transmembrane region
    :param t: tuple, data from xml file
    :param pdb: pdb id
    :return: binary
    """
    tm1 = get_first_tm(t)
    tm2 = get_last_tm(t)
    st1 = get_first_st(t, pdb)
    st2 = get_last_st(t, pdb)
    if tm1 >= st1 and tm2 <= st2:
        return True
    else:
        return False


def compare_all(list_tmr, list_por):
    """
    finds if a structure contains the complete TM part for all the structures in the list
    :param list_tmr: list of tuples from xml files
    :param list_por: list of pdbids
    :return: dictionary pdbid:Binary, says if the structure has all the TM region
    """
    result = {}
    pupper = list_upper(list_por)
    for i in range(len(pupper)):
        for j in range(len(list_tmr)):
            if pupper[i] in list_tmr[j][2] and list_tmr[j][1] != []:  # and len(re.findall('\d+', list_tmr[j][2][pupper[i]])) < 3:
                if pupper[i] not in result:
                    result[pupper[i]] = compare_tm_structure(list_tmr[j], pupper[i])
                elif not result[pupper[i]] or result[pupper[i]] == 'No TM':
                    result[pupper[i]] = compare_tm_structure(list_tmr[j], pupper[i])
            elif pupper[i] in list_tmr[j][2] and list_tmr[j][1] == []:
                if pupper[i] not in result:
                    result[pupper[i]] = 'No TM'
    return result


def compare_family(d):
    """
    compares if the pdb structure satisfies the family comparison (a CDB structure with the associated Uniprot)
    :param d: dictionary
    :return: binary
    """
    list_cdb = get_pores_from_cdb()
    list_cdb_up = []
    for i in range(len(list_cdb)):
        list_cdb_up.append(list_cdb[i].upper())
    for id in list_cdb_up:
        if id in d:
            return True
    return False


def compare_family_all(list_tmr, list_por):
    """
    finds if a structure has a family CDB structure for all the structures in the list
    :param list_tmr: list of tuples from xml files
    :param list_por: list of pdbids
    :return: dictionary pdbid:Binary, says if the structure has a family ChannelsDB structure
    """
    result = {}
    pupper = []
    for i in range(len(list_por)):
        pupper.append(list_por[i].upper())
        for j in range(len(list_tmr)):
            if list_por[i] in list_tmr[j][2]:
                if pupper[i] not in result:
                    result.update({pupper[i]: compare_family(list_tmr[j][2])})
                    # print(pupper[i])
    return result


def search_uniprot_line(list_tmr, pdbid):
    """
    finds the line with the uniprot containing the pdb structure
    :param list_tmr: list of tuples from xml files
    :param pdbid: pdb id
    :return: tuple
    """
    mapp = find_uniprot_from_sifts(pdbid.lower())
    l = []
    for t in list_tmr:
        if t[0] in mapp:
            l.append(t)
    return l


def filter_tm(t, pdbid, per, abs):
    """
    verifies if the structure approaches the entire TM part of the protein
    :param t: tuple
    :param pdbid: pdb id
    :return: binary
    """
    try:
        tm1 = get_first_tm(t)
        tm2 = get_last_tm(t)
        length = tm2 - tm1
        st1 = get_first_st(t, pdbid)
        st2 = get_last_st(t, pdbid)
        diff = 0
        if (st1 - tm1) > 0:
            diff += (st1 - tm1)
        if (tm2 - st2) > 0:
            diff += (tm2 - st2)
        if diff < (per*length + abs):
            return True
    except:
        print(str(t) + pdbid)
    return False


def compare_family_tm(t, pdbid):
    """
    compares the coordinates of the structure with a CDB structure
    :param t: tuple
    :param pdbid: pdb id
    :return: binary
    """
    res = False
    cdb_list = list_upper(get_pores_from_cdb())
    this_cdb = {}
    for key in t[2]:
        if key in cdb_list:
            this_cdb[key] = t[2][key]
    try:
        length = int(re.findall('\d+', t[2][pdbid])[-1]) - int(re.findall('\d+', t[2][pdbid])[0])
        len_cdb = find_min_position(this_cdb)
        if length >= len_cdb:
            res = True
    except KeyError:
        print(str(t) + pdbid)
    if res:
        this_str = len(find_chain(t[2][pdbid]))
        str_cdb = len(find_min_chain(this_cdb))
        if this_str >= str_cdb:
            return True
    return False


def find_min_position(d):
    """
    finds the smallest sequence of a residues
    :param d: dictionary
    :return: length of the sequence
    """
    diff_min = 0
    for key in d:
        st1 = int(re.findall('\d+', d[key])[0])
        st2 = int(re.findall('\d+', d[key])[-1])
        if (st2 - st1) < diff_min:
            diff_min = (st2 - st1)
        elif diff_min == 0:
            diff_min = (st2 - st1)
    return diff_min


def find_chain(s):
    """
    finds the letters (protein chain) in the string
    :param s: string
    :return: integer, chains of the protein
    """
    if s.find(',') == -1:
        return s[:s.find('=')]
    else:
        return find_chain(s[:s.find(',')]) + '/' + find_chain(s[(s.find(',')+2):])


def find_min_chain(d):
    """
    finds the smallest chain of a structure
    :param d: dictionary
    :return: string, chain sequence
    """
    str_min = ''
    for key in d:
        chain = find_chain(d[key])
        if len(chain) < len(str_min):
            str_min = chain
        elif str_min == '':
            str_min = chain
    return str_min


def check_tm(l, per, abs):
    """
    finds the pdb ids apt to have pores
    :param l: list of pdbids to check
    :return: list of pdbids to send to mole
    """
    unpr_list = se.get_unpr_list(get_all_uniprotid(l))

    download_xml(unpr_list)
    pdbid_to_mole = []; pdbid_to_mole_filter = []; pdbid_to_mole_filter_famtm = []
    pdbid_no_filter = []; pdbid_no_tm = []; pdbid_no_family = []; pdbid_no_uni = []
    list_data = load_all(unpr_list)
    print('Data loaded')
    l_upper = list_upper(l)
    tm_comparison = compare_all(list_data, l_upper)
    print(tm_comparison)
    family_comparison = compare_family_all(list_data, l_upper)
    print(family_comparison)
    for pdbid in l_upper:
        if pdbid not in tm_comparison:
            pdbid_no_uni.append(pdbid)
        elif tm_comparison[pdbid] == 'No TM':
            pdbid_no_tm.append(pdbid)
        else:
            if family_comparison[pdbid] and tm_comparison[pdbid]:
                pdbid_to_mole.append(pdbid)
            elif family_comparison[pdbid] and not tm_comparison[pdbid]:
                print(pdbid)
                lines = search_uniprot_line(list_data, pdbid)
                for line in lines:
                    if filter_tm(line, pdbid, per, abs):
                        pdbid_to_mole_filter.append(pdbid)  # possible plurality
                        break
                    else:
                        if compare_family_tm(line, pdbid):
                            pdbid_to_mole_filter_famtm.append(pdbid)  # possible plurality
                            break
                        else:
                            if pdbid not in pdbid_no_filter:
                                pdbid_no_filter.append(pdbid)  # possible plurality
            else:
                pdbid_no_family.append(pdbid)
    return pdbid_to_mole, pdbid_to_mole_filter, pdbid_to_mole_filter_famtm, pdbid_no_filter, pdbid_no_family, \
           pdbid_no_tm, pdbid_no_uni

# DB selection


def get_pores_from_cdb():
    """
    gets the current list of pores in the ChannelsDB
    :return: list of pores in ChannelsDB
    """
    pores = []
    if not os.path.isfile('content.json'):
        urllib.request.urlretrieve("https://webchem.ncbr.muni.cz/API/ChannelsDB/Content", 'content.json')
    with open('content.json', 'r') as f:
        content = json.load(f)
        for pdbid in content:
            if content[pdbid]['counts'][2] != 0:
                pores.append(pdbid)
    return pores


def get_data_from_cdb(l):
    """
    compares the channels db entries with the pdb entries in the list
    :param l: list of pdb ids
    :return: list of commun entries
    """
    pores = []
    if not os.path.isfile('content.json'):
        urllib.request.urlretrieve("https://webchem.ncbr.muni.cz/API/ChannelsDB/Content", 'content.json')
    with open('content.json', 'r') as f:
        content = json.load(f)
        for pdbid in l:
            if pdbid in content:
                pores.append(pdbid)
    return pores


def get_all_memprotmd_references(database):
    """
    commands all the references from the given database from the MemProtMD site
    :return: json file (list)
    """
    MEMPROTMD_ROOT_URI = "http://memprotmd.bioch.ox.ac.uk/"
    return requests.post(MEMPROTMD_ROOT_URI + "api/references/all/" + database).json()


def get_dict_classes(database):
    """
    returns a dictionary containing classes with its pdb ids from a given database
    :param database: database
    :return: dictionary of classes in a database
    """
    new_d = {}
    l = get_all_memprotmd_references(database)
    for d in l:
        new_d[d['accession']] = []  # 'accession/title'
        for item in d['simulations']:
            new_d[d['accession']].append(item[:4])  # 'accession/title'
    return new_d


def get_pores_from_db(d, db):
    """
    finds the pdbids with pores for a database
    :param d: dictionary
    :param db: database
    :return: dictionary of classes:pores
    """
    new_d = {}
    config = load_json("config.json")
    list_accession = []
    if db == 'TCDB' or db == 'mpm':
        list_accession = config['pores']['database'][db]['classes']
    elif db == 'mpstruc':
        list_accession = config['pores']['database'][db]['all']['classes']
    elif db == 'mpstruc-alpha':
        list_accession = config['pores']['database']['mpstruc']['alpha']['classes']
    elif db == 'mpstruc-beta':
        list_accession = config['pores']['database']['mpstruc']['beta']['classes']
    for key in d:
        if key in list_accession:
            new_d[key] = d[key]
    return new_d


def check_db():
    """
    finds the pdb ids apt to have pores in the mpm, mpstruc and tcdb databases
    :return: tuple: pdb ids apt to be send to MOLE, all the structures in the databases
    """
    structures = []
    list_db = ['mpm', 'TCDB', 'mpstruc']
    for db in list_db:
        list_pores = an.list_from_dict(get_pores_from_db(get_dict_classes(db), db))
        for pdbid in list_pores:
            if pdbid not in structures:
                structures.append(pdbid)
    cdb = get_pores_from_cdb()
    return an.no_intersection_list(structures, cdb), structures

# Additional


def list_upper(l):
    """
    transforms every string in the list on a uppercase string
    :param l: list of pdb ids
    :return: list of pdb ids
    """
    new_l = []
    for i in range(len(l)):
        new_l.append(l[i].upper())
    return new_l


def list_lower(l):
    """
    transforms every string in the list on a lowercase string
    :param l: list of pdb ids
    :return: list of pdb ids
    """
    new_l = []
    for i in range(len(l)):
        new_l.append(l[i].lower())
    return new_l


if __name__ == "__main__":

    with open('pores.txt') as f:  # get the list of pores
        my_pores = f.readlines()
        for i in range(len(my_pores)):
            my_pores[i] = my_pores[i].rstrip('\n')
    print(len(my_pores))
    # download_xml(my_pores) # works
    # my_uni  = [*get_all_uniprotid(my_pores)]

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
    print(len(PORES_UP))
    print(get_first_st(list_tm[0], '3S3W'))
    print(get_last_st(list_tm[0], '3S3W'))
    print(get_first_tm(list_tm[0]))
    print(get_last_tm(list_tm[0]))
    print(compare_tm_structure(list_tm[0], '3S3W'))

    no_complet_TM_cdb = []
    no_complet_TM_family = []

    my_di_str = compare_all(list_tm, PORES_UP)
    print(my_di_str)
    # print('nepreruseny chain:' + str(len(my_di_str)))
    tm_true_cdb = []
    for key in my_di_str:
        if my_di_str[key]:
            tm_true_cdb.append(key)
        else:
            no_complet_TM_cdb.append(key)
    print('komplet tm:' + str(len(tm_true_cdb)))

    pdb_keys = []
    for key in my_di_str:
        for k in list_pores:
            if key in k[2]:
                for pdbid in k[2]:
                    if pdbid not in pdb_keys:
                        pdb_keys.append(pdbid)
    print(pdb_keys)
    print('vsechny pdbid od komplet tm:' + str(len(pdb_keys)))
    print('pouze spriznene pdbid od komplet tm:' + str(len(an.no_intersection_list(pdb_keys, tm_true_cdb))))
    print('')
    pdbid_to_mole = check_db()[1]
    print(len(check_db()[0]), len(pdbid_to_mole))
    print(len(get_pores_from_cdb()))
    """
    with open('pdbid_to_mole.txt', 'w') as f:
        for item in check_db()[0]:
            f.write("%s\n" % item)
    """
    pores_tocdb_UPPER = []
    for i in range(len(pdbid_to_mole)):
        pores_tocdb_UPPER.append(pdbid_to_mole[i].upper())
    print(an.no_intersection_list(pdb_keys, pores_tocdb_UPPER))
    pdbid_to_mole_fromTM = an.no_intersection_list(pdb_keys, pores_tocdb_UPPER)
    """
    with open('pdbid_to_mole_fromTM.txt', 'w') as f:
        for item in an.no_intersection_list(pdb_keys, pores_tocdb_UPPER):
            f.write("%s\n" % item.lower())
    """
    family_st = an.no_intersection_list(pdb_keys, tm_true_cdb)
    my_di_str2 = compare_all(list_tm, family_st)
    print(my_di_str2)
    # print('pouze family pdbid od komplet tm s neprerusenym chainem:' + str(len(my_di_str2)))
    tm_to_mole = []
    for key in my_di_str2:
        if my_di_str2[key]:
            tm_to_mole.append(key)
        else:
            no_complet_TM_family.append(key)
    print('pouze family pdbid s komplet tm:' + str(len(tm_to_mole)))
    print('pouze family pdbid s komplet tm co nejsou v pores db:' + str(len(an.no_intersection_list(tm_to_mole, pores_tocdb_UPPER))))
    print(len(no_complet_TM_cdb))
    print(len(no_complet_TM_family))
    print(no_complet_TM_cdb)
    print(no_complet_TM_family)
    print(an.no_intersection_list(no_complet_TM_cdb, pores_tocdb_UPPER))
    print(len(an.no_intersection_list(no_complet_TM_family, pores_tocdb_UPPER)))


    ana_no_complet_TM_family = []
    for pdbid in an.no_intersection_list(no_complet_TM_family, pores_tocdb_UPPER):
        for i in list_pores:
            if pdbid in i[2]:
                if i not in ana_no_complet_TM_family:
                    ana_no_complet_TM_family.append(i)
    for i in range(len(ana_no_complet_TM_family)):
        print(ana_no_complet_TM_family[i])
    print(len(ana_no_complet_TM_family))
    

    print(len(tm_to_mole))
    with open('pores.txt') as f:  # get the list of pores
        my_pores_new = f.readlines()
        for i in range(len(my_pores_new)):
            my_pores_new[i] = my_pores_new[i].rstrip('\n')
    tm_comp = check_tm(list_lower(tm_to_mole), 0.1, 5)
    print(len(tm_to_mole))
    print('MOLE - TM + Family: ' + str(len(tm_comp[0])))
    print('MOLE - Filter + Family: ' + str(len(tm_comp[1])))
    print('MOLE - FamilyTM: ' + str(len(tm_comp[2])))
    print('No Filter or FamilyTM: ' + str(len(tm_comp[3])))
    print('No Family: ' + str(len(tm_comp[4])))
    print('No TM: : ' + str(len(tm_comp[5])))
    print('No UniProt: ' + str(len(tm_comp[6])))
    print(tm_comp)
    print(an.no_intersection_list(tm_comp[0], list_upper(pdbid_to_mole)))




