#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 2019

@author: jruzicka
"""
import os

key_words = {"PORE", "CHANNEL", "MEMBRANE", "TRANSMEMBRANE"}

def check_header(file):
    with open(file) as input_file:
        lines = input_file.readlines()
    res = ""
    header = lines[0]
    classification = header[10:49]
    words = classification.split()
    for kw in key_words:
        for w in words:
            if kw == w:
                res += w + " "
    return res

def check_title(file):
    with open(file) as input_file:
        lines = input_file.readlines()
    res = ""
    i = 1
    while lines[i][0:5] == "TITLE":
        title = lines[i]
        classification = title[10:79]
        words = classification.split()
        for kw in key_words:
            for w in words:
                if kw == w:
                    res += w + " "
        i += 1
    return res

def check_keywds(file):
    with open(file) as input_file:
        lines = input_file.readlines()
    res = ""
    i = 0
    if i < len(lines):
        while lines[i][0:6] != "KEYWDS":
            i += 1
    if i < len(lines):
        while lines[i][0:6] == "KEYWDS":
            title = lines[i]
            classification = title[10:78]
            words = classification.split()
            for j in range(len(words)):
                if words[j][-1] == ",":
                    words[j] = words[j][:-1]
            for kw in key_words:
                for w in words:
                    if kw == w:
                        res += w + " "
            i += 1
    return res

def not_empty(str):
    if str != "":
        return True
    return False

def check_all(file):
    score = 0
    header = check_header(file)
    if not_empty(header):
        score += 20
    title = check_title(file)
    if not_empty(title):
        score += 20
    keywds = check_keywds(file)
    if not_empty(keywds):
        score += 20
    if score >= 20:
        return True
    return False

def check_release(path):
    new_release = os.listdir(path)
    dict = {}
    #i = 1
    for pdb_file in new_release:
        value = check_all(pdb_file)
        dict[pdb_file] = value
        #print(i)
        #i += 1
    return dict


# Tests
print(check_header("6g8z.pdb"))
print(check_title("6g8z.pdb"))
print(check_keywds("6g8z.pdb"))
print(check_all("6g8z.pdb"))
#print(check_all("6ji0-pdb-bundle1"))

# New release
new_release = os.listdir("Releases/1707/unzipped/pdb")
print(new_release)

# dictionary with pdbid.pdb : boolean
# value determines if the research was successful or not
rel = check_release("Releases/1707/unzipped/pdb")
print(rel)
nb_pores = 0
for item in rel:
    if rel[item]:
        nb_pores += 1
print(nb_pores)

# writes a file with the result of the research
with open('release0717.txt', 'w') as f:
    for item in rel:
        f.write(item)
        f.write(" ")
        f.write("%s\n" % str(rel[item]))
    f.write(str(nb_pores))
