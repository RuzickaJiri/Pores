import os
import shutil


def mole_basic(start, end):
    with open('pores.txt') as f:  # get the list of pores
        my_pores_all = f.readlines()
        for i in range(len(my_pores_all)):
            my_pores_all[i] = my_pores_all[i].rstrip('\n')
    print(my_pores_all)

    pores_trunc = my_pores_all[start:end]

    for pdbid in pores_trunc:
        try:
            os.mkdir("C:\Computations\Calc\\" + pdbid)

            os.popen("C:\Computations\Pores\Pores.exe ""C:\Computations\Calc\\" + pdbid + ".json")

        except FileExistsError:
            print(pdbid)


def mole_tm(start, end):
    with open('pores.txt') as f:  # get the list of pores
        my_pores_all = f.readlines()
        for i in range(len(my_pores_all)):
            my_pores_all[i] = my_pores_all[i].rstrip('\n')
    print(my_pores_all)

    pores_trunc = my_pores_all[start:end]

    for pdbid in pores_trunc:
        try:
            os.mkdir("C:\Computations\Calc\\" + pdbid + "_TM")

            os.popen("C:\Computations\Pores\Pores.exe ""C:\Computations\Calc\\" + pdbid + "_TM.json")

        except FileExistsError:
            print(pdbid + "_TM")


def get_exist(path):
    return os.path.exists(path + "\\json")


def get_list_exist(l):
    does_not = 0
    pdb_list = []
    for pdbid in l:
        if not get_exist("C:\Computations\Calc\\" + pdbid):
            does_not += 1
            pdb_list.append(pdbid)
    return does_not, pdb_list


def get_list_exist_tm(l):
    does_not = 0
    pdb_list = []
    for pdbid in l:
        if not get_exist("C:\Computations\Calc\\" + pdbid + "_TM"):
            does_not += 1
            pdb_list.append(pdbid)
    return does_not, pdb_list


with open('pores.txt') as f:  # get the list of pores
    my_pores_all = f.readlines()
    for i in range(len(my_pores_all)):
        my_pores_all[i] = my_pores_all[i].rstrip('\n')
print(my_pores_all)
print(get_exist("C:\Computations\Calc\\3s3w"))
print(get_list_exist(my_pores_all))
print(get_list_exist_tm(my_pores_all))
none_basic = get_list_exist(my_pores_all)[1]
none_tm = get_list_exist_tm(my_pores_all)[1]
i = 0
for p in none_basic:
    if p in none_tm:
        print(p)


