import os
import shutil

with open('pores.txt') as f:  # get the list of pores
    my_pores_all = f.readlines()
    for i in range(len(my_pores_all)):
        my_pores_all[i] = my_pores_all[i].rstrip('\n')
print(my_pores_all)

pores_trunc = my_pores_all[30:40]

for pdbid in pores_trunc:
    try:
        os.mkdir("C:\Computations\Calc\\" + pdbid)

        os.popen("C:\Computations\Pores\Pores.exe ""C:\Computations\Calc\\" + pdbid + ".json")

        """
        shutil.move("C:\Computations\Calc\csv", "C:\Computations\Calc\\" + pdbid + "\csv")
        shutil.move("C:\Computations\Calc\chimera", "C:\Computations\Calc\\" + pdbid + "\chimera")
        shutil.move("C:\Computations\Calc\json", "C:\Computations\Calc\\" + pdbid + "\json")
        shutil.move("C:\Computations\Calc\pdb", "C:\Computations\Calc\\" + pdbid + "\pdb")
        shutil.move("C:\Computations\Calc\pymol", "C:\Computations\Calc\\" + pdbid + "\pymol")
        shutil.move("C:\Computations\Calc\\vmd", "C:\Computations\Calc\\" + pdbid + "\\vmd")
        shutil.move("C:\Computations\Calc\\xml", "C:\Computations\Calc\\" + pdbid + "\\xml")
        """
    except FileExistsError:
        print(pdbid)
