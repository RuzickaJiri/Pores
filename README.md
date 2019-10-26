# Pores

## About
Internship research work on detecting and anatomy analysis of pores in transmembrane proteins. Based on Mole 2.5 computations and channelsDB database. 

### Selection
A priori search for pore structures took place on MemProt API and in its 3 databases - tcdb, mpmstruc and memprotDB. This gave a larger collection of pore structures and its UniProt ids. Selection parses a new PDB release ans looks for structures with the same UniProt ids. These structures pass controls (TM part) and if given conditions are satisfied, the structure is evaluated as susceptible to be added in channelsDB.

### Anatomy 
Anatomy analysis explores the statistics behind transmembrane pores - properties (such as polarity, charge, hydrophobicity etc.), length, residues (all, bottleneck) and detailed analysis of some proteins pore classes


## Files description

### Python files

<ul>
<li>analysis.py - contains methods for anatomy analysis</li>
<li>batch.py - contains methods to execute and get Mole computations </li>
<li>content.py - gets the content of channelsDB</li>
<li>data_analysis.py - main script for anatomy analysis</li>
<li>script.py - contains methods for pores selection</li>
<li>selection.py 	- main script for pores selection</li>
</ul> 

### Other files

<ul>
<li>Content.txt - the content of channelsDB </li>
<li>config.json - configuration file with a list of classes from MemProtDB containing pores </li>
<li>pdbid_to_mole.txt - pdbids of structures to by added in the database</li>
<li>pdbid_to_mole_fromTM.txt - pdbids of structures to by added in the database</li>
<li>template.json - template of json file that defines Mole 2.5 computations</li>
</ul> 
