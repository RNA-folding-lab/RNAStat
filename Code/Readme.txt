RNAStat can be used to calculate structural information of RNA 3D structure(s) include four modules: (1) the radius of gyration (i.e., size) and shape; (2) the secondary structure motifs; (3) the geometry of base-pairing and base-stacking; (4) the distances between atoms;

operating environment：windows and python3
note：It can only run on Windows, Linux version will be updated in the future
If you don't have a PYTHon3 environment, run the exe file directly. 


1.Calculate the radius of gyration and shape parameter for RNA structures.
Note: PDB format, e.g., .cif
    1.)Calculate the size and shape of RNA structure.
run size_shape.exe file ---> Enter 0 to select a single file---> Enter the path of the RNA structure, such as D:\dataset\1b7f.cif.
    2.)Calculate the size and shape of RNA structures.
run size_shape.exe file ---> Enter 1 to select dataset---> Enter the path of the RNA structure dataset, such as D:\ dataset\---> Enter the path of results, such as 
D:\ save_path\---> output size_shape.csv, meaning of document see ReadMe of the file of Statistical results.

2.Obtain the secondary structure motifs.
    1.)Run cmd.exe to call the DSSR---> Please put DSSR and RNA structures in the same folder---> Output json file, which is analyzed by DSSR(http://x3dna.org/)
    2.)Please put the JSON file in a folder (no other file), such as, D:\json\---> Run motifs_stat.exe to obtain the secondary structure motifs---> Enter the json file path, such as D:\json\---> Enter the save_file path, such as D:\ ---> Output basepairs.csv and motif_stat.csv, meaning of document see ReadMe of the file of Statistical results.

3.Calculate the geometry of base-pairing and base-stacking.
Note: Please run motifs_stat.exe to get basepairs.csv file first
    1.)run base-pairing.exe ---> Enter the path of dataset, such as D:\ dataset\ ---> Enter the save_file path, such as D:\ ---> Enter the path of basepairs.csv---> output base-pairing.csv, meaning of document see ReadMe of the file of Statistical results.
    2.)run base-stacking.exe---> Enter the path of dataset, such as D:\ dataset\ ---> Enter the save_file path, such as D:\ ---> Enter the path of basepairs.csv---> output base-pairing.csv, meaning of document see ReadMe of the file of Statistical results.


4.Calculate the distances between atoms of RNA structure(s).
    1.)Calculate the distance between any two atoms
Run atom_distance.exe ---> Enter 1 to Calculate the distance between any two atoms ---> Enter 0 to select single file (or Enter 1 to select dataset) ---> Enter the path to store the file(.txt), such as D:/distance.txt---> Input the spacing of atomic distances(1-5), such as, 1: The distance between two atoms in the nearest neighbor nucleotides, 2 :The distance between two atoms in the second-nearest neighbor nucleotides and so on---> Enter the atoms name ---> Enter the path of RNA structure file(s), such as D:/dataset/1a9n.cif(D:/dataset/) ---> Output xx.txt, meaning of document see ReadMe of the file of Statistical results.
    2.)Calculate the distance between two atoms in base (A, U, G, C)
Run atom_distance.exe ---> Refer to the previous step
3.)Calculate the distance between atoms in different base pairs (A-U, G-C…)
Run atom_distance.exe ---> Refer to the previous step
