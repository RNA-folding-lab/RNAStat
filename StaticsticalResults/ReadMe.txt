The description of each file in this folder is as follows.
1. size and shape.csv: The calculation of the radius of gyration, asphericity parameter and shape parameter for RNA structures in the dataset.
1th column : PDB ID of  RNA structure ; 
2th column : the all length of RNA structure;
3th column : the radius of gyration(Rg) of RNA structure; 
4th column : the asphericity parameter ∆;
5th column :  the shape parameter S

2.motif_stat.csv: The statistics of secondary structural elements(stems,all loops,hairpin loops,bulge loops,internal loops and junction loops)
 for RNA structures in the dataset.
1th column : PDB ID of  RNA structure ; 
2th column : the num of all stems in RNA structure; 
3th column : the length of all stems(bp) in RNA structure;
4th column : the num of all hairpin loops  in RNA structure; 
5th column : the length of all hairpin loops(nt)  in RNA structure; 
6th column : the num of bulge loops  in RNA structure;  
7th column : the length of bulge loops(nt)  in RNA structure;  
8th column : the num of internal loops  in RNA structure; 
9th column : the length of internal loops(nt)  in RNA structure;
10th column : the num of junction loops in RNA structure；
11th column : the length of junction loops(nt)  in RNA structure.

3.basepairs.csv: The statistics information of base-pairing for RNA structures in the dataset.
1th column : PDB ID of  RNA structure ;  
2th column : the first base;
3th column : the second base;
4th column : the name of base-pairing; 
5th column and 6th column : the type of base-pairing.
.
4.Geometrical properties of base-pairing.csv:The calculation of the coordinates (ρ,θ,z) in base-pairing. The define of  (ρ,θ,z) see Figure 4(A).		
1th column : PDB ID of  RNA structure ; 
2th column and 3th column：Residue name of base pairing in PDB file ;
4th column and 5th column：Base name of base pairing ;
6th column：The distance between the centers of the two base planes ;
7th column : |z| for base in the coordinate system of its paired base  ;
8th column : ρ ;  
9th column : θ;
10th,11th and 12th:The types of base pairing.

5.Geometrical properties of base-stacking.csv:The calculation of the coordinates (ρ,θ,z) in base-stacking.
1th column : PDB ID of  RNA structure ; 
2th column and 3th column：Residue name of base-stacking ;
4th column and 5th column：Base name of base-stacking ;
6th column：The distance between the centers of the two base planes ;
7th column : |z| for base in the coordinate system of its paired base  ;
8th column : ρ ;  
9th column : θ;
10th,11th and 12th:The types of base pairing.

6.Ap-Up.txt: The distance statistics between P atoms in two paired nucleotides (A-U).

8.base-P-P.txt: The distance statistics between P atoms in  all types of base pairs.

9.N1C-N3G-basepairs.txt：The distance statistics between N3 and N1 atoms in two paired nucleotides (C-G).

10.N1C-N3G.txt:The distance statistics between N3 atoms in cytosine and N1 atoms in guanine.

11.P-2P.txt:The distance statistics between two P atoms in the second-nearest neighbor nucleotides.

12.PA-PG.txt:The distance statistics between P atoms in nucleotides with bases of A and G.

13.PA-PU.txt:The distance statistics between P atoms in nucleotides with bases of A and U.

14.PC4-.txt: The distance statistics between P and C4' atoms.

15.PP.txt:The distance statistics between P atoms.

16.P-P.txt:The distance statistics between two P atoms in the second-nearest neighbor nucleotides.

17.N1G-N3U.txt:The distance statistics between N1 and N3 atoms in two paired nucleotides (G-U).

18.O2'G-C3'U.txt:The distance statistics between O2' atoms in guanine and C3' atoms in uracil.

19.N1A-N3U.txt:The distance statistics between N1 atoms in adenine and N3 atoms in uracil.

20.C2A-O2U.txt:The distance statistics between C2 atoms in adenine and O2 atoms in uracil.

21.Parameter file of base-pairing.csv:Parameter file of Gaussian density distribution of ρ and theta
