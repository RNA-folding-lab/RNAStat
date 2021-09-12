import numpy as np
import os
import csv
import re


A = ['P', 'OP1', 'OP2', '\"O5\'"', '\"C5\'"', '\"C4\'"', '\"O4\'"', '\"C3\'"', '\"O3\'"', '\"C2\'"', '\"O2\'"',
     'N9', 'C8', 'N7', 'C5', 'C6', 'N6', 'N1', 'C2', 'N3', 'C4']
U = ['P', 'OP1', 'OP2', '\"O5\'"', '\"C5\'"', '\"C4\'"', '\"O4\'"', '\"C3\'"', '\"O3\'"', '\"C2\'"', '\"O2\'"',
     'N1', 'C2', 'O2', 'N3', 'C4', 'O4', 'C5', 'C6']
G = ['P', 'OP1', 'OP2', '\"O5\'"', '\"C5\'"', '\"C4\'"', '\"O4\'"', '\"C3\'"', '\"O3\'"', '\"C2\'"', '\"O2\'"',
     'N9', 'C8', 'N7', 'C5', 'C6', 'O6', 'N1', 'C2', 'N2', 'N3', 'C4']
C = ['P', 'OP1', 'OP2', '\"O5\'"', '\"C5\'"', '\"C4\'"', '\"O4\'"', '\"C3\'"', '\"O3\'"', '\"C2\'"', '\"O2\'"',
     'N1', 'C2', 'O2', 'N3', 'C4', 'N4', 'C5', 'C6']
element = ['A', 'U', 'G', 'C']


def get_basepair_coord(FilePath, chain, id, residue, Atom):
    print(FilePath, chain, id, residue, Atom)
    coord = np.array([0, 0, 0])
    with open(FilePath, 'r') as f:
        for line in f:
            start = line.split()
            try:
                if start[-3] == chain and start[-5] ==id and start[-4] ==residue and start[3] == Atom:
                    print(line)
                    x = round(float(start[10]), 2)
                    y = round(float(start[11]), 2)
                    z = round(float(start[12]), 2)
                    coord = np.array([x, y, z], dtype=float)
            finally:
                continue
    return coord


def get_BasePair_distance(save_path, csv_path):
    base1 = input('输入碱基: ')
    atom1 = input('原子名称: ')
    base2 = input('')
    atom2 = input('')
    DatasetPath = input('数据集存放路径: ')
    csv_file = open(csv_path, 'r')
    reader = csv.reader(csv_file)
    print(reader)
    for item in reader:
        print(item)
        pdb_id = item[0]
        base_1 = item[1]
        base_2 = item[2]
        bp = item[3]
        if bp == base1 + base2:
            i = pdb_id + '.cif'
            FilePath = os.path.join(DatasetPath, i)
            chain1 = base_1.split('.')[0]
            chain2 = base_2.split('.')[0]

            id1 = re.findall(r"\d+\.?\d*", base_1.split('.')[1])[0]
            id2 = re.findall(r"\d+\.?\d*", base_2.split('.')[1])[0]

            residue1 = base_1.split('.')[1].replace(id1, '')
            residue2 = base_2.split('.')[1].replace(id2, '')

            coord1 = get_basepair_coord(FilePath, chain1, id1, residue1, atom1)
            coord2 = get_basepair_coord(FilePath, chain2, id2, residue2, atom2)
            if all(coord1 != [0, 0, 0]) and all(coord2 != [0, 0, 0]):
                distance = round(np.linalg.norm(coord1 - coord2), 2)
                print(distance)
                with open(save_path, 'a') as f:
                    if distance != 0:
                        print(distance)
                        f.write(str(distance))
                        f.write('\n')


def get_BaseAtom_coord(file_path, residue, Atom):
    Coords = []
    with open(file_path, 'r') as file:
        for line in file:
            start = line.split()
            if len(start) < 1:
                continue
            if start[0] == 'ATOM' and start[-4] ==residue and start[3] == Atom:
                # atom_id = start[6] + '-' + start[5] + '-' + start[2]
                x = round(float(start[10]), 2)
                y = round(float(start[11]), 2)
                z = round(float(start[12]), 2)
                coord = np.array([x, y, z], dtype=float)
                # print(coord)
                Coords.append(coord)
    CoordsArray = np.array(Coords)
    return CoordsArray



def get_coord(file_path, name):
    Coords = []
    with open(file_path, 'r') as file:
        for line in file:
            start = line.split()
            if len(start) < 1:
                continue
            if start[0] == 'ATOM' and start[3] == name and start[-4] in element:
                x = round(float(start[10]), 2)
                y = round(float(start[11]), 2)
                z = round(float(start[12]), 2)
                coord = np.array([x, y, z], dtype=float)
                # print(coord)
                Coords.append(coord)
    CoordsArray = np.array(Coords)
    return CoordsArray


def calc_adjacent_distance(save_path, coord, Count_interval):
    c_list = list(range(len(coord)))
    for i in c_list:
        if i + Count_interval < len(coord):
            distance = round(np.linalg.norm(coord[i + Count_interval]-coord[i]), 2)
            with open(save_path, 'a') as f:
                if distance != 0:
                    print(distance)
                    f.write(str(distance))
                    f.write('\n')


def calc_distance(save_path, coord1, coord2, Count_interval):
    c1_list = list(range(len(coord1)))
    c2_list = list(range(len(coord2)))

    for i in c1_list:
        for j in c2_list:
            if j + Count_interval < len(coord2):
                distance = round(np.linalg.norm(coord2[j + Count_interval]-coord1[i]), 2)
                with open(save_path, 'a') as f:
                    if distance != 0:
                        print(distance)
                        f.write(str(distance))
                        f.write('\n')


def ChoiceCont():
    ChoiceCout = input('coutinue? y/n ').upper()
    if ChoiceCout == 'Y':
        main()
    if ChoiceCout == 'N':
        quit()


def main():
    while True:
        print('---------------------------------')
        print('1.Calculate the distance between any two atoms')
        print('2.Calculate the distance between two atoms in base(A, U, G, C)')
        print('3.Calculate the distance  between  atoms in different base pairs(A-U, G-C...)')
        print('''
        A : P, OP1, OP2, "O5'", "C4'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'", N9, C8, N7, C5, C6, N6, N1, C2, N3, C4
        U : P, OP1, OP2, "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'", N1, C2, O2, N3, C4, O4, C5, C6
        G : P, OP1, OP2, "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", “O2'”, "C1'", N9, C8, N7, C5, C6, O6, N1, C2, N2, N3, C4
        C : P, OP1, OP2, "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", ”O2'", "C1'", N1, C2, O2, N3, C4, N4, C5, C6
            ''')
        print('---------------------------------')

        choice = int(input('Enter 1,2 or 3 to Calculate the distance between different types of atoms---> :  '))
        choice_file = int(input('Enter 0 to select a single file/Enter 1 to select a dataset---> : '))
        SavePath = input('Enter the path to store the file(.txt), such as D:/distance.txt---> :')

        if choice == 1:
            interval_count = int(input('''Input the spacing of atomic distances(1-5), such as,
                                       1 : The distance  between two atoms in the nearest neighbor nucleotides
                                       2 :  The distance  between two atoms in the second-nearest neighbor nucleotides
                                       and so on---> : '''))
            atom1 = input('Enter the first atom name---> : ').upper()
            atom2 = input('Enter the second atom name---> : ').upper()

            if choice_file == 0:
                FilePath = input('Enter the path of RNA structure file, such as D:/dataset/1a9n.cif---> : ')
                if atom1 == atom2:
                    Coord = get_coord(FilePath, atom1)
                    if Coord.size > 0:
                        calc_adjacent_distance(SavePath, Coord, interval_count)
                        ChoiceCont()
                    else:
                        print('Please whether the entered atom name is correct ')
                        print('Check whether the entered file path is correct ')
                        main()
                else:
                    Coord1 = get_coord(FilePath, atom1)
                    Coord2 = get_coord(FilePath, atom2)
                    print(Coord1, Coord2)
                    if Coord1.size > 0 and Coord2.size > 0:
                        calc_distance(SavePath, Coord1, Coord2, interval_count)
                        ChoiceCont()
                    else:
                        print('Please whether the entered atom name is correct ')
                        print('Check whether the entered file path is correct ')
                        main()

            if choice_file == 1:
                DatasetPath = input('Enter the path of dataset, such as D:/dataset---> : ')
                file_list = os.listdir(DatasetPath)
                for i in file_list:
                    FilePath = os.path.join(DatasetPath, i, )
                    print(FilePath)
                    if atom1 == atom2:
                        Coord = get_coord(FilePath, atom1)
                        calc_adjacent_distance(SavePath, Coord, interval_count)
                    else:
                        Coord1 = get_coord(FilePath, atom1)
                        Coord2 = get_coord(FilePath, atom2)
                        calc_distance(SavePath, Coord1, Coord2, interval_count)
                ChoiceCont()

        if choice == 2:
            interval_count = int(input('''Input the spacing of atomic distances(1-5), such as,
                                       1 : The distance  between two atoms in the nearest neighbor nucleotides
                                       2 :  The distance  between two atoms in the second-nearest neighbor nucleotides
                                       and so on---> : '''))
            residue1 = input('Enter the first base(A, U, G, C)---> : ').upper()
            Atom1 = input('Enter the first atom name---> : ').upper()
            residue2 = input('Enter the second base(A, U, G, C)---> : ').upper()
            Atom2 = input('Enter the second atom name---> : ').upper()
            if choice_file == 0:
                FilePath = input('Enter the path of RNA structure file, such as D:/dataset/1a9n.cif---> : ')
                coord1 = get_BaseAtom_coord(FilePath, residue1, Atom1)
                coord2 = get_BaseAtom_coord(FilePath, residue2, Atom2)
                calc_distance(SavePath, coord1, coord2, interval_count)
                ChoiceCont()
            if choice_file == 1:
                DatasetPath = input('Enter the path of dataset, such as D:/dataset---> : ')
                file_list = os.listdir(DatasetPath)
                for i in file_list:
                    FilePath = os.path.join(DatasetPath, i)
                    coord1 = get_BaseAtom_coord(FilePath, residue1, Atom1)
                    coord2 = get_BaseAtom_coord(FilePath, residue2, Atom2)
                    calc_distance(SavePath, coord1, coord2, interval_count)
                ChoiceCont()

        if choice == 3:
            if choice_file == 1:
                print('note: Please run motifs_stat.exe to obtain the statistics data of base-pairing')
                basepairs_save_path = input('Enter the path of statistics file of base-pairing---> :  ')
                get_BasePair_distance(SavePath, basepairs_save_path)
                ChoiceCont()
            if choice_file == 1:
                print('Calculating the distance between atoms in base-pairing of individual file is not supported')
                print('Please select again')
                ChoiceCont()

        else:
            break


if __name__ == '__main__':
    main()
