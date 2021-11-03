import numpy as np
import os
import csv


element = ['A', 'U', 'G', 'C']
atom = ['P', 'OP1', 'OP2', '\"O5\'"', '\"C5\'"', '\"C4\'"', '\"O4\'"', '\"C3\'"', '\"O3\'"', '\"C2\'"', '\"O2\'"',
        '\"C1\'"', 'N9', 'C8', 'N7', 'C5', 'C6', 'O6', 'N1', 'C2', 'N2', 'N3', 'C4', 'O2', 'N4', 'O4', 'N6']


def chains(file_name):
    """
    How many chains are there in the structure?
    """
    with open(file_name, 'r') as file:
        chain_name = -1
        chains_list = []
        for line in file:
            start = line.split()
            if len(start) < 1:
                continue
            if start[0] in ['ATOM', 'HETATM'] and start[-4] in element:
                if start[-3] != chain_name:
                    # print(line)
                    chain_name = start[-3]
                    chains_list.append(chain_name)
    chain_count = len(chains_list)

    return chains_list


def sequence(chain, file_name):
    """
    The sequence
    Modifying bases cannot be identified
    """
    residues_list = []
    residues_num = -1
    with open(file_name, 'r') as file:
        for line in file:
            start = line.split()
            if len(start) < 1:
                continue
            if start[0] in ['ATOM', 'HETATM'] and start[-3] == chain:
                if start[5] in element:
                    if start[8] != residues_num:
                        residues_list.append(start[5])
                        residues_num = start[8]
                elif 'A' in start[5] and start[3] in atom:
                    if start[8] != residues_num:
                        residues_list.append('a')
                        residues_num = start[8]
                elif 'U' in start[5] and start[3] in atom:
                    if start[8] != residues_num:
                        residues_list.append('u')
                        residues_num = start[8]
                elif 'G' in start[5] and start[3] in atom:
                    if start[8] != residues_num:
                        residues_list.append('g')
                        residues_num = start[8]
                elif 'C' in start[5] and start[3] in atom:
                    if start[8] != residues_num:
                        residues_list.append('c')
                        residues_num = start[8]
            seq = ''.join(residues_list)
    chain_nt = len(seq)

    return chain_nt


def get_coord(chain, file_name):
    coor = []
    with open(file_name, 'r') as file:
        for line in file:
            start = line.split()
            if len(start) < 1:
                continue
            if start[0] in ['ATOM', 'HETATM'] and start[-3] == chain and start[3] in atom:
                # print(line)
                x_coor = float(start[10])
                y_coor = float(start[11])
                z_coor = float(start[12])
                ATOM_coor = np.array([x_coor, y_coor, z_coor])
                coor.append(ATOM_coor)
    return coor


def calc_Rg_func(coordlist):
    coordlist = np.array(coordlist)
    N = coordlist.size / 3
    # print(N)
    Rg, derta, s = 0, 0, 0
    try:
        cent_coor = coordlist.mean(axis=0)
        D = 0
        D_x = 0
        D_y = 0
        D_z = 0
        for atomcoord in coordlist:
            D += np.square(atomcoord[0] - cent_coor[0]) + np.square(atomcoord[1] - cent_coor[1]) + \
                 np.square(atomcoord[2] - cent_coor[2])
            D_x += np.square(atomcoord[0] - cent_coor[0])
            D_y += np.square(atomcoord[1] - cent_coor[1])
            D_z += np.square(atomcoord[2] - cent_coor[2])
        Rg2 = D / N
        Rg = np.sqrt(Rg2)
        l1 = round((D_x / N), 3)
        l2 = round((D_y / N), 3)
        l3 = round((D_z / N), 3)
        li = (l1 + l2 + l3) / 3
        derta = ((np.square(l1 - li) + np.square(l2 - li) + np.square(l3 - li)) / (Rg2 ** 2)) * (3 / 2)
        s = (((l1 - li) * (l2 - li) * (l3 - li)) / (Rg2 ** 3)) * 27
        # print(Rg, Rg2, l1 + l2 + l3, derta, s)

        return round(Rg, 3), round(derta, 3), round(s, 3)
    finally:
        return round(Rg, 3), round(derta, 3), round(s, 3)


def calc(File_path):
    nt = 0
    coord = []
    Chain_list = chains(File_path)
    for c in Chain_list:
        nt += sequence(c, File_path)
        coord += get_coord(c, File_path)
    RG, Derta, S = calc_Rg_func(coord)
    print('Rg = {}'.format(RG), 'Derta = {}'.format(Derta), 'S = {}'.format(S))
    return nt, RG, Derta, S


def main():
    choice_file = int(input('Enter 0 to select a single file/Enter 1 to select a dataset-----: '))
    if choice_file == 0:
        file_path = input('Enter the path of the RNA structure : ')
        calc(file_path)
        choice = input('continue? Y/N').upper()
        print('-----------------')
        if choice == 'Y':
            main()
        if choice == 'N':
            quit()
    if choice_file == 1:
        dataset_path = input('Enter the path of the RNA structure dataset :  ')
        save_path = input('Enter the path of results :  ')
        file_list = os.listdir(dataset_path)
        for i in file_list:
            file_path = os.path.join(dataset_path, i)
            N, RG, Derta, S = calc(file_path)
            if RG != 0:
                data = [i[0:4], N, RG, Derta, S]
                path = os.path.join(save_path, 'size_shape.csv')
                with open(path, "a") as csvfile:
                    CsvWriter = csv.writer(csvfile, lineterminator='\n')
                    CsvWriter.writerow(data)
        choice = input('continue? Y/N').upper()
        print('-----------------')
        if choice == 'Y':
            main()
        if choice == 'N':
            quit()


if __name__ == '__main__':
    main()





