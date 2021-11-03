import pandas as pd
import re
import numpy as np
import os
import csv
from pathlib import Path


pur = ['N1', 'C2', 'N3', 'C4', 'C5', 'C6', 'N7', 'C8', 'N9']
pyr = ['N1', 'C2', 'N3', 'C4', 'C5', 'C6']



def get_coord(filename, chain_name, base_name, id_name):
    """Extract the atomic coordinates of each base"""
    # print(file)
    # print(residue.get_chains())
    # print(residue.get_residue_id())
    coords = {}
    with open(filename, 'r') as f:
        for line in f:
            start = line.split()
            try:
                if start[-3] == chain_name and start[-4] == base_name and start[-5] == id_name:
                    # print(line)
                    x = round(float(start[10]), 2)
                    y = round(float(start[11]), 2)
                    z = round(float(start[12]), 2)
                    coord = np.array([x, y, z], dtype=float)
                    coords[start[3]] = coord
            finally:
                continue
    # print(coords)
    return coords


def Centroid_coord(base, coords):
    """
    Find the centroid coordinates of the base plane
    Determine purine or pyrimidine
    """
    cent_coor = np.array([])
    cent_coord_list = []

    if base == 'A' or base == 'G' or 'A' in base or 'G' in base:
        for i in pur:
            if i in coords:
                cent_coord_list.append(coords[i])
    elif base == 'C' or base == 'U' or 'C' in base or 'U' in base:
        for i in pyr:
            if i in coords:
                cent_coord_list.append(coords[i])
    cent_coord_list = np.array(cent_coord_list)
    if cent_coord_list.size > 0:
        cent_coor = cent_coord_list.mean(axis=0)

    return cent_coor


def base_stack(base, coords, cent_coor1, cent_coor2):
    L, bases_Z, bases_D, theta = 0, 0, 0, 0

    if base == 'A' or base == 'G':
        if 'N1' in coords.keys() and 'C8' in coords.keys() and 'N3' in coords.keys():
            N1 = coords['N1']
            C8 = coords['C8']
            N3 = coords['N3']
            A = N1[1] * (C8[2] - N3[2]) + C8[1] * (N3[2] - N1[2]) + N3[1] * (N1[2] - C8[2])
            B = N1[2] * (C8[0] - N3[0]) + C8[2] * (N3[0] - N1[0]) + N3[2] * (N1[0] - C8[0])
            C = N1[0] * (C8[1] - N3[1]) + C8[0] * (N3[1] - N1[1]) + N3[0] * (N1[1] - C8[1])
            plane_D = N1[0] * (N3[1] * C8[2] - C8[1] * N3[2]) + C8[0] * (N1[1] * N3[2] - N3[1] * N1[2]) + \
                      N3[0] * (C8[1] * N1[2] - N1[1] * C8[2])

            n = np.array([A, B, C])  
            t = (A * cent_coor2[0] + B * cent_coor2[1] + C * cent_coor2[2] + plane_D) / (A * A + B * B + C * C)
            point_coord = np.array([cent_coor2[0] - A * t, cent_coor2[1] - B * t, cent_coor2[2] - C * t])
            bases_Z = np.linalg.norm(cent_coor2 - point_coord)
            bases_D = np.linalg.norm(cent_coor1 - point_coord)
            L = np.linalg.norm(cent_coor1 - cent_coor2)
            V1 = np.array([N1[0] - cent_coor1[0], N1[1] - cent_coor1[1], N1[2] - cent_coor1[2]])
            V2 = np.array(
                [point_coord[0] - cent_coor1[0], point_coord[1] - cent_coor1[1], point_coord[2] - cent_coor1[2]])
            V1_mo = np.sqrt(V1.dot(V1))
            V2_mo = np.sqrt(V2.dot(V2))
            vector_product = V1.dot(V2)
            cos_ = vector_product / (V1_mo * V2_mo)
            angle_hu = np.arccos(cos_)
            angle = angle_hu * 180 / np.pi
            cross_product = np.cross(V1, V2)
            direction = cross_product.dot(n)
            if direction >= 0:
                theta = angle
            else:
                theta = 2 * 180 - angle

    elif base == 'U' or base == 'C':
        if 'N1' in coords.keys() and 'N3' in coords.keys() and 'C4' in coords.keys():
            N1 = coords['N1']
            N3 = coords['N3']
            C4 = coords['C4']
            A = N1[1] * (N3[2] - C4[2]) + N3[1] * (C4[2] - N1[2]) + C4[1] * (N1[2] - N3[2])
            B = N1[2] * (N3[0] - C4[0]) + N3[2] * (C4[0] - N1[0]) + C4[2] * (N1[0] - N3[0])
            C = N1[0] * (N3[1] - C4[1]) + N3[0] * (C4[1] - N1[1]) + C4[0] * (N1[1] - N3[1])
            plane_D = N1[0] * (C4[1] * N3[2] - N3[1] * C4[2]) + \
                N3[0] * (N1[1] * C4[2] - C4[1] * N1[2]) + \
                C4[0] * (N3[1] * N1[2] - N1[1] * N3[2])
            n = np.array([A, B, C])
            t = (A * cent_coor2[0] + B * cent_coor2[1] + C * cent_coor2[2] + plane_D) / (A * A + B * B + C * C)
            point_coord = np.array([cent_coor2[0] - A * t, cent_coor2[1] - B * t, cent_coor2[2] - C * t])
            bases_Z = np.linalg.norm(cent_coor2 - point_coord)
            bases_D = np.linalg.norm(cent_coor1 - point_coord)
            L = np.linalg.norm(cent_coor1 - cent_coor2)
            V1 = np.array([N3[0] - cent_coor1[0], N3[1] - cent_coor1[1], N3[2] - cent_coor1[2]])
            V2 = np.array(
                [point_coord[0] - cent_coor1[0], point_coord[1] - cent_coor1[1], point_coord[2] - cent_coor1[2]])
            V1_mo = np.sqrt(V1.dot(V1))
            V2_mo = np.sqrt(V2.dot(V2))
            vector_product = V1.dot(V2)
            cos_ = vector_product / (V1_mo * V2_mo)
            angle_hu = np.arccos(cos_)
            angle = angle_hu * 180 / np.pi
            cross_product = np.cross(V1, V2)
            direction = cross_product.dot(n)
            if direction >= 0:
                theta = angle
            else:
                theta = 2 * 180 - angle

    # print(L, bases_Z, bases_D, theta)
    return round(L, 3), round(bases_Z, 3), round(bases_D, 3), round(theta, 3)


def get_base_name(base):
    if base == 'A' or 'A' in base:
        base = 'A'
    elif base == 'U' or 'U' in base:
        base = 'U'
    elif base == 'G' or 'G' in base:
        base = 'G'
    elif base == 'C' or 'C' in base:
        base = 'C'
    return base


def calc_basestacking(root, csv_path, save_path):
    path = Path(csv_path)
    data = pd.read_csv(path, names=['name', 'base1', 'base2', 'bp', 'bp_name', 'LW'])

    base_1 = data[['name', 'base1', 'bp_name']]
    base_2 = data[['name', 'base2', 'bp_name']]

    def calc(base_data, base_col):
        ranges = list(range(0, base_data.size - 1))
        for i in ranges:
            if i < (base_data.size) / 3 - 1:
                d1 = base_data.loc[[i]]
                d2 = base_data.loc[[i + 1]]

                name1 = d1['name'][i]
                name2 = d2['name'][i + 1]

                r1 = d1[base_col][i]
                r2 = d2[base_col][i + 1]

                bp_name1 = d1['bp_name'][i]
                bp_name2 = d2['bp_name'][i + 1]

                C1 = r1.split('.')[0]
                C2 = r2.split('.')[0]
                try:
                    id1 = re.findall('\d+$', r1.split('.')[1])[0]
                    id2 = re.findall('\d+$', r2.split('.')[1])[0]
                    if name1 == name2 and C1 == C2 and int(id1) == int(id2) - 1 and bp_name1 in ['WC', 'Wobble'] \
                            and bp_name2 in ['WC', 'Wobble']:
                        file = os.path.join(root, name1 + '.cif')
                        b1 = r1.split('.')[1].replace(id1, '')
                        b2 = r2.split('.')[1].replace(id2, '')
                        # print(name1, r1, C1, id1, b1)
                        coor1 = get_coord(file, C1, b1, id1)
                        coor2 = get_coord(file, C2, b2, id2)

                        centcoor1 = Centroid_coord(b1, coor1)
                        centcoor2 = Centroid_coord(b2, coor2)
                        # print(centcoor1, centcoor2)

                        b_1 = get_base_name(b1)
                        b_2 = get_base_name(b2)

                        L, Z, D, Theta = base_stack(b_1, coor1, centcoor1, centcoor2)
                        data = [name1, r1, r2, b_1, b_2, L, Z, D, Theta]
                        if L != 0 and Z != 0 and D != 0 and Theta != 0:
                            csvPath = os.path.join(save_path, 'base-stacking.csv')
                            with open(csvPath, "a") as csvfile:
                                csv_writer = csv.writer(csvfile, lineterminator='\n')
                                csv_writer.writerow(data)
                finally:
                    continue

    calc(base_1, 'base1')
    calc(base_2, 'base2')


def main():
    dataset_path = input('Enter the path of dataset-----: ')
    basepairs_save_path = input('Enter the path of statistics file of base-pairing---> :  ')
    save_path = input('Enter the path of results---> :  ')
    calc_basestacking(dataset_path, basepairs_save_path, save_path)
    choice = input('continue? Y/N  ').upper()
    print('-----------------')
    if choice == 'Y':
        main()
    if choice == 'N':
        quit()


if __name__ == '__main__':
    main()