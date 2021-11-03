import numpy as np
import csv
import re
import os


pur = ['N1', 'C2', 'N3', 'C4', 'C5', 'C6', 'N7', 'C8', 'N9']
pyr = ['N1', 'C2', 'N3', 'C4', 'C5', 'C6']


def get_coord(file, chain, base, base_id):
    """Extract the atomic coordinates of each base"""
    coords = {}
    with open(file, 'r') as f:
        for line in f:
            start = line.split()
            try:
                if start[-3] == chain and start[-5] == base_id and start[-4] == base:
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


def Base_plane(base, coords, cent_coor1, cent_coor2):
    """
    1. Determine purine or pyrimidine
    2. Calculate the base plane expression
    3. Calculate the coordinates of projection points
    4. Calculate the distance between bases, the Angle, and the distance between base pairs
    """
    A, B, C = 0, 0, 0
    L, Z, D, theta = 0, 0, 0, 0
    N1 = np.array([])
    try:
        if base == 'A' or base == 'G':
            N1 = coords['N1']
            C8 = coords['C8']
            N3 = coords['N3']
            A = N1[1] * (C8[2] - N3[2]) + C8[1] * (N3[2] - N1[2]) + N3[1] * (N1[2] - C8[2])
            B = N1[2] * (C8[0] - N3[0]) + C8[2] * (N3[0] - N1[0]) + N3[2] * (N1[0] - C8[0])
            C = N1[0] * (C8[1] - N3[1]) + C8[0] * (N3[1] - N1[1]) + N3[0] * (N1[1] - C8[1])
            D = N1[0] * (N3[1] * C8[2] - C8[1] * N3[2]) + \
                C8[0] * (N1[1] * N3[2] - N3[1] * N1[2]) + \
                N3[0] * (C8[1] * N1[2] - N1[1] * C8[2])
            # N = N1

        elif base == 'C' or base == 'U':
            N1 = coords['N1']
            N3 = coords['N3']
            C4 = coords['C4']
            A = N1[1] * (N3[2] - C4[2]) + N3[1] * (C4[2] - N1[2]) + C4[1] * (N1[2] - N3[2])
            B = N1[2] * (N3[0] - C4[0]) + N3[2] * (C4[0] - N1[0]) + C4[2] * (N1[0] - N3[0])
            C = N1[0] * (N3[1] - C4[1]) + N3[0] * (C4[1] - N1[1]) + C4[0] * (N1[1] - N3[1])
            D = N1[0] * (C4[1] * N3[2] - N3[1] * C4[2]) + \
                N3[0] * (N1[1] * C4[2] - C4[1] * N1[2]) + \
                C4[0] * (N3[1] * N1[2] - N1[1] * N3[2])
            # N = N3

        n = np.array([A, B, C])   # The normal vector
        t = (A * cent_coor2[0] + B * cent_coor2[1] + C * cent_coor2[2] + D) / (A * A + B * B + C * C)
        point_coord = np.array([cent_coor2[0] - A * t, cent_coor2[1] - B * t, cent_coor2[2] - C * t])  # Projection point coordinates
        Z = np.linalg.norm(cent_coor2 - point_coord)   # Vertical distance, coplanar or not
        D = np.linalg.norm(cent_coor1 - point_coord)  # The distance between the bases
        L = np.linalg.norm(cent_coor1 - cent_coor2)   # The distance between the centers of mass
        V1 = np.array([N1[0] - cent_coor1[0], N1[1] - cent_coor1[1], N1[2] - cent_coor1[2]])
        V2 = np.array([point_coord[0] - cent_coor1[0], point_coord[1] - cent_coor1[1], point_coord[2] - cent_coor1[2]])
        V1_mo = np.sqrt(V1.dot(V1))
        V2_mo = np.sqrt(V2.dot(V2))
        vector_product = V1.dot(V2)
        cos_ = vector_product / (V1_mo * V2_mo)
        angle_hu = np.arccos(cos_)
        angle = angle_hu * 180 / np.pi
        cross_product = np.cross(V1, V2)
        direction = cross_product.dot(n)
        print('direction：', direction)
        print('Angle：', angle)
        if direction >= 0:
            theta = angle
        else:
            theta = 2 * 180 - angle
        return round(L, 3), round(Z, 3), round(D, 3), round(theta, 3)

    finally:
        return round(L, 3), round(Z, 3), round(D, 3), round(theta, 3)


def cal_basepairing(root, csv_path, savepath):
    csv_file = open(csv_path, "r")
    reader = csv.reader(csv_file)
    for item in reader:
        pdb_id = item[0]
        base_1 = item[1]
        base_2 = item[2]
        bp, bp_name, LW = item[3], item[4], item[5]

        C1 = base_1.split('.')[0]  #
        C2 = base_2.split('.')[0]

        File = os.path.join(root, pdb_id + '.cif')  # root: dataset path

        try:
            id1 = re.findall('\d+$', base_1.split('.')[1])[0]
            id2 = re.findall('\d+$', base_2.split('.')[1])[0]

            b1 = base_1.split('.')[1].replace(id1, '')
            b2 = base_2.split('.')[1].replace(id2, '')

            coor1 = get_coord(File, C1, b1, id1)
            coor2 = get_coord(File, C2, b2, id2)

            centcoor1 = Centroid_coord(b1, coor1)
            centcoor2 = Centroid_coord(b2, coor2)

            b_1 = get_base_name(b1)
            b_2 = get_base_name(b2)

            l, Z, D, Theta = Base_plane(b_1, coor1, centcoor1, centcoor2)
            if l != 0 and Z != 0:
                data = [pdb_id, base_1, base_2, b_1, b_2, l, Z, D, Theta, bp, bp_name, LW]
                print(data)
                csvPath = os.path.join(savepath, 'base-pairing.csv')
                with open(csvPath, "a") as csvfile:
                    csv_writer = csv.writer(csvfile, lineterminator='\n')
                    csv_writer.writerow(data)
        finally:
            continue


def main():
    dataset_path = input('Enter the path of dataset-----: ')
    basepairs_save_path = input('Enter the path of statistics file of base-pairing---> :  ')
    save_path = input('Enter the path of results---> :  ')
    cal_basepairing(dataset_path, basepairs_save_path, save_path)
    choice = input('continue? Y/N  ').upper()
    print('-----------------')
    if choice == 'Y':
        main()
    if choice == 'N':
        quit()


if __name__ == '__main__':
    main()
