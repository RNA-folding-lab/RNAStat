import json
import os
import csv


def motif(folderpath, file, name):
    with open(file) as json_data:
        dic = json.load(json_data)

    # 茎
    num_stems = 0
    nts_stems = []
    if 'num_stems' in dic:
        num_stems = dic['num_stems']
        for a in dic['stems']:
            nts_stems_num = int(a['num_pairs'])
            nts_stems.append(nts_stems_num)

    # 发夹环
    num_hairpins = 0
    hairpins_num_nts = []
    if 'num_hairpins' in dic:
        num_hairpins = dic['num_hairpins']
        for b in dic['hairpins']:
            nts_hairpins_num = int(b['num_nts']) - 2
            hairpins_num_nts.append(nts_hairpins_num)

    # 突环
    num_bulges = 0
    bulges_num_nts = []
    if 'num_bulges' in dic:
        num_bulges = dic['num_bulges']
        for c in dic['bulges']:
            nts_bulges_num = int(c['num_nts']) - 4
            bulges_num_nts.append(nts_bulges_num)

    # 内部环
    num_iloops = 0
    iloops_num_nts = []
    if 'num_iloops' in dic:
        num_iloops = dic['num_iloops']
        for d in dic['iloops']:
            nts_iloops_num = int(d['num_nts']) - 4
            iloops_num_nts.append(nts_iloops_num)

    # 多分支环
    num_junctions = 0
    junctions_num_nts = []
    if 'num_junctions' in dic:
        num_junctions = dic['num_junctions']
        for e in dic['junctions']:
            stems_num = int(e['num_stems'])
            nts_junctions_num = int(e['num_nts']) - stems_num * 2
            junctions_num_nts.append(nts_junctions_num)

    data = [name[0:4], num_stems, nts_stems, num_hairpins, hairpins_num_nts, num_bulges, bulges_num_nts, num_iloops,
            iloops_num_nts, num_junctions, junctions_num_nts]
    file_path = os.path.join(folderpath, 'motif.csv')
    with open(file_path, "a") as csvfile:
        # 参数a，追加数据
        csv_writer = csv.writer(csvfile, lineterminator='\n')
        csv_writer.writerow(data)


def base_pairs(folderpath, file, name):
    with open(file) as json_data:
        dic = json.load(json_data)
    if 'num_pairs' in dic.keys():
        num = dic['num_pairs']
        for i in dic['pairs']:
            nt1 = i['nt1']
            nt2 = i['nt2']
            bp = (i['bp'][0] + i['bp'][2]).upper()
            bp_name = i['name']
            LW = i['LW']
            data = [name[0:4], nt1, nt2, bp, bp_name, LW]
            print(data)
            csvPath = os.path.join(folderpath, 'basepairs.csv')
            with open(csvPath, "a") as csvfile:
                # 参数a，追加数据
                csv_writer = csv.writer(csvfile, lineterminator='\n')
                csv_writer.writerow(data)


def Stat():
    json_root = input("Enter the json file path, such as ' D:\dataset '------ : ")
    path_root = input("Enter the  save_file path, such as ' D:\dataset '------ : ")
    file_list = os.listdir(json_root)
    for i in file_list:
        print(i)
        FilePath = os.path.join(json_root, i)
        base_pairs(path_root, FilePath, i)
        motif(path_root, FilePath, i)


def main():
    Stat()
    choice = input('continue? Y/N').upper()
    print('-----------------')
    if choice == 'Y':
        Stat()
    if choice == 'N':
        quit()


if __name__ == '__main__':
    main()