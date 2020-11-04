import glob

# 清除从ncbi virus下载的sequences.csv中存在，但是all.fna.tar.gz里不存在的项
def filter(all_fnas):
    print('filter')
    f = open('sequences_filtered.csv', 'r')
    lines = f.readlines()

    out = open('./sequences_only_all_fna.csv', 'w')

    for line in lines:
        items = line.strip().split(',')
        if items[0] not in all_fnas:
            continue
        else:
            del all_fnas[all_fnas.index(items[0])]
            out.write(line)
    print(all_fnas)


if __name__ == "__main__":
    file_names = glob.glob('./all_fnas/*.fna')
    all_fnas = []
    for f in file_names:
        all_fnas.append(f[11:-4])
    # print(all_fnas)
    filter(all_fnas)

