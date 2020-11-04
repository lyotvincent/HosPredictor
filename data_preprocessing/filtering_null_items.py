
# 清除从ncbi virus下载的sequences.csv中host值为空的项
def filter(file_name):
    print('filter')
    f = open(file_name, 'r')
    lines = f.readlines()

    out = open('./sequences_filtered.csv', 'w')

    for line in lines:
        items = line.strip().split(',')
        if items[2] == '':
            continue
        out.write(line)


if __name__ == "__main__":
    filter('./sequences_waiting_for_filtering.csv')
