
def add_class():
    hosts_class = open('distinct_hosts.csv', 'r')
    class_lines = hosts_class.readlines()
    hosts_class.close()

    class_dict = dict()

    for l in class_lines:
        l = l.strip().split(',')
        class_dict[l[0]] = l[2]

    only_all_fna = open('sequences_only_all_fna.csv', 'r')
    only_lines = only_all_fna.readlines()
    only_all_fna.close()

    out = open('sequences_only_all_fna_with_class.csv', 'w')

    for l in only_lines:
        if l == 'Accession,Length,Host\n':
            out.write('Accession,Length,Host,Class\n')
        else:
            if class_dict[l.strip().split(',')[2]] == 'others':
                out.write(l.strip()+',others\n')
            else:
                out.write(l.strip()+',Mammalia|Aves\n')

    out.close()

if __name__ == "__main__":
    add_class()
