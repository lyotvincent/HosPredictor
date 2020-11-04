import sys


def count_gc(k):
    print('filter')
    f = open('./sequences_only_all_fna_with_class.csv', 'r')
    lines = f.readlines()

    out = open('./sequences_with_GC_%s.csv' % k, 'w')

    title = 'Accession,Length,Host,Class,GC_content'
    bases = ['A', 'T', 'C', 'G']

    k = int(k)
    for num in range(4**k):
        rest = num
        field = ''
        for i in range(k):
            divisor = 4**(k-1-i)
            if rest < divisor:
                root = 0
            else:
                root = rest // divisor
                rest -= root * divisor
            field += bases[root]
        title += ',' + field +'_num'

    title += ',Sequence\n'

    for line in lines:
        if line == 'Accession,Length,Host,Class\n':
            out.write(title)
            continue
        items = line.strip().split(',')
        fna = open('./all_fnas/%s.fna' % items[0], 'r')
        fna_lines = fna.readlines()
        fna_seq = ''
        for fna_line in fna_lines:
            if fna_line.startswith('>'):
                continue
            fna_seq += fna_line.strip()
        fna_seq.upper()
        gc_content = str( round( (fna_seq.count('G')+fna_seq.count('C'))/len(fna_seq), 2) )
        result_line = line.strip()+','+gc_content
        for num in range(4**k):
            rest = num
            field = ''
            for i in range(k):
                divisor = 4**(k-1-i)
                if rest < divisor:
                    root = 0
                else:
                    root = rest // divisor
                    rest -= root * divisor
                field += bases[root]
            result_line += ',' + str(fna_seq.count(field))
        result_line += ',' + fna_seq+'\n'
        out.write(result_line)


if __name__ == "__main__":
    k=sys.argv[1]
    print(k)
    count_gc(k)
