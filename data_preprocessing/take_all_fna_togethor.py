import glob
import subprocess


if __name__ == "__main__":
    file_names = glob.glob('./all.fna/*/*.fna')

    for f in file_names:
        print(f)
        subprocess.run('mv %s ./all_fnas/%s' % (f, f.split('/')[-1]), shell=True, check=True)
