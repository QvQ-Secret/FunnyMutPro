# -*- coding: UTF-8 -*-
import os
import platform
import re
import subprocess


def get_usalign_path(exe="USalign"):
    if platform.system().lower().startswith("win"):
        exe += ".exe"
    filename = os.path.join(os.path.dirname(os.path.abspath(__file__)), exe)
    if os.path.isfile(filename):
        return filename
    else:
        for p in os.getenv("PATH").split(os.pathsep):
            filename = os.path.join(p, exe)
            if os.path.isfile(filename):
                return filename
    print("ERROR! Cannot locate %s at %s or at %s" % (exe,
                                                      os.path.dirname(os.path.abspath(__file__)), os.getenv("PATH")))
    print("Please put the USalign executable at one of the aforementioned paths")
    return exe


def us_align(mobile_filename, target_filename, args='', exe=''):
    if len(exe) == 0:
        exe = get_usalign_path("USalign")
    if args == '""':
        args = ''
    if len(args) > 2 and args[0] == '"' and args[-1] == '"':
        args = args[1:-1]
    if "-outfmt" not in args:
        args += " -outfmt -1"
    args = ' '.join([exe, mobile_filename, target_filename, args, '-m -'])
    process = subprocess.Popen(args, stdout=subprocess.PIPE, shell=True,
                               universal_newlines=True)
    lines = process.stdout.readlines()
    alignLength, rmsd, ratio, tmScoreQ, d0Q, tmScoreM, d0M = '/', '/', '/', '/', '/', '/', '/'
    for line in iter(lines):
        line = line.rstrip()
        if 'length=' in line:
            alignLength = re.search(r'Aligned length=(.*?),', line).group(1).strip()
            rmsd = re.search(r'RMSD=(.*?),', line).group(1).strip()
            ratio = re.search(r'Seq_ID=n_identical/n_aligned=(.*)', line).group(1).strip()
        if 'TM-score=' in line:
            if 'Structure_1' in line:
                tmScoreQ = re.search(r'TM-score=(.*?)\(', line).group(1).strip()
                d0Q = re.search(r'd0=(.*?)\)', line).group(1).strip()
            if 'Structure_2' in line:
                tmScoreM = re.search(r'TM-score=(.*?)\(', line).group(1).strip()
                d0M = re.search(r'd0=(.*?)\)', line).group(1).strip()
    returnDict = {
        'alignLength': alignLength,
        'rmsd': rmsd,
        'n_identical/n_aligned': ratio,
        'tmScoreQ': tmScoreQ,
        'd0Q': d0Q,
        'tmScoreM': tmScoreM,
        'd0M': d0M,
    }
    return returnDict


us_align(r'G:\RCSB_PDB\PDB_FILES\package_023\1ENH.pdb', r'G:\RCSB_PDB\PDB_FILES\package_026\1ZTR.pdb')
