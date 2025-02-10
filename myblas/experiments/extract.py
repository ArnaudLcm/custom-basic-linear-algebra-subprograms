#!/bin/python3
import sys
from pathlib import Path

param_set = False

def create_csv(file):
    global param_set
    perfFirstLine = ""
    for l in file.readlines():
        if l == "\n":
            continue
        if l[:5]=="Test:":
            continue
        if l[1]=='(':
            continue

        perfCurLine = ""
        for w in l.split()[:-1]:
            if w == ':':
                continue
            s = w.split('=')
            
            if len(s) > 1:
                if not param_set:
                    perfFirstLine = perfFirstLine + s[0] +','
                if s[-1]=='':
                    continue
                w = s[-1]
            
            if w[-1] == ':':
                w = w[:-1]

            perfCurLine = perfCurLine + w + ','
        
        if not param_set:
            perfFirstLine = perfFirstLine + 'GFlops\n'
            perfFile.write(perfFirstLine)
            param_set = True
        perfCurLine = perfCurLine[:-1]
        perfCurLine = perfCurLine + "\n"
        perfFile.write(perfCurLine)


nbargs = len(sys.argv)
if  nbargs < 2:
    print("./extract.py <input> <output>\n")
    exit(1)

logPath = Path(sys.argv[1])

if nbargs == 2:
    perfPath = Path.cwd().joinpath(logPath.stem + ".csv")
else:
    perfPath = Path(sys.argv[2])
    if perfPath.is_dir():
        perfPath = perfPath.joinpath(logPath.stem + ".csv")

try:
    perfFile = perfPath.open("a")
except:
    print("Could not open/create csv file")
    exit(1)



logList = ""

if logPath.is_dir():
    logList = logPath.iterdir()
else:
    logList = [logPath]

for p in logList:
    if p.is_file():
            try:
                logFile = p.open("r")
            except:
                print("Could not open %s".format(p))
                exit(1)
            create_csv(logFile)
            logFile.close()

perfFile.close()
