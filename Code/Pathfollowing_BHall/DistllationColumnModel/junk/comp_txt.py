#!/opt/local/bin/python
# -*- encoding: ascii -*-
"""
    @purpose: Compare two text files
    @author: Brittany Hall
    @date: 30.10.2017
    @version: 0.1
    @updates:
"""
File1 = 'dMxdt_matlab.txt'
File2 = 'dMxdt_py.txt'
filename = 'file_dMxdt.txt'

def compare(File1,File2, filename):
    f1 = open(File1, "r")
    f2 = open(File2, "r")

    fileOne = f1.readlines()
    fileTwo = f2.readlines()

    f1.close()
    f2.close()
    outFile = open(filename, "w")
    x = 0
    for i in fileOne:
#        print i
#        print fileTwo[x]
#        raw_input()
        if i != fileTwo[x]:
#if float("{:.8f}".format(float(i))) != float("{:.8f}".format(float(fileTwo[x]))):
            outFile.write(i+" <> "+ fileTwo[x])
        x += 1
    outFile.close()

compare(File1,File2,filename)
