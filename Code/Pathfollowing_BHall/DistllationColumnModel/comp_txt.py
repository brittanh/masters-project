#!/opt/local/bin/python
# -*- encoding: ascii -*-
"""
    @purpose: Compare two text files
    @author: Brittany Hall
    @date: 30.10.2017
    @version: 0.1
    @updates:
"""
File1 = 'lbw.txt'
File2 = 'lbw_py.txt'

#def compare(File1,File2):
#    with open(File1,'r') as f:
#        d=set(f.readlines())
#
#
#    with open(File2,'r') as f:
#        e=set(f.readlines())
#
#    open('file3.txt','w').close() #Create the file
#
#    with open('differences_lbw.txt','a') as f:
#        for line in list(d-e):
#            f.write(line)
#
#compare(File1,File2)
def compare(File1,File2):
    f1 = open(File1, "r")
    f2 = open(File2, "r")

    fileOne = f1.readlines()
    fileTwo = f2.readlines()

    f1.close()
    f2.close()
    outFile = open("file_lbw.txt", "w")
    x = 0
    for i in fileOne:
        if i != fileTwo[x]:
            outFile.write(i+" <> "+ fileTwo[x])
        x += 1
    outFile.close()

compare(File1,File2)
