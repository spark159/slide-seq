#!/usr/bin/env python
import sys
'''
This script expects a text file with pairs of rows, first is name, second has data.
The data row for each should have the same number of items, separated by commands

This script will convert these to columns, separated by spaces.

usage:
python thisscript.py textfile.txt
'''

rawfile = open(sys.argv[1],'r')
filelist = rawfile.readlines()
rawfile.close()

mylists = []
for i in range(len(filelist)):
    if filelist[i][0]==">":
        continue ### go to next line
    mylists.append([]) ### record only rows with data

listtransposedyet=0

for i in range(len(filelist)):
    if filelist[i][0]==">":
        continue ### go to next line

    ### create listtransposed if empty
    if listtransposedyet==0:
        listtransposedyet=1
        listtransposed=[]
        templist = filelist[i].split(",")
        for j in range(len(templist)):
            listtransposed.append([])
    
    ### listtransposed should now exist
    templist = filelist[i].split(",")
    #print(templist)
    for j in range(len(templist)):
        listtransposed[j].append(templist[j].strip())

for mylist in listtransposed:
    print(mylist)
outputfilename = sys.argv[1]+"_transposed"
outfile = open(outputfilename,'w')

for j in range(len(listtransposed)):
    for i in range(len(listtransposed[j])):
        outfile.write("{} ".format(listtransposed[j][i]))
    outfile.write("\n")

outfile.close()

print("")
print("written file: {}".format(outputfilename))
print("")










