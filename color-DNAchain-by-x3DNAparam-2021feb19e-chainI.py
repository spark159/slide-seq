import sys

'''
This script should generate a pymol script that colors DNA residues based on a 2-column text file, where
the second column is a float that represents a basestep parameter (or difference from a reference).

This should be run by python3

Inputs are
- molecule name in pymol
- chain ID
- start residue (assumption is that the parameters for the first residue are labeled 1)
- min and max values of the parameter of interest
- color to use; for now, going to just use RGB varying from pure red [1,0,0] to white [1,1,1]
    values below min are [1,1,1] and values more than max are [1,0,0]

'''
### set variables
molname="chd1apo"
chainID="I"
resiOFFSET = 0 ### this is to adjust second DNA chain so color matches with first; here we will make 0 for TApoor and 1 for TArich
parametervalueMAX=20
parametervalueMIN=-20
parametervalueMID=0
colorMAX=[1,0,0] ### red
colorMID=[1,1,1] ### white
colorMIN=[0,0,1] ### blue
colorBLACK=[0,0,0]
scriptversion="3"

### read in 2-column file given on commandline 
try:
    parameterfileNAME = sys.argv[1]
    file = open(parameterfileNAME,'r')
    parameterfileINPUT=file.readlines()
    file.close()
except IOError:  
    print("\nCannot read ", sys.argv[1])
    print("please supply a PDB file.\n")
    sys.exit("EXITING...")
except IndexError:
    print("need to supply pdb filename")
    sys.exit("EXITING...")

### for each line of the inputfile, write a line of pymol script
pymolscriptname="pymolscriptoutput-"+molname+"-chain"+chainID+"-v"+str(scriptversion)+".pml"
pymoloutput = open(pymolscriptname,"w")

### header for pymol script
pymoloutput.write("### this script should color DNA residues according to parameters given in file\n")
pymoloutput.write("# {}\n\n".format(parameterfileNAME))
pymoloutput.write("### this script was run with the following parameters: \n")
pymoloutput.write("# molname = {}\n".format(molname))
pymoloutput.write("# chainID = {}\n".format(chainID))
pymoloutput.write("# resiOFFSET = {}\n".format(resiOFFSET))
pymoloutput.write("# parametervalueMAX = {}\n".format(parametervalueMAX))
pymoloutput.write("# parametervalueMIN = {}\n".format(parametervalueMIN))
pymoloutput.write("# parametervalueMID = {}\n".format(parametervalueMID))
pymoloutput.write("# scriptversion = {}\n".format(scriptversion))

### for each line in file, where first and second item are numbers, set the resi to first
#    and use second to generate the color
#    if the second parameter is above MAX, then make color MAX
#    if the second parameter is bewteen MID and MAX, then divide param by MAX and use that fraction
#    to scale all values of colorMAX
#    if the second parameter is below MIN, then make color MIN
#    if the second parameter is between MIN and MID, then divide param by min and use that fraction
#
#    in the case that MID is not simply zero, should get the ranges of MAX-MID and MID-MIN for the 
#    division steps 

for line in parameterfileINPUT:
    try:  ### only using lines that have two numbers; skipping those with words
        inputDATA = line.split()
        resiCURRENT=int(inputDATA[0]) + resiOFFSET
        valueRAW=float(inputDATA[1])
        print("resi = {}; value = {}".format(resiCURRENT, valueRAW),end="")
        ### here scaling the input parameter (second column) relative to MAX, MID, MIN values
        #   and creating 3-number list for rgb color
        if valueRAW > parametervalueMAX: 
            tempCOLOR = colorMAX
        elif valueRAW < -1000: ### this is the unique case where bpstep is missing, and set parameters to 2000 for reference
            tempCOLOR = colorBLACK
        elif valueRAW < parametervalueMIN:
            tempCOLOR = colorMIN
        elif valueRAW == parametervalueMID:
            tempCOLOR = colorMID
        elif valueRAW > parametervalueMID: # Between MID and MAX
            colorFRACTION = (parametervalueMAX - valueRAW)/(parametervalueMAX - parametervalueMID)
            ### here, only want to change rgb values that are otherwise zero, so making mask of [0,1,1] for red
            tempCOLOR = [colorMAX[0]+(colorMID[0]-colorMAX[0])*colorFRACTION,colorMAX[1]+(colorMID[1]-colorMAX[1])*colorFRACTION,colorMAX[2]+(colorMID[2]-colorMAX[2])*colorFRACTION]
        elif valueRAW < parametervalueMID: #between MID and MIN
            colorFRACTION = (valueRAW - parametervalueMIN)/(parametervalueMID - parametervalueMIN)
            tempCOLOR = [colorMIN[0]+(1-colorMIN[0])*colorFRACTION,colorMIN[1]+(colorMID[1]-colorMIN[1])*colorFRACTION,colorMIN[2]+(colorMID[2]-colorMIN[2])*colorFRACTION]
        else: ### this should never happen
            tempCOLOR = [0,0,0]
        print("tempCOLOR = {}".format(tempCOLOR))
        print("")
        ### here need to have a unique color for each residue, chain, mol
        pymoloutput.write("set_color tempcolor{}{}{},[{},{},{}]\n".format(molname, chainID,resiCURRENT,tempCOLOR[0],tempCOLOR[1],tempCOLOR[2]))
        pymoloutput.write("color tempcolor{}{}{},{} and chain {} and resi {}\n".format(molname, chainID,resiCURRENT,molname, chainID,resiCURRENT))
    except:
        print("skipping line {}".format(inputDATA))
        continue






pymoloutput.close()
