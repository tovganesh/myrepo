#
# extractConvergedGeometry.py - extract final geometry from output file
# currently supports only GAMESS
# 
# @author V. Ganesh
# @date 2 Aug 2017
# @license GPL v3
#
# (c) V. Ganesh
#

import sys
import string

# extractor for GAMESS
def extractGAMESSGeometry(fileName):
   fl = open(fileName)
   lines = map(string.strip, fl.readlines())
   fl.close()

   idx = 0
   str = ""
   for line in lines:
      if (line.find("***** EQUILIBRIUM GEOMETRY LOCATED *****") >= 0):
         for jdx in range(idx+4, len(lines)):
           if (lines[jdx] == '' or lines[jdx] == None): break
           str += lines[jdx] + '\n'

      idx += 1   

   return str

# mapping
extractors = {}
extractors["GAMESS"] = extractGAMESSGeometry

# driver
def extractConvergedGeometry(fileName, program="GAMESS"):
   global extractors 

   print(extractors[program](fileName))   

if __name__=="__main__":
   if (len(sys.argv) != 2):
     print("usage: python " + sys.argv[0] + " <output file>")
     sys.exit(10)

   extractConvergedGeometry(sys.argv[1])
