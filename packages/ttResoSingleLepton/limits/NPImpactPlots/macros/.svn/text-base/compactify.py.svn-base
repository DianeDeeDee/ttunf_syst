import glob, os, sys, re, operator
from math import sqrt
from optparse import OptionParser
"""
  Need as input a card name + fit result
"""

###############################################
# Edit this block, the convention is
# pattern in the systematics : new name
# 
# the systematics to be grouped have to *start* with the pattern
#

group_name        = {
  "BTAGBREAK":"b-tagging efficiency",
  "CTAGBREAK":"c-tagging efficiency",
  "LTAGBREAK":"Light jet-tagging efficiency",
  "^Jet":"Jet energy scale",
  "ttbar-DataRw":"tt reweighting systematics",
  "QCD":"Multijet normalisation",
  "MUON|ELE":"Lepton systematics",
  "Wjets|Zjets":"V+jets normalisation",
  "ttH":"ttH modeling",
  "MG":"tt modeling",
  "PartonShower":"Parton shower",
  "XS":"Theoretical cross sections",
  "bb|cc":"tt heavy-flavour fractions",
  "JEFF" : "Jet reconstructin efficiency",
  "JER" : "Jet energy resolution",
  "JVF" : "Jet vertex fraction",
  "LUMI" : "Luminosity",
  "SigXsecOverSM" : "Signal strength",
  }

# renaming but not grouping some systematics
rename = {
  }

#
##############################################


def main():
  parser = OptionParser() 
  parser.add_option("-f", "--fit", 
                    dest="fit", 
                    help="Fit result to be used",
                    default="")
  parser.add_option("-i", "--in", 
                    dest="incard", 
                    help="In card name",
                    default="")
  parser.add_option("-o", "--out", 
                    dest="outcard", 
                    help="Out card name",
                    default="")
  options,args = parser.parse_args() 
  if not options.fit or not options.incard or not options.outcard:
    print "Usage: python compactify.py -i incardname -o outcardname -f fitresult"
    sys.exit(1)
  
  correlationMatrix = readCorrMatrix(options.fit)
  ascii_input = readAscii(options.incard)

  group_sampleerror = {}
  group_sys         = {}
  group_sys2         = {}
  for key in group_name.keys():
    group_sampleerror[key] = []
    group_sys[key] = []
    group_sys2[key] = []

  start = False
  written = False
  for sysname, errors in ascii_input.iteritems():
    for key,g_errors in group_sampleerror.iteritems():
      if re.search(key,sysname):
        if len(g_errors)==0:
          g_errors = [errors[0]]
          for err in errors[1:]:
            g_errors.append(err*err)
          group_sampleerror[key] = g_errors
        else:
          assert len(errors) == len(g_errors),(line,errors,group_sampleerror)
          for i,gerr in enumerate(g_errors):
            if i==0: continue
            val = gerr + errors[i]* errors[i]
            for corrsys, correrrors in group_sys[key]:
              val += 2*correlationMatrix[(corrsys,sysname)]*errors[i]*correrrors[i]
            g_errors[i] = val
        group_sys[key].append( (sysname,errors) )
        group_sys2[key].append( sysname )
        break
    else:
      print "non recognized systematic",sysname
      print group_name,group_sampleerror
      sys.exit(1)
  print group_sampleerror
  writeAscii(options.incard,options.outcard,group_sampleerror, group_name)

def readCorrMatrix(fitfile):

  sys = []
  corrMatrix = {}
  with open(fitfile) as f:
    read = False
    i = 0
    for line in f:
      if "alpha" in line or "SigX" in line: 
        name = line.split()[0].replace("alpha_","")
        sys.append(name)
        continue
      if "CORRELATION_MATRIX" in line: read = True
      if not read: continue
      tokens = line.split()
      if len(tokens) < 3: continue
      for j, corr in enumerate(tokens):
        corrMatrix[(sys[i],sys[j])] = float(corr)
      i += 1
  return corrMatrix
    
def readAscii(incard):

  names = []
  with open("../ascii/"+incard+"_pulls_id.txt") as f:
    for line in f:
      names.append(line.replace("alpha_","").replace("\n",""))
  errors = []
  with open("../ascii/"+incard+"_pulls.txt") as f:
    for line in f:
      id, pull, errup, errdown, mu, postup, postdown, preup, predown = [float(x) for x in line.split()]
      error = (mu, mu-postup, mu-postdown, mu-preup, mu-predown )
      errors.append(error)

  return dict(zip(names,errors))

def  writeAscii(incard, outcard, group_sampleerror, group_name):

  todel = []
  for key,errors in group_sampleerror.iteritems():
    if not errors:
      todel.append(key)
  for key in todel: del group_sampleerror[key]
  for key in todel: del group_name[key]

  with open("../ascii/"+outcard+"_pulls_id.txt","w") as f:
    for key in group_sampleerror.iterkeys():
      value = group_name[key].replace(" ","_")
      f.write(value+"\n")
  with open("../ascii/"+outcard+"_pulls.txt","w") as f:
    for i,(key,errors) in enumerate(group_sampleerror.iteritems()):
      mu, postup, postdown, preup, predown = errors
      line = "%d 0. 1. 1. %f %f %f %f %f" % (i,mu,mu-postup,mu+postdown,mu-preup,mu-predown)
      f.write(line+"\n")

  for i in ("breakdown_add","breakdown_add_nf","pulls_nf"):
    filein = "../ascii/"+incard+"_"+i+".txt"
    fileout = "../ascii/"+outcard+"_"+i+".txt"
    os.system("cp "+filein+" "+fileout)
    filein = filein.replace(".txt","_id.txt")
    fileout = fileout.replace(".txt","_id.txt")
    os.system("cp "+filein+" "+fileout)

if __name__ == "__main__":
  main()
