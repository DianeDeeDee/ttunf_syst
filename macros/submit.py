import string
import os

from optparse import OptionParser

print('Usage: python macros/submit.py [options]')

parser = OptionParser()
parser.add_option('-t', '--tag', dest='versionTag', help = 'String appended to the end of the output dataset', metavar='TAG', default='v0')
parser.add_option('-i', '--input', dest='inDSFile', help = 'Name of the text file that contains the name of the datasets on which this should run (one per line).',metavar='INPUT', default = 'datasets/ds_mc.txt') 
parser.add_option('-u', '--user', dest='user', help='AFS username used to build the ouput dataset name.', metavar='USER', default='dshoaleh')
parser.add_option('-x', '--suffix', dest='suffix', help='Suffix to be added in command line.', metavar='CMD', default='')
parser.add_option('-j', '--nFilesPerJob', dest='nFilesPerJob', help='Number of files per job', metavar='NJOBS', default='2')
parser.add_option('-s', '--site', dest='site', help='Send jobs to site.', metavar='SITE', default='')
parser.add_option('-S', '--save', dest='save', help='Save which output file?', metavar='OUTPUT', default='tree.root')
parser.add_option('-q', '--data', dest='isData', help='Is it data?', metavar='BOOL', default='False')
parser.add_option('-l', '--doLoose', dest='loose', help='Do QCD systematic variations? (0/1)', metavar='BOOL', default='False')
parser.add_option('-p', '--pdfName', dest='', help='the ttbar PDDF name?', metavar='BOOL', default='False')
parser.add_option('-n', '--nPdf', dest='', help='the number of PDF error sets?', metavar='BOOL', default='False')
parser.add_option('-b', '--destSE', dest='destSE', help='Send jobs to site.', metavar='SITE', default='')
(options, args) = parser.parse_args()

isData = (options.isData.capitalize() == 'True')

print "Options: ",options
print "Arguments: ", args

import sys
import random
import uuid
random.seed()

def getListOfFiles(infile, exts, badExt):
    l = []
    for i in infile:
        hasBad = False
        for ext in badExt:
            if ext in i:
                hasBad = True
                break
        if hasBad: 
            continue
        try:
            x = [i+'/'+k for k in os.listdir(i)]
            l.extend(getListOfFiles(x, exts, badExt))
        except:
            hasExt = False
            for ext in exts:
                if ext in i:
                    hasExt = True
                    break
            if hasExt:
                l.append(i)
    return l


def runJob(options):

  outPrefix = 'user.'+options.user+'.ttunf8TeV_'
  suffix = options.suffix
  extraTag = ''

  #simpletarFile = 'submit.tar.gz'
  #tarCmd2 = 'tar cfz '+simpletarFile+' down_grid.sh'
  #os.system(tarCmd2)
  simpletarFile = 'lastsubmit.tar.gz'

  cmd_line='grid_exec.sh --files input.txt '+suffix
 
#  strFileList = ''
#  first = True
#  for fi in fileList:
#    if fi != '':
#      if first:
#          strFileList = strFileList + fi
#          first = False
#      else:
#          strFileList = strFileList + ',' + fi

  file = open(options.inDSFile, 'r')

  singleOut = False
  tagSingleOut = ''
  if 'dataElectron' in options.inDSFile:
    tagSingleOut = 'DataElectron'
    singleOut = True
  elif 'dataMuon' in options.inDSFile:
    tagSingleOut = 'DataMuon'
    singleOut = True
  elif 'dataJet' in options.inDSFile:
    tagSingleOut = 'DataJet'
    singleOut = True

  print 'The command line will be:'
  print cmd_line

#  if singleOut:
#      outDS = outPrefix+tagSingleOut+extraTag+'.'+options.versionTag
#
#      inDSList = []
#      strInDSList = ''
#      first = True
#      for inDS in file.readlines():
#          inDS = inDS[0:-1]  # remove '\n'
#          print('Input  dataset: '+inDS)
#          inDSList.append(inDS)
#          if first:
#              first = False
#              strInDSList = strInDSList + inDS
#          else:
#              strInDSList = strInDSList + ',' + inDS
#
#      print('Output dataset: '+outDS)
#      siteStr = ''
#      if options.site != '':
#          siteStr = ' --site '+options.site
#
#      siteStr = siteStr + ' --excludedSite=ANALY_CONNECT '
#
#      strFilesPerJ = ' --nGBPerJob=MAX '
#      if int(options.nFilesPerJob) != -1:
#          strFilesPerJ = ' --nFilesPerJob '+(options.nFilesPerJob)
#      prun_line = 'prun --exec "echo %IN | sed \'s/,/\\n/g\' > input.txt ; cat input.txt ; '+cmd_line+'" --inDS '+strInDSList+' --outDS '+outDS+' --outputs tree.root --rootVer=5.34.19 --cmtConfig=x86_64-slc6-gcc47-opt '+strFilesPerJ+' --inTarBall '+simpletarFile+' --noCompile '+siteStr
#      print 'The prun command line will be: '
#      print prun_line
#      os.system(prun_line)

  if (isData):
   for inDS in file.readlines():
        inDS = inDS[0:-1]  # remove '\n'
        pos = string.find(inDS, 'data12_8TeV')
        #i = inDS[pos+12:pos+20] # the inDS name will be in the format 'mcXX_7TeV.ABCDEF.*' and i = ABCDEF
        i = inDS[pos+12:pos+19] #data12_8TeV.=12 caracters and data12_8TeV.periodJ = 19 caracters
        if 'a159' in inDS: # check for AtlFast tag
            i = i + 'a'
        outDS = outPrefix+i+extraTag+'.'+options.versionTag
        print('Input  dataset: '+inDS)
        print('Output dataset: '+outDS)

        siteStr = ''
        if options.site != '':
            #siteStr = ' --site '+options.site
             siteStr = ' --destSE '+options.site

        #siteStr = siteStr + ' --excludedSite=GRIF-IRFU_SCRATCHDISK CA-MCGILL-CLUMEQ-T2_LOCALGROUPDISK CA-MCGILL-CLUMEQ-T2_SCRATCHDISK'
        siteStr = siteStr + ' --excludedSite=GOEGRID_SCRATCHDISK ' #El_periodC
        #siteStr = siteStr + ' --excludedSite=ANALY_LIV_SL6 '
        strOut = options.save
        strFilesPerJ = ' --nGBPerJob=MAX '
        if int(options.nFilesPerJob) != -1:
             strFilesPerJ = ' --nFilesPerJob '+(options.nFilesPerJob)
        prun_line = 'prun --exec "echo %IN | sed \'s/,/\\n/g\' > input.txt ; cat input.txt ; '+cmd_line+'" --inDS '+inDS+' --outDS '+outDS+' --outputs tree.root --rootVer 5.34/19 --cmtConfig=x86_64-slc6-gcc47-opt '+strFilesPerJ+' --inTarBall '+simpletarFile+' --noCompile '+siteStr
        print 'The prun command line will be: '
        print prun_line
        os.system(prun_line)
  else:
   for inDS in file.readlines():
        inDS = inDS[0:-1]  # remove '\n'
        pos = string.find(inDS, 'mc12_8TeV')
        i = inDS[pos+10:pos+16] # mc12_8TeV.117050
        if 'a159' in inDS: # check for AtlFast tag
            i = i + 'a'
        outDS = outPrefix+i+extraTag+'.'+options.versionTag
        print('Input  dataset: '+inDS)
        print('Output dataset: '+outDS)

        siteStr = ''
        if options.site != '':
            #siteStr = ' --site '+options.site
            siteStr = ' --destSE '+options.site
            siteStr = siteStr + ' --excludedSite=CA-MCGILL-CLUMEQ-T2_LOCALGROUPDISK CA-MCGILL-CLUMEQ-T2_SCRATCHDISK'

        strOut = options.save
        strFilesPerJ = ' --nGBPerJob=MAX '
        if int(options.nFilesPerJob) != -1:
             strFilesPerJ = ' --nFilesPerJob '+(options.nFilesPerJob)
        prun_line = 'prun --exec "echo %IN | sed \'s/,/\\n/g\' > input.txt ; cat input.txt ; '+cmd_line+'" --inDS '+inDS+' --outDS '+outDS+' --outputs tree.root --rootVer 5.34.19 --cmtConfig=x86_64-slc6-gcc47-opt '+strFilesPerJ+' --inTarBall '+simpletarFile+' --noCompile '+siteStr
        print 'The prun command line will be: '
        print prun_line
        os.system(prun_line)
  file.close()      
#  else:
#      for inDS in file.readlines():
#        inDS = inDS[0:-1]  # remove '\n'
#        pos = string.find(inDS, 'mc12_8TeV')
#        i = inDS[pos+10:pos+16] # the inDS name will be in the format 'mcXX_7TeV.ABCDEF.*' and i = ABCDEF
#        if 'a159' in inDS: # check for AtlFast tag
#            i = i + 'a'
#        outDS = outPrefix+i+extraTag+'.'+options.versionTag
#        print('Input  dataset: '+inDS)
#        print('Output dataset: '+outDS)
#
#        siteStr = ''
#        if options.site != '':
#            siteStr = ' --site '+options.site
#
#        siteStr = siteStr + ' --excludedSite=ANALY_CONNECT '
#
#        strFilesPerJ = ' --nGBPerJob=MAX '
#        if int(options.nFilesPerJob) != -1:
#             strFilesPerJ = ' --nFilesPerJob '+(options.nFilesPerJob)
#        prun_line = 'prun --exec "echo %IN | sed \'s/,/\\n/g\' > input.txt ; cat input.txt ; '+cmd_line+'" --inDS '+inDS+' --outDS '+outDS+' --outputs tree.root --rootVer=5.34.19 --cmtConfig=x86_64-slc6-gcc47-opt '+strFilesPerJ+' --inTarBall '+simpletarFile+' --noCompile '+siteStr
#        print 'The prun command line will be: '
#        print prun_line
#        os.system(prun_line)
#  file.close()

runJob(options)

