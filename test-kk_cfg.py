

import FWCore.ParameterSet.Config as cms

########################################
# Command line argument parsing
########################################
import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing ('analysis')
options.register ('crab',
                  0, # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.int,          # string, int, or float
                  "Set to 1 to run on CRAB.")
steps_options = ['skim', 'analyze'] # Valid options for steps 
options.register ('steps',
                  [],
                  VarParsing.VarParsing.multiplicity.list, # singleton or list
                  VarParsing.VarParsing.varType.string,          # string, int, or float
                  "Steps to execute. Possible values: skim, analyze.")
options.steps = ['skim', 'analyze'] # default value
# Get and parse the command line arguments
options.parseArguments()
# Check validity of command line arguments
for step in options.steps:
    if step not in steps_options:
        print "Skipping invalid steps: %s" % step
        options.steps.remove(step)

print options.steps

from EmergingJetGenAnalysis.Configuration.GenTools import *

process = cms.Process('QQQ')
if 'skim' in options.steps and len(options.steps)==1:
    # If only running skim, add AOD/AODSIM and jetFilter/wJetFilter to output
    process.setName('SKIM')

########################################
# Stable configuration
########################################
# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')


## Geometry and Detector Conditions (needed for a few patTuple production steps)
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

## Options and Output Report
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False),
        # SkipEvent = cms.untracked.vstring('ProductNotFound')
)

# Unscheduled execution
# process.options.allowUnscheduled = cms.untracked.bool(False)
process.options.allowUnscheduled = cms.untracked.bool(True)

########################################
# Skim
########################################
import os
cmssw_version = os.environ['CMSSW_VERSION']
skimStep = cms.Sequence()
if 'skim' in options.steps:
    print ''
    print '####################'
    print 'Adding Skim step'
    print '####################'
    print ''
#    skimStep = addSkim(process)
########################################
# Analyze
########################################
analyzeStep = cms.Sequence()
if 'analyze' in options.steps:
    print ''
    print '####################'
    print 'Adding Analyze step'
    print '####################'
    print ''
    analyzeStep = addAnalyze(process)


process.p = cms.Path( skimStep * analyzeStep )
process.out = cms.OutputModule("PoolOutputModule",
        fileName = cms.untracked.string('output.root'),
        outputCommands = cms.untracked.vstring('drop *'),
    )


########################################
# Generic configuration
########################################
if 'CMSSW_7_4_12' in cmssw_version:
    globalTags=['74X_mcRun2_design_v2']
elif 'CMSSW_7_4_1_patch4' in cmssw_version:
    globalTags=['MCRUN2_74_V9']
elif 'CMSSW_7_6_3' in cmssw_version:
    globalTags=['76X_mcRun2_asymptotic_RunIIFall15DR76_v1']
print 'CMSSW_VERSION is %s' % cmssw_version
print 'Using the following global tags [MC, DATA]:'
print globalTags
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, globalTags[0], '')

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1)
process.MessageLogger.cerr.FwkReport.limit = 20
process.MessageLogger.cerr.default.limit = 100

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20) )

process.source = cms.Source("PoolSource",
    # eventsToProcess = cms.untracked.VEventRange("1:36:3523-1:36:3523"),
    fileNames = cms.untracked.vstring('file:/data/users/hon-martyn/CMSSW_7_6_3/src/GEN1KEventsAt6TeV.root'
)
)
process.TFileService = cms.Service("TFileService", fileName = cms.string("ntuple.root") )

# storage
process.outpath = cms.EndPath(process.out)
