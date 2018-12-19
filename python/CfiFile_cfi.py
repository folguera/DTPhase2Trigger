import FWCore.ParameterSet.Config as cms

demo = cms.EDAnalyzer('DTPhase2Trigger'
     ,tracks = cms.untracked.InputTag('ctfWithMaterialTracks')
)
