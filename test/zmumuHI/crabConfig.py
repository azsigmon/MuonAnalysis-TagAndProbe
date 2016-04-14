from WMCore.Configuration import Configuration

config = Configuration()

config.section_("General")
#config.General.requestName   = 'tnp_pp_data_v1'
#config.General.requestName   = 'tnp_pp_data_trk_v4'
#config.General.requestName   = 'tnp_PbPb_data_v1'
#config.General.requestName   = 'tnp_PbPb_data_trk_v4'
#config.General.requestName   = 'tnp_pp_mc_v3'
#config.General.requestName   = 'tnp_pp_mc_trk_v4'
#config.General.requestName   = 'tnp_PbPb_mc_v2'
config.General.requestName   = 'tnp_PbPb_mc_trk_v4'

config.section_("JobType")
config.JobType.pluginName  = 'Analysis'
#config.JobType.psetName    = 'tnp_pp_data.py'
#config.JobType.psetName    = 'tnp_pp_data_tracking.py'
#config.JobType.psetName    = 'tnp_PbPb_data.py'
#config.JobType.psetName    = 'tnp_PbPb_data_tracking.py'
#config.JobType.psetName    = 'tnp_pp_mc.py'
#config.JobType.psetName    = 'tnp_pp_mc_tracking.py'
#config.JobType.psetName    = 'tnp_PbPb_mc.py'
config.JobType.psetName    = 'tnp_PbPb_mc_tracking.py'

config.section_("Data")
#config.Data.inputDataset = '/SingleMuHighPt/Run2015E-PromptReco-v1/AOD'
#config.Data.inputDataset = '/HIEWQExo/HIRun2015-PromptReco-v1/AOD'
#config.Data.inputDataset = '/Pythia8_Zmu10mu10/HINppWinter16DR-75X_mcRun2_asymptotic_ppAt5TeV_v3-v1/AODSIM'
config.Data.inputDataset = '/Pythia8_Zmu10mu10_Hydjet_MB/HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1/AODSIM'
#config.Data.splitting = 'LumiBased'
config.Data.splitting = 'FileBased'
#config.Data.unitsPerJob = 30
config.Data.unitsPerJob = 1
#config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/5TeV/Cert_262081-262328_5TeV_PromptReco_Collisions15_25ns_JSON.txt'
#config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/HI/Cert_262548-263757_PromptReco_HICollisions15_JSON.txt'
config.Data.publication = False
config.Data.outLFNDirBase = '/store/group/phys_heavyions/azsigmon/TagAndProbe2015/'

config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'
