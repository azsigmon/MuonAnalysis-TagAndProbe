import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

FITFUNC = "BWResCBExp"
#FITFUNC = "BWResCBCheb"

process.TnP_TrackSelection = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
    ## Input, output 
    InputFileNames = cms.vstring("file:TnP_Z_pp5TeV_data_trk_v4.root"),
    #InputFileNames = cms.vstring("file:TnP_Z_pp5TeV_MC_Zmu10mu10_trk_v4.root"),
    OutputFileName = cms.string(str("fits_pp5TeV_data_TrackSelection1_%s_v1.root" % FITFUNC)),
    #OutputFileName = cms.string(str("fits_pp5TeV_MC_Zmu10mu10_TrackSelection1_%s_v1.root" % FITFUNC)),
    InputTreeName = cms.string("fitter_tree"),
    InputDirectoryName = cms.string("tpTree"),
    # WeightVariable = cms.string("weight"),
    ## Variables for binning
    Variables = cms.PSet(
        mass   = cms.vstring("Tag-Probe Mass", "60", "120", "GeV/c^{2}"),
        pt     = cms.vstring("Probe p_{T}", "0", "1000", "GeV/c"),
        eta    = cms.vstring("Probe #eta", "-2.5", "2.5", ""),
        abseta = cms.vstring("Probe |#eta|", "0", "2.5", ""),
	tkDxyOverDxyError = cms.vstring("Track Dxy/DxyError", "0", "1000", ""),
	tkDzOverDzError = cms.vstring("Track Dz/DzError", "0", "10000", ""),
    ),
    ## Flags you want to use to define numerator and possibly denominator
    Categories = cms.PSet(
	Track_HP = cms.vstring("Track high purity", "dummy[pass=1,fail=0]"),
	Track_HISel_v1 = cms.vstring("Track selection v1", "dummy[pass=1,fail=0]"),
	Track_HISel_v2 = cms.vstring("Track selection v2", "dummy[pass=1,fail=0]"),
	tag_HIL2Mu15 = cms.vstring("HLT_HIL2Mu15", "dummy[pass=1,fail=0]"),
    ),
    Cuts = cms.PSet(
	dxy_cut = cms.vstring("dxy/dxyerror < 3", "tkDxyOverDxyError", "3.0"),
	dz_cut = cms.vstring("dz/dzerror < 3", "tkDzOverDzError", "3.0"),
    ),
    ## What to fit
    Efficiencies = cms.PSet(
        TrackSel_pt = cms.PSet(
            UnbinnedVariables = cms.vstring("mass"),
            EfficiencyCategoryAndState = cms.vstring("Track_HISel_v1", "pass", "dxy_cut", "below", "dz_cut", "below"), ## Numerator definition
            BinnedVariables = cms.PSet(
                ## Binning in continuous variables
                eta = cms.vdouble(-1.0, 1.0),
                pt = cms.vdouble( 10, 20, 30, 40, 50, 70, 100 ),
                ## flags and conditions required at the denominator, 
                tag_HIL2Mu15 = cms.vstring("pass"), ## i.e. use only events for which this flag is true
            ),
            BinToPDFmap = cms.vstring(FITFUNC), ## PDF to use, as defined below
        ),
        TrackSel_eta = cms.PSet(
            UnbinnedVariables = cms.vstring("mass"),
            EfficiencyCategoryAndState = cms.vstring("Track_HISel_v1", "pass", "dxy_cut", "below", "dz_cut", "below"),
            BinnedVariables = cms.PSet(
                eta = cms.vdouble(-2.4, -2.1, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2.1, 2.4),
                pt = cms.vdouble(20, 100),
                tag_HIL2Mu15 = cms.vstring("pass"),
            ),
            BinToPDFmap = cms.vstring(FITFUNC),
        ),
        TrackSel_1bin = cms.PSet(
            UnbinnedVariables = cms.vstring("mass"),
            EfficiencyCategoryAndState = cms.vstring("Track_HISel_v1", "pass", "dxy_cut", "below", "dz_cut", "below"),
            BinnedVariables = cms.PSet(
                eta = cms.vdouble(-1.0, 1.0),
                pt = cms.vdouble(20, 100),
                tag_HIL2Mu15 = cms.vstring("pass"),
            ),
            BinToPDFmap = cms.vstring(FITFUNC),
        ),
    ),

    ## PDF for signal and background (double voigtian + exponential background)
    PDFs = cms.PSet(
	VoigtExp = cms.vstring(
		"Voigtian::signal(mass, mean[91,85,95], width[3,1,10], sigma[3,1,10])",
		"Exponential::backgroundPass(mass, lp[0,-5,5])",
		"Exponential::backgroundFail(mass, lf[0,-5,5])",
		"efficiency[0.9,0,1]",
		"signalFractionInPassing[0.9]"
	),
	BWResCBExp = cms.vstring(
		"BreitWigner::bw(mass, m0[91.2,81.2,101.2], width[2.495,1,10])",
		"RooCBShape::res(mass, peak[0], sigma[1.7,0.01,10], alpha[1.8,0,3], n[0.8,0,10])",
		"FCONV::signal(mass, bw, res)",
		"Exponential::backgroundPass(mass, lp[0,-5,5])",
		"Exponential::backgroundFail(mass, lf[0,-5,5])",
		"efficiency[0.9,0.5,1]",
		"signalFractionInPassing[0.9]",
	),
	BWResCBCheb = cms.vstring(
		"BreitWigner::bw(mass, m0[91.2,81.2,101.2], width[2.495,1,10])",
		"RooCBShape::res(mass, peak[0], sigma[1.7,0.01,10], alpha[1.8,0,3], n[0.8,0,10])",
		"FCONV::signal(mass, bw, res)",
		"Chebychev::backgroundPass(mass, {c1p[0,-10,10], c2p[0,-10,10], c3p[0,-10,10]})",
		"Chebychev::backgroundFail(mass, {c1f[0,-10,10], c2f[0,-10,10], c3f[0,-10,10]})",
		"efficiency[0.9,0.5,1]",
		"signalFractionInPassing[0.9]",
	),
    ),

    ## How to do the fit
    binnedFit = cms.bool(True),
    binsForFit = cms.uint32(40),
    saveDistributionsPlot = cms.bool(True),
    NumCPU = cms.uint32(1), ## leave to 1 for now, RooFit gives funny results otherwise
    SaveWorkspace = cms.bool(True),
)

process.p = cms.Path(process.TnP_TrackSelection)
