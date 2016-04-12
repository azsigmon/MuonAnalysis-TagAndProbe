// system include files
#include <memory>
#include <cmath>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include <DataFormats/MuonReco/interface/Muon.h>

#include "DataFormats/CLHEP/interface/AlgebraicObjects.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"


//
// class declaration
//

class MuonImpactParameter : public edm::EDProducer {
public:
  explicit MuonImpactParameter(const edm::ParameterSet&);
  ~MuonImpactParameter();

private:
  virtual void produce(edm::Event&, const edm::EventSetup&);

      const edm::EDGetTokenT<edm::View<reco::Muon>> probes_;
      const edm::EDGetTokenT<std::vector<reco::Vertex>> vertices_;

      void writeValueMap(edm::Event &iEvent,
                const edm::Handle<edm::View<reco::Muon> > & handle,
                const std::vector<float> & values,
                const std::string    & label) const ;
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
MuonImpactParameter::MuonImpactParameter(const edm::ParameterSet& iConfig):
   probes_(consumes<edm::View<reco::Muon>>(iConfig.getParameter<edm::InputTag>("probes"))),
   vertices_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertices")))
{

  produces<edm::ValueMap<float> >("Dxy");
  produces<edm::ValueMap<float> >("DxyError");
  produces<edm::ValueMap<float> >("DxyOverDxyError");

  produces<edm::ValueMap<float> >("Dz");
  produces<edm::ValueMap<float> >("DzError");
  produces<edm::ValueMap<float> >("DzOverDzError");

}


MuonImpactParameter::~MuonImpactParameter()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
MuonImpactParameter::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;

  vector<float> muon_Dxy;
  vector<float> muon_DxyError;
  vector<float> muon_DxyOverDxyError;
  vector<float> muon_Dz;
  vector<float> muon_DzError;
  vector<float> muon_DzOverDzError;
  
  edm::Handle<vector<reco::Vertex> > primaryVerticesHandle;
  iEvent.getByToken(vertices_, primaryVerticesHandle);

  int bestVtx = -1;
  float maxSumPt = 0;

  //vertex loop
  for(unsigned int iVtx=0; iVtx<primaryVerticesHandle->size(); iVtx++){

    float vtxSumPt = 0;
    for (reco::Vertex::trackRef_iterator it = primaryVerticesHandle->at(iVtx).tracks_begin(); it != primaryVerticesHandle->at(iVtx).tracks_end(); it++) vtxSumPt += (**it).pt();

    if(vtxSumPt > maxSumPt) {
      maxSumPt = vtxSumPt;
      bestVtx = iVtx;
    }

  } //end of vertex loop

  Handle<View<reco::Muon> > probes;
  iEvent.getByToken(probes_,  probes);

  //muon loop
  for (View<reco::Muon>::const_iterator probe = probes->begin(); probe != probes->end(); ++probe) {
    
    float Dxy = 9999;
    float DxyError = 999;
    float Dz = 9999;
    float DzError = 999;

    if(probe->innerTrack().isNonnull() && bestVtx>=0) {
      
      edm::ESHandle<TransientTrackBuilder> ttBuilder;
      iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", ttBuilder);
      
      reco::TrackRef muonTrack = probe->innerTrack();

      Dxy = muonTrack->dxy(primaryVerticesHandle->at(bestVtx).position());
      DxyError = sqrt(muonTrack->dxyError()*muonTrack->dxyError() + primaryVerticesHandle->at(bestVtx).xError()*primaryVerticesHandle->at(bestVtx).yError());
      Dz = muonTrack->dz(primaryVerticesHandle->at(bestVtx).position());
      DzError = sqrt(muonTrack->dzError()*muonTrack->dzError() + primaryVerticesHandle->at(bestVtx).zError()*primaryVerticesHandle->at(bestVtx).zError() );

    }

    muon_Dxy.push_back(Dxy);
    muon_DxyError.push_back(DxyError);
    muon_DxyOverDxyError.push_back(fabs(Dxy/DxyError));
    muon_Dz.push_back(Dz);
    muon_DzError.push_back(DzError);
    muon_DzOverDzError.push_back(fabs(Dz/DzError));

  } //end of muon loop

  writeValueMap(iEvent,	probes,	muon_Dxy,		"Dxy");
  writeValueMap(iEvent,	probes,	muon_DxyError,		"DxyError");
  writeValueMap(iEvent,	probes,	muon_DxyOverDxyError,	"DxyOverDxyError");
  writeValueMap(iEvent,	probes,	muon_Dz,		"Dz");
  writeValueMap(iEvent,	probes,	muon_DzError,		"DzError");
  writeValueMap(iEvent,	probes,	muon_DzOverDzError,	"DzOverDzError");

}

void
MuonImpactParameter::writeValueMap(edm::Event &iEvent,
        const edm::Handle<edm::View<reco::Muon> > & handle,
        const std::vector<float> & values,
        const std::string    & label) const 
{
    using namespace edm; 
    using namespace std;
    auto_ptr<ValueMap<float> > valMap(new ValueMap<float>());
    edm::ValueMap<float>::Filler filler(*valMap);
    filler.insert(handle, values.begin(), values.end());
    filler.fill();
    iEvent.put(valMap, label);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonImpactParameter);
