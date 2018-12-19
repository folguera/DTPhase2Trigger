// -*- C++ -*-
//
// Package:    L1Trigger/DTPhase2Trigger
// Class:      DTPhase2Trigger
//
/**\class DTPhase2Trigger DTPhase2Trigger.cc L1Trigger/DTPhase2Trigger/plugins/DTPhase2Trigger.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Santiago Folgueras
//         Created:  Thu, 06 Dec 2018 09:44:33 GMT
//
//


// system include files
#include <memory>

// user include files (from CMSSW)
#include "L1Trigger/DTPhase2Trigger/plugins/DTPhase2TriggerTest.h"

// ROOT headers 
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"

//
// constants, enums and typedefs
//
using namespace std;
using namespace edm;
//
// static data member definitions
//

//
// constructors and destructor
//
DTPhase2Trigger::DTPhase2Trigger(const edm::ParameterSet& pset)
 :
  
  
{
   //now do what ever initialization is needed
  my_debug = pset.getUntrackedParameter<bool>("debug");
  my_DTTFnum = pset.getParameter<bool>("DTTFSectorNumbering");
  my_params = pset;
    
  digiLabel_ = pset.getParameter<edm::InputTag>("digiTag");
  dt4DSegments = consumes<DTRecSegment4DCollection>(pset.getParameter < edm::InputTag > ("dt4DSegments"));
  
  string outputfile = pset.getUntrackedParameter<string>("outputFileName");
  my_rootfile = new TFile(outputfile.c_str(),"RECREATE");
  my_tree = new TTree("tree","L1T",0);

  // get the tTrigDBInfo
  theSync = DTTTrigSyncFactory::get()->create(pset.getUntrackedParameter<std::string>("tTrigMode"),
					      pset.getUntrackedParameter<edm::ParameterSet>("tTrigModeConfig"));
  
}


DTPhase2Trigger::~DTPhase2Trigger()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  delete my_rootfile;
  if (my_debug) 
    cout << "[DTTrigTest] Destructor executed!!!" << endl;
  
}


//
// member functions
//

// ------------ method called once each job just before starting event loop  ------------
void
DTPhase2Trigger::beginJob()
{
  h_digiTDC       = new TH1F("h_digiTDC","h_digiTDC",1601,-0.5,1600.5);
  h_digiTDCPhase2 = new TH1F("h_digiTDCPhase2","h_digiTDCPhase2",3563*32,-0.5,3563*32+1);

  h_digiTime = new TH1F("h_digiTime","h_digiTime",1275,-0.5,1274.5);
  h_digiTimePhase2 = new TH1F("h_digiTimePhase2","h_digiTimePhase2",8907,-0.5,89075.5);    
     

}


// ------------ method called for each event  ------------
void
DTPhase2Trigger::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   Handle<DTDigiCollection> dtdigis;
   iEvent.getByLabel(digiLabel_, dtdigis);
   
   Handle<DTRecSegment4DCollection> all4DSegments;
   iEvent.getByToken(dt4DSegments, all4DSegments);


   // The code (for now), will focus in one chamber and
   // look at segments / digis there
   DTChamberId selected_chamber_ID(0,1,6);//ciemat_chamber

   // Select event if there is a 4D segment in desired chamber
   int DTSegmentCounterInChamber=0;
   vector<const DTRecSegment4D*> my_segments4D;
   for (DTRecSegment4DCollection::const_iterator segm = all4DSegments->begin(); segm!=all4DSegments->end(); ++segm){
     if (segm->chamberId()!=selected_chamber_ID) continue;
     if (!segm->hasPhi()) continue;
     if (!segm->hasZed()) continue;

     DTSegmentCounterInChamber++;
     my_segments4D.push_back(&(*segment1));
   }

   //   DTRecSegment4DCollection::const_iterator segment;
   //   for (segment = all4DSegments->begin();segment!=all4DSegments->end(); ++segment){}     
   
   /// FOCUS ON SELECTED CHAMBER:
   DTDigiCollection::DigiRangeIterator dtLayerId_It;
   for (dtLayerId_It=dtdigis->begin(); dtLayerId_It!=dtdigis->end(); ++dtLayerId_It){
     for (DTDigiCollection::const_iterator digiIt = ((*dtLayerId_It).second).first;digiIt!=((*dtLayerId_It).second).second; ++digiIt){
       const DTLayerId dtLId = (*dtLayerId_It).first;
      
       if (dtLId.wheel()!=selected_chamber_ID.wheel()) continue;
       if (dtLId.sector()!=selected_chamber_ID.sector()) continue;
       if (dtLId.station()!=selected_chamber_ID.station()) continue;

       int l  = dtLId.layer();
       int sl = dtLId.superlayer();
       
       int digiTDC = (*digiIt).countsTDC();
       int digiTDCPhase2 =  (*digiIt).countsTDC()+ iEvent.eventAuxiliary().bunchCrossing()*32;

       float ttrig = theSync->offset((*digiIt.wireId()));
       float digiTIME = (*digiIt).time();
       float digiTIMEPhase2 =  (*digiIt).time()+ iEvent.eventAuxiliary().bunchCrossing()*25-ttrig;//how to get the value of other station/chamber?
       
       
       h_digiTDC->Fill(digiTDC); 
       h_digiTDCPhase2->Fill(digiTDCPhase2);  
       
       h_digiTime->Fill(digiTIME); 
       h_digiTimePhase2->Fill(digiTIMEPhase2);
     }
   } // end DTDigiCollection
   
   
   /// NOW build the L1
}



// ------------ method called once each job just after ending the event loop  ------------
void
DTPhase2Trigger::endJob()
{
  
  my_rootfile->cd();

  h_digiTDC->Write();
  h_digiTDCPhase2->Write();
  
  h_digiTime->Write();
  h_digiTimePhase2->Write();  
    
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DTPhase2Trigger::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DTPhase2Trigger);
