#ifndef DTPHASE2TRIGGERTEST_H
#define DTPHASE2TRIGGERTEST_H
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

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


class DTPhase2Trigger : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit DTPhase2Trigger(const edm::ParameterSet&);
      ~DTPhase2Trigger();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      
  
};

#endif
