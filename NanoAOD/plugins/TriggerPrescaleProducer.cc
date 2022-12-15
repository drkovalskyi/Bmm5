#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"

#include <boost/regex.hpp>
// 
// TriggerPrescaleProducer that stores prescale values for each event
// for a set of triggers. If a trigger is not present prescale will be
// set to zero.
//

using namespace std;

///////////////////////////////////////////////////////////////////////////
///                             P L U G I N
///////////////////////////////////////////////////////////////////////////

class TriggerPrescaleProducer : public edm::stream::EDProducer<> {
    
public:
    
  explicit TriggerPrescaleProducer(const edm::ParameterSet &iConfig);
    
  ~TriggerPrescaleProducer() override {};
    
    
private:
    
  virtual void produce(edm::Event&, const edm::EventSetup&);
  void trigger_init(const edm::TriggerNames & triggerNames);
  
  // ----------member data ---------------------------
  edm::EDGetTokenT<edm::TriggerResults>         trigger_token;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> prescale_token;
    
  vector<string> triggers;
  edm::ParameterSetID trigger_names_id;
  std::map<string, unsigned int> trigger_map;

};

TriggerPrescaleProducer::TriggerPrescaleProducer(const edm::ParameterSet &iConfig):
  trigger_token(  consumes<edm::TriggerResults>(         iConfig.getParameter<edm::InputTag>("trigger"))),
  prescale_token( consumes<pat::PackedTriggerPrescales>( iConfig.getParameter<edm::InputTag>("prescales"))),
  triggers( iConfig.getParameter<std::vector<std::string> > ("triggerNames"))
{
  produces<nanoaod::FlatTable>();
} 

void TriggerPrescaleProducer::trigger_init(const edm::TriggerNames & triggerNames){
  trigger_map.clear();
  for (unsigned int trigger_id = 0; trigger_id < triggerNames.size(); ++trigger_id) {
    std::string triggerName = triggerNames.triggerName(trigger_id);
    for (auto const & required_trigger_name: triggers){
      if (boost::regex_match(triggerName, boost::regex("^" + required_trigger_name + "_v\\d+$"))){
	trigger_map[required_trigger_name] = trigger_id;
      }
    }
  }
}

void TriggerPrescaleProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  edm::Handle<edm::TriggerResults> trigger_handle;
  iEvent.getByToken(trigger_token, trigger_handle);

  edm::Handle<pat::PackedTriggerPrescales> prescale_handle;
  iEvent.getByToken(prescale_token,prescale_handle);

  
  if (not trigger_handle.isValid() or not prescale_handle.isValid())
    throw cms::Exception("Cannot get trigger information");
  
  edm::TriggerNames triggerNames = iEvent.triggerNames(*trigger_handle);
  if (trigger_names_id != triggerNames.parameterSetID()) {
    trigger_names_id = triggerNames.parameterSetID();
    trigger_init(triggerNames);
  }

  auto out  = std::make_unique<nanoaod::FlatTable>(1, "prescale", true);
  for (auto const & trigger: triggers){
    auto trigger_id = trigger_map.find(trigger);
    int prescale = 0;
    if (trigger_id != trigger_map.end()){
      prescale = prescale_handle->getPrescaleForIndex(trigger_id->second);
    }
    out->addColumnValue<int>(trigger, prescale, "HLT prescale");
  }
  iEvent.put(std::move(out));
}

DEFINE_FWK_MODULE(TriggerPrescaleProducer);
