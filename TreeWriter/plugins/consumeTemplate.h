/*
 *   consumeTemplate.h
 * 
 *  Created on: Nov 26, 2015
 *  Author: tarndt
 */


#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "version.h"

//This class allows compatibility between CMSSW76 and CMSSW5, by templating
//the consumes function. The functions are meant to do nothing for CMSSW5
//in the moment.
//These classes have to be inserted into the inheritance tree of the framework modules
//in order to work: The framework module has to inherit from the respective template class
//instead of the edm::ED* class.
//That allows to call the consumeTemplate function in the constructor.

#ifndef TEMPLATE_H_
#define TEMPLATE_H_

class ProducerTemplate : public edm::EDProducer{

public:
template<typename type> void consumeTemplate(edm::InputTag tag ){ 
#ifndef CMSSW_LEQ_5
consumes<type>(tag);
#endif

}
};

class AnalyzerTemplate : public edm::EDAnalyzer{

public:
template<typename type> void consumeTemplate(edm::InputTag tag){
#ifndef CMSSW_LEQ_5
consumes<type>(tag);
#endif

}
};

class FilterTemplate : public edm::EDFilter{

public:
template<typename type> void consumeTemplate(edm::InputTag tag){
#ifndef CMSSW_LEQ_5
consumes<type>(tag);
#endif

}
};


#endif


