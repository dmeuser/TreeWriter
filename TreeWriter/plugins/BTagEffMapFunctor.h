/*
 * BTagEffMapFunctor.h
 *
 *  Created on: 04.12.2015
 *      Author: cschomak
 */

#ifndef BTAGEFFMAPFUNCTOR_H_
#define BTAGEFFMAPFUNCTOR_H_

#include <FWCore/Framework/interface/Event.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>

#include <iostream>
#include <boost/shared_ptr.hpp>

#include "TFile.h"
#include "TH2.h"

class BTagEffMapFunctor
{
public:
  BTagEffMapFunctor( ){}
  BTagEffMapFunctor( edm::ParameterSet const & params ){
            
            std::string fileName_ = params.getParameter<std::string>("fullSimFile");
	    std::string bEffName_ = params.getParameter<std::string>("bEffFullSimName");
	    std::string cEffName_ = params.getParameter<std::string>("cEffFullSimName");
	    std::string lightEffName_ = params.getParameter<std::string>("lightEffFullSimName");
	    getBTagEffHistos(fileName_,bEffName_,cEffName_,lightEffName_);
    
  }
  
  
  
  const double operator()(int &flavor, double pt, double eta){
	  
	  float result = 1.0;
	  
	  if (flavor == 5){
		  result = bEfficiencies_->GetBinContent(bEfficiencies_->GetXaxis()->FindBin(pt),bEfficiencies_->GetYaxis()->FindBin(eta));
	  }
	  else if (flavor == 4){
		  result = cEfficiencies_->GetBinContent(cEfficiencies_->GetXaxis()->FindBin(pt),cEfficiencies_->GetYaxis()->FindBin(eta));
	  }
	  else {
		  result = lightEfficiencies_->GetBinContent(lightEfficiencies_->GetXaxis()->FindBin(pt),lightEfficiencies_->GetYaxis()->FindBin(eta));
	  }
		
    return result;
  }
  
  
  
  

private:
  
  void getBTagEffHistos(std::string fileName, std::string bEffName, std::string cEffName, std::string lightEffName);

  //~ boost::shared_ptr<TFile> file_;
  TH2F * bEfficiencies_;
  TH2F * cEfficiencies_;
  TH2F * lightEfficiencies_;
};

#endif /* BTAGEFFMAPFUNCTOR_H_ */
