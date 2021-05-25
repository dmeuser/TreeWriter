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
	    std::string bEffName_ee_ = params.getParameter<std::string>("bEffFullSimName_ee");
	    std::string cEffName_ee_ = params.getParameter<std::string>("cEffFullSimName_ee");
	    std::string lightEffName_ee_ = params.getParameter<std::string>("lightEffFullSimName_ee");
	    std::string bEffName_mumu_ = params.getParameter<std::string>("bEffFullSimName_mumu");
	    std::string cEffName_mumu_ = params.getParameter<std::string>("cEffFullSimName_mumu");
	    std::string lightEffName_mumu_ = params.getParameter<std::string>("lightEffFullSimName_mumu");
	    std::string bEffName_emu_ = params.getParameter<std::string>("bEffFullSimName_emu");
	    std::string cEffName_emu_ = params.getParameter<std::string>("cEffFullSimName_emu");
	    std::string lightEffName_emu_ = params.getParameter<std::string>("lightEffFullSimName_emu");
	    getBTagEffHistos(fileName_,bEffName_ee_,cEffName_ee_,lightEffName_ee_,bEffName_mumu_,cEffName_mumu_,lightEffName_mumu_,bEffName_emu_,cEffName_emu_,lightEffName_emu_);
    
  }
  
  
  
  const double operator()(int &flavor, double pt, double eta, int &channel){
	  
	  float result = 1.0;
	  
	  if (flavor == 5){
		  if(channel==1)result = bEfficiencies_ee_->GetBinContent(bEfficiencies_ee_->GetXaxis()->FindBin(pt),bEfficiencies_ee_->GetYaxis()->FindBin(eta));
		  else if(channel==2)result = bEfficiencies_mumu_->GetBinContent(bEfficiencies_mumu_->GetXaxis()->FindBin(pt),bEfficiencies_mumu_->GetYaxis()->FindBin(eta));
		  else result = bEfficiencies_emu_->GetBinContent(bEfficiencies_emu_->GetXaxis()->FindBin(pt),bEfficiencies_emu_->GetYaxis()->FindBin(eta));
	  }
	  else if (flavor == 4){
		  if(channel==1)result = cEfficiencies_ee_->GetBinContent(cEfficiencies_ee_->GetXaxis()->FindBin(pt),cEfficiencies_ee_->GetYaxis()->FindBin(eta));
		  else if(channel==2)result = cEfficiencies_mumu_->GetBinContent(cEfficiencies_mumu_->GetXaxis()->FindBin(pt),cEfficiencies_mumu_->GetYaxis()->FindBin(eta));
		  else result = cEfficiencies_emu_->GetBinContent(cEfficiencies_emu_->GetXaxis()->FindBin(pt),cEfficiencies_emu_->GetYaxis()->FindBin(eta));
	  }
	  else {
		  if(channel==1)result = lightEfficiencies_ee_->GetBinContent(lightEfficiencies_ee_->GetXaxis()->FindBin(pt),lightEfficiencies_ee_->GetYaxis()->FindBin(eta));
		  else if(channel==2)result = lightEfficiencies_mumu_->GetBinContent(lightEfficiencies_mumu_->GetXaxis()->FindBin(pt),lightEfficiencies_mumu_->GetYaxis()->FindBin(eta));
		  else result = lightEfficiencies_emu_->GetBinContent(lightEfficiencies_emu_->GetXaxis()->FindBin(pt),lightEfficiencies_emu_->GetYaxis()->FindBin(eta));
	  }
		
    return result;
  }
  
  
  
  

private:
  
  void getBTagEffHistos(std::string fileName, std::string bEffName_ee, std::string cEffName_ee, std::string lightEffName_ee, std::string bEffName_mumu, std::string cEffName_mumu, std::string lightEffName_mumu, std::string bEffName_emu, std::string cEffName_emu, std::string lightEffName_emu);

  //~ boost::shared_ptr<TFile> file_;
  TH2F * bEfficiencies_ee_;
  TH2F * cEfficiencies_ee_;
  TH2F * lightEfficiencies_ee_;
  
  TH2F * bEfficiencies_mumu_;
  TH2F * cEfficiencies_mumu_;
  TH2F * lightEfficiencies_mumu_;
  
  TH2F * bEfficiencies_emu_;
  TH2F * cEfficiencies_emu_;
  TH2F * lightEfficiencies_emu_;
};

#endif /* BTAGEFFMAPFUNCTOR_H_ */
