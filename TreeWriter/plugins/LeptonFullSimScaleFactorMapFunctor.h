/*
 * LeptonFullSimScaleFactorMapFunctor.h
 *
 *  Created on: 04.12.2015
 *      Author: cschomak
 */

#ifndef LEPTONFULLSIMSCALEFACTORMAPFUNCTOR_H_
#define LEPTONFULLSIMSCALEFACTORMAPFUNCTOR_H_

#include <FWCore/Framework/interface/Event.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>

#include <DataFormats/PatCandidates/interface/Lepton.h>
#include <DataFormats/PatCandidates/interface/Electron.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/Tau.h>

#include <iostream>

#include "TFile.h"
#include "TH2.h"
#include "TH3.h"

#include "TreeParticles.hpp"

class LeptonFullSimScaleFactorMapFunctor
{
public:
  LeptonFullSimScaleFactorMapFunctor( ){}
  LeptonFullSimScaleFactorMapFunctor( edm::ParameterSet const & params ){
            
            std::string dataMCScaleFactorFile_mu_ID_ = params.getParameter<std::string>("dataMCScaleFactorFile_mu_ID");
            std::string dataMCScaleFactorFile_mu_Iso_ = params.getParameter<std::string>("dataMCScaleFactorFile_mu_Iso");
            
            std::string dataMCScaleFactorFile_ele_ = params.getParameter<std::string>("dataMCScaleFactorFile_ele");
            
            std::string dataMCScaleFactorHisto_mu_ID_ = params.getParameter<std::string>("dataMCScaleFactorHisto_mu_ID");
            std::string dataMCScaleFactorHisto_mu_Iso_ = params.getParameter<std::string>("dataMCScaleFactorHisto_mu_Iso");
            
            std::string dataMCScaleFactorHisto_ele_ID_ = params.getParameter<std::string>("dataMCScaleFactorHisto_ele_ID");            
            
      
      getFullSimScaleFactorHistos(dataMCScaleFactorFile_mu_ID_, 
                dataMCScaleFactorFile_mu_Iso_, 
                
                dataMCScaleFactorFile_ele_, 
                
                dataMCScaleFactorHisto_mu_ID_, 
                dataMCScaleFactorHisto_mu_Iso_, 
                
                dataMCScaleFactorHisto_ele_ID_ 
                      );
    
  }
  
  
  
  const double * operator()(const  tree::Electron &ele, double pt, double eta){

    static double result[2];
    result[0] = 1.;
    
    if(pt > 500.) pt= 499.;
    
    result[0] *= dataMCScaleFactorHisto_ele_ID_->GetBinContent(dataMCScaleFactorHisto_ele_ID_->GetXaxis()->FindBin(eta), dataMCScaleFactorHisto_ele_ID_->GetYaxis()->FindBin(pt));
    result[1] = dataMCScaleFactorHisto_ele_ID_->GetBinError(dataMCScaleFactorHisto_ele_ID_->GetXaxis()->FindBin(eta), dataMCScaleFactorHisto_ele_ID_->GetYaxis()->FindBin(pt));
    
    return result;
  }
  
  const double * operator()(const  tree::Muon &mu, double pt, double eta){

    static double result[2];
    result[0] = 1.;
    
    if(pt > 120.) pt= 119.;
    else if(pt < 20.) pt= 21.;
    
    double ID = dataMCScaleFactorHisto_mu_ID_->GetBinContent(dataMCScaleFactorHisto_mu_ID_->GetXaxis()->FindBin(fabs(eta)),dataMCScaleFactorHisto_mu_ID_->GetYaxis()->FindBin(fabs(pt)));
    double ISO = dataMCScaleFactorHisto_mu_Iso_->GetBinContent(dataMCScaleFactorHisto_mu_Iso_->GetXaxis()->FindBin(fabs(eta)),dataMCScaleFactorHisto_mu_Iso_->GetYaxis()->FindBin(fabs(pt)));
    result[0] = ID*ISO;
    
    double err_ID = dataMCScaleFactorHisto_mu_ID_->GetBinError(dataMCScaleFactorHisto_mu_ID_->GetXaxis()->FindBin(fabs(eta)),dataMCScaleFactorHisto_mu_ID_->GetYaxis()->FindBin(fabs(pt)));
    double err_ISO = dataMCScaleFactorHisto_mu_Iso_->GetBinError(dataMCScaleFactorHisto_mu_Iso_->GetXaxis()->FindBin(fabs(eta)),dataMCScaleFactorHisto_mu_Iso_->GetYaxis()->FindBin(fabs(pt)));
    result[1] = sqrt((ISO*err_ID)*(ISO*err_ID)+(ID*err_ISO)*(ID*err_ISO));
    
    return result;
  }
 
  

private:

  bool useFastSim_;
  
  void getFullSimScaleFactorHistos(std::string dataMCScaleFactorFile_mu_ID, 
              std::string dataMCScaleFactorFile_mu_Iso, 
              
              std::string dataMCScaleFactorFile_ele, 
              
              std::string dataMCScaleFactorHisto_mu_ID, 
              std::string dataMCScaleFactorHisto_mu_Iso, 
              
              std::string dataMCScaleFactorHisto_ele_ID
              );

  TH2D * dataMCScaleFactorHisto_mu_ID_;
  TH2D * dataMCScaleFactorHisto_mu_Iso_;
  
  TH2D * dataMCScaleFactorHisto_ele_ID_;
};

#endif /* LEPTONFULLSIMSCALEFACTORMAPFUNCTOR_H_ */
