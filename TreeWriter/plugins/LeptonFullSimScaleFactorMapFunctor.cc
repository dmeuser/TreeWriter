#include "LeptonFullSimScaleFactorMapFunctor.h"
#include "TH2.h"
#include "TH2F.h"
#include "TH3.h"
#include "TH3F.h"
#include "TFile.h"
#include <string>

#include <DataFormats/PatCandidates/interface/Lepton.h>
#include <DataFormats/PatCandidates/interface/Electron.h>
#include <DataFormats/PatCandidates/interface/Muon.h>


void
LeptonFullSimScaleFactorMapFunctor::getFullSimScaleFactorHistos(std::string dataMCScaleFactorFile_mu_ID, 
                                                                                        std::string dataMCScaleFactorFile_mu_Iso, 
                                                                                        
                                                                                        std::string dataMCScaleFactorFile_ele,  
                                                                                        
                                                                                        std::string dataMCScaleFactorHisto_mu_ID, 
                                                                                        std::string dataMCScaleFactorHisto_mu_Iso, 
                                                                                        
                                                                                        std::string dataMCScaleFactorHisto_ele_ID 
                                                                                        )
{
        TFile dataMCScaleFactorFile_mu_ID_(dataMCScaleFactorFile_mu_ID.c_str());
        TFile dataMCScaleFactorFile_mu_Iso_(dataMCScaleFactorFile_mu_Iso.c_str());
        
        TFile dataMCScaleFactorFile_ele_(dataMCScaleFactorFile_ele.c_str());
        
        dataMCScaleFactorHisto_mu_ID_ = (TH2D*)dataMCScaleFactorFile_mu_ID_.Get( dataMCScaleFactorHisto_mu_ID.c_str() );
        dataMCScaleFactorHisto_mu_Iso_ = (TH2D*)dataMCScaleFactorFile_mu_Iso_.Get( dataMCScaleFactorHisto_mu_Iso.c_str() );
        
        dataMCScaleFactorHisto_ele_ID_ = (TH2D*)dataMCScaleFactorFile_ele_.Get( dataMCScaleFactorHisto_ele_ID.c_str() );
}

