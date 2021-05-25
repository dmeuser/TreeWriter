#include "BTagEffMapFunctor.h"
#include "TH2.h"
#include "TH2F.h"
#include "TFile.h"
#include <boost/shared_ptr.hpp>
#include <string>



void
BTagEffMapFunctor::getBTagEffHistos(std::string fileName, std::string bEffName_ee, std::string cEffName_ee, std::string lightEffName_ee, std::string bEffName_mumu, std::string cEffName_mumu, std::string lightEffName_mumu, std::string bEffName_emu, std::string cEffName_emu, std::string lightEffName_emu)
{
	TFile file_(fileName.c_str());
	
	bEfficiencies_ee_ = (TH2F*)file_.Get( bEffName_ee.c_str() );
	cEfficiencies_ee_ = (TH2F*)file_.Get( cEffName_ee.c_str() );
	lightEfficiencies_ee_ = (TH2F*)file_.Get( lightEffName_ee.c_str() );
	
	bEfficiencies_mumu_ = (TH2F*)file_.Get( bEffName_mumu.c_str() );
	cEfficiencies_mumu_ = (TH2F*)file_.Get( cEffName_mumu.c_str() );
	lightEfficiencies_mumu_ = (TH2F*)file_.Get( lightEffName_mumu.c_str() );
	
	bEfficiencies_emu_ = (TH2F*)file_.Get( bEffName_emu.c_str() );
	cEfficiencies_emu_ = (TH2F*)file_.Get( cEffName_emu.c_str() );
	lightEfficiencies_emu_ = (TH2F*)file_.Get( lightEffName_emu.c_str() );
}



