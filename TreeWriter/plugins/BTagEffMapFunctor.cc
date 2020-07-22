#include "BTagEffMapFunctor.h"
#include "TH2.h"
#include "TH2F.h"
#include "TFile.h"
#include <boost/shared_ptr.hpp>
#include <string>



void
BTagEffMapFunctor::getBTagEffHistos(std::string fileName, std::string bEffName, std::string cEffName, std::string lightEffName)
{
	TFile file_(fileName.c_str());
	
	bEfficiencies_ = (TH2F*)file_.Get( bEffName.c_str() );
	cEfficiencies_ = (TH2F*)file_.Get( cEffName.c_str() );
	lightEfficiencies_ = (TH2F*)file_.Get( lightEffName.c_str() );
}



