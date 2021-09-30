// To make a new filter, copy this file under a new name in the "filter" directory.
// The "FilterCondition" function should return a boolean value.
// The boolean indicates whether the event should be kept or not.

#include "TRuntimeObjects.h"

#include <iostream>
#include <map>
#include <cstdio>

#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TRandom.h>
#include <TObject.h>
#include <TLine.h>

#include "TGretina.h"
#include "TS800.h"
#include "TBank29.h"
#include "TS800.h"
#include "TFastScint.h"
#include "GCutG.h"

#include "TChannel.h"
#include "GValue.h"

GCutG* incoming = 0;
GCutG* outgoing = 0;
GCutG* ic_corr = 0;
GCutG* time_energy = 0;
GCutG* xToF = 0;

int gates_loaded=0;

extern "C"
bool FilterCondition(TRuntimeObjects& obj) {
  
  TS800 *s800 = obj.GetDetector<TS800>();
  TBank29 *bank29 = obj.GetDetector<TBank29>(); 

  if(!s800 || !bank29)
    {return false;}  
   
  TList *gates = &(obj.GetGates());  
  
  if(gates_loaded!=gates->GetSize()) {
    TIter iter(gates);
    while(TObject *obj = iter.Next()) {
      GCutG *gate = (GCutG*)obj;
      std::string tag = gate->GetTag();
      if(!tag.compare("incoming")) {
        incoming = gate;
        std::cout << "incoming: << " << gate->GetName() << std::endl;
      } else if(!tag.compare("outgoing")) {
        outgoing = gate;
        std::cout << "outgoing: << " << gate->GetName() << std::endl;
      } else if(!tag.compare("time_energy")) {
        time_energy = gate;
        std::cout << "time energy: << " << gate->GetName() << std::endl;
      } else if(!tag.compare("x_tof")) {
        xToF = gate;
        std::cout << "xToF: << " << gate->GetName() << std::endl;
      } else if(!tag.compare("ic_corr")) {
        ic_corr = gate;
        std::cout << "IC_Corr: << " << gate->GetName() << std::endl;
      }
      gates_loaded++;
    }
  }
  
  //////////////////////////////////////////////////////////
  //Incoming Gates on Xfp - E1 vs. Object - E1 Uncorrected//
  //////////////////////////////////////////////////////////
/*
  if(incoming) {  
    if(!incoming->IsInside(s800->GetMTof().GetCorrelatedObjE1(),s800->GetMTof().GetCorrelatedXfpE1()))
      {return false;}
  }
*/
  ///////////////////////////////////////////////////////////////////
  //Outgoing Gates on Ion Chamber Charge vs Object - E1 (Corrected)//
  ///////////////////////////////////////////////////////////////////
  if(outgoing) {
    if(!outgoing->IsInside(s800->GetCorrTOF_OBJ(),s800->GetIonChamber().Charge()))
      {return false;}
  }
/*
  //////////////////////////////////////////////////////////
  //Correlation Gates on CRDC1X vs Object - E1 (Corrected)//
  //////////////////////////////////////////////////////////
  if(xToF) {
    if(!xToF->IsInside(s800->GetMTofObjE1(),s800->GetXFP()))
      {return false;}
  }

  //////////////////////////////////////////////////////////
  //PID Gates on Ion Chamber dE vs Object - E1 (Corrected)//
  //////////////////////////////////////////////////////////
  if(ic_corr) {
    if(!ic_corr->IsInside(s800->GetMTofObjE1(),s800->GetIonChamber().GetdE(s800->GetXFP(),s800->GetYFP())))
      {return false;}
  }
*/  
  return true;

}
