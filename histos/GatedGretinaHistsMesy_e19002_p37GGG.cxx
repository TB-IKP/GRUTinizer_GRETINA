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

class BinCenters {
    public:
        BinCenters(int bins, double low, double high)
            : bins(bins), low(low), high(high) { }

        class iterator {
            public:
                iterator(BinCenters& axis, int binnum)
                    : axis(axis), binnum(binnum) { }

                double operator*() const;

                bool operator==(const iterator& other) const {
                    return
                        (&axis == &other.axis) &&
                        (binnum == other.binnum);
                }

                bool operator!=(const iterator& other) const {
                    return !(*this == other);
                }

                iterator& operator++() {
                    binnum++;
                    return *this;
                }

            private:
                BinCenters& axis;
                int binnum;
        };
        friend class BinCenters::iterator;

        iterator begin() {
            return iterator(*this, 0);
        }

        iterator end() {
            return iterator(*this, bins);
        }

    private:
        int bins;
        double low;
        double high;
};

double BinCenters::iterator::operator*() const {
    return axis.low + (binnum+0.5) * (axis.high-axis.low)/axis.bins;
}

std::map<int,int> detMap = {{26,0}, {30,1}, {34,2}, {38,3}, {25,4}, {29,5}, {33,6}, {37,7}, {27,8}, {31,9}, {35,10}, {39,11},
			    {24,12},{28,13},{32,14},{36,15},{47,16},{51,17},{59,18},{63,19},{67,20},{71,21}, {79,22}, {44,23},
			    {50,24}, {58,25}, {60,26}, {66,27}, {68,28}, {76,29}, {46,30}, {48,31}, {56,32}, {62,33},
			    {64,34}, {70,35}, {78,36}, {45,37}, {49,38}, {57,39}, {61,40}, {65,41}, {69,42}, {77,43}};
/*
std::map<int,int> detMap = {{26,1}, {30,2}, {34,3}, {38,4}, {25,5}, {29,6}, {33,7}, {37,8}, {27,9}, {31,10}, {35,11}, {39,12},
			    {24,13},{28,14},{32,15},{36,16},{63,17},{71,18},{79,19},{59,20},{67,21},{58,22}, {66,23}, {60,24},
			    {68,25}, {76,26}, {62,27}, {70,28}, {78,29}, {56,30}, {64,31}, {57,32}, {65,33}, {61,34},
			    {69,35}, {77,36}};
*/

std::map<int,int> quadMap = {{5,1}, {6,2}, {7,3}, {8,4}, {14,5}, {16,6}, {18,7}, {13,8}, {15,9}};

std::map<int,int> crysThetaMap = {{26,1}, {30,1}, {34,1}, {38,1}, {25,2}, {29,2}, {33,2}, {37,2}, {27,3}, {31,3}, {35,3}, {39,3},
				  {24,4}, {28,4}, {32,4}, {36,4}, {47,5}, {51,6}, {59,6}, {63,5}, {67,6}, {71,5}, {79,5}, {44,8},
				  {50,7}, {58,7}, {60,8}, {66,7}, {68,8}, {76,8}, {46,9}, {48,10}, {56,10}, {62,9}, {64,10},
				  {70,9}, {78,9}, {45,12}, {49,11}, {57,11}, {61,12}, {65,11}, {69,12}, {77,11}};

/*
bool HandleTiming_Gated(TRuntimeObjects &obj, TCutG *incoming, TCutG* outgoing, TCutG* xToF) {
  TS800 *s800  = obj.GetDetector<TS800>();
  
  if(!s800)
    {return false;}

  std::string dirname = "Timing";

  //////////////////////////////////////////////////////////
  //incoming gates on Xfp - E1 vs. Object - E1 Uncorrected//
  //////////////////////////////////////////////////////////
  if(incoming) {  
    if(!incoming->IsInside(s800->MCorrelatedOBJ_E1(),s800->MCorrelatedXFP_E1()))
      {return false;}
    dirname += Form("_%s",incoming->GetName());
  }

  std::string suf = "";

  ///////////////////////////////////////////////////////////////////
  //Outgoing Gates on Ion Chamber Charge vs Object - E1 (Corrected)//
  ///////////////////////////////////////////////////////////////////
  if(outgoing) {
    if(!outgoing->IsInside(s800->GetMTOF_ObjE1(),s800->GetIonChamber().Charge()))
      {return false;}
    suf += Form("_%s",outgoing->GetName());
  }

  //////////////////////////////////////////////////////////
  //Correlation Gates on CRDC1X vs Object - E1 (Corrected)//
  //////////////////////////////////////////////////////////
  if(xToF) {
    if(!xToF->IsInside(s800->GetMTOF_ObjE1(),s800->GetXFP()))
      {return false;}
    suf += Form("_%s",xToF->GetName());
  }

  //int E1UpSize = s800->GetTof().E1UpSize();
  //int XfpSize = s800->GetTof().XfpSize();
  //int ObjSize = s800->GetTof().ObjSize();

  //obj.FillHistogram(dirname,Form("E1UpSize%s",suf.c_str()),10,0,10,E1UpSize);

  //obj.FillHistogram(dirname,Form("XfpSize%s",suf.c_str()),10,0,10,XfpSize);

  //obj.FillHistogram(dirname,Form("ObjSize%s",suf.c_str()),10,0,10,ObjSize);

  obj.FillHistogram(dirname,Form("E1_raw%s",suf.c_str()),4000,0,64000,s800->GetTof().GetRF());

  obj.FillHistogram(dirname,Form("Xfp_raw%s",suf.c_str()),4000,0,64000,s800->GetTof().GetXFP());

  obj.FillHistogram(dirname,Form("Obj_raw%s",suf.c_str()),4000,0,64000,s800->GetTof().GetOBJ());

  obj.FillHistogram(dirname,Form("Xfp-E1%s",suf.c_str()),
		    4000,-10000,10000,s800->MCorrelatedXFP_E1());

  obj.FillHistogram(dirname,Form("Obj-E1%s",suf.c_str()),
		    4000,-10000,10000,s800->MCorrelatedOBJ_E1());

  obj.FillHistogram(dirname,Form("Obj-Xfp%s",suf.c_str()),
		    800,-5200,-4400,s800->MCorrelatedOBJ_E1() -
		    s800->MCorrelatedXFP_E1());

  obj.FillHistogram(dirname,Form("Xfp-E1_raw%s",suf.c_str()),
                    4000,-10000,10000,s800->GetTof().GetXFP() - s800->GetTof().GetRF());
  
  obj.FillHistogram(dirname,Form("Obj-E1_raw%s",suf.c_str()),
                    4000,-10000,10000,s800->GetTof().GetOBJ() - s800->GetTof().GetRF());

  return true;

}
*/
bool HandleS800_Gated(TRuntimeObjects &obj,TCutG *incoming, TCutG* outgoing, TCutG* xToF) {
  TS800 *s800  = obj.GetDetector<TS800>();
  TGretina *gretina  = obj.GetDetector<TGretina>();
  
  if(!s800 || !incoming || !outgoing)
    {return false;}

  if(!gretina)
    {return false;}

  //if(!incoming || !outgoing || !xToF)
  //{return false;}

  std::string dsuf = "";

  ///////////////////////////////////////////////////////////
  //Incoming  Gates on Xfp - E1 vs. Object - E1 Uncorrected//
  ///////////////////////////////////////////////////////////
  if(incoming) {  
    if(!incoming->IsInside(s800->GetMTof().GetCorrelatedObjE1(),s800->GetMTof().GetCorrelatedXfpE1()))
      {return false;}
    dsuf += Form("_%s",incoming->GetName());
  }

  std::string hsuf = "";
  ///////////////////////////////////////////////////////////////////
  //Outgoing Gates on Ion Chamber Charge vs Object - E1 (Corrected)//
  ///////////////////////////////////////////////////////////////////
  if(outgoing) {
    if(!outgoing->IsInside(s800->GetMTofObjE1(),s800->GetIonChamber().Charge()))
    //if(!outgoing->IsInside(s800->GetOBJRaw_TAC(),s800->GetIonChamber().Charge()))
      {return false;}
    hsuf += Form("_%s",outgoing->GetName());
  }

  //////////////////////////////////////////////////////////
  //Correlation Gates on CRDC1X vs Object - E1 (Corrected)//
  //////////////////////////////////////////////////////////
  //if(xToF) {
  //  if(!xToF->IsInside(s800->GetMTofXfpE1(),s800->GetXFP()))
  //    {return false;}
  //  hsuf += Form("_%s",xToF->GetName());
  //}

  std::string dirname = "";
  
  //Incoming Beam and Outgoing PIDs
  dirname = "S800" + dsuf;
  
//Histogram ToF vs XTof
  obj.FillHistogram(dirname,Form("Xfp-E1_v_Obj-E1_Uncorrected%s",hsuf.c_str()),
		    1000,-2500,2500,s800->GetMTof().GetCorrelatedObjE1(),
		    1000,-2500,5000,s800->GetMTof().GetCorrelatedXfpE1());

  obj.FillHistogram(dirname,Form("Xfp-E1_v_Obj-E1_Corrected%s",hsuf.c_str()),
                    1000,-2500,2500,s800->GetMTofObjE1(),
                    1000,-2500,5000,s800->GetMTofXfpE1());

//obj.FillHistogram(dirname,Form("Take2_Xfp-E1_v_Obj-E1_Uncorrected%s",hsuf.c_str()),
//		    2000,-2000,2000,s800->GetTof().GetOBJ() - s800->GetTof().GetRF(),
//		    2000,-4000,4000,s800->GetTof().GetXFP() - s800->GetTof().GetRF());

//Histogram Charge vs Tof
  obj.FillHistogram(dirname,Form("Charge_v_XFP-E1Raw%s",hsuf.c_str()),
                  1000,-2500,5000,s800->GetMTof().GetCorrelatedXfpE1(),
                  1000,0,30000,s800->GetIonChamber().Charge());

  obj.FillHistogram(dirname,Form("Charge_v_Obj-E1Raw%s",hsuf.c_str()),
                  1000,-2500,2500,s800->GetMTof().GetCorrelatedObjE1(),
                  1000,0,30000,s800->GetIonChamber().Charge());

  obj.FillHistogram(dirname,Form("Charge_v_Obj-E1%s",hsuf.c_str()),
                  1000,-2500,2500,s800->GetMTofObjE1(),
                  1000,0,30000,s800->GetIonChamber().Charge());

  obj.FillHistogram(dirname,Form("Charge_v_XToF-E1%s",hsuf.c_str()),
                  1000,-2500,5000,s800->GetMTofXfpE1(),
                  1000,0,30000,s800->GetIonChamber().Charge());
/*
  obj.FillHistogram(dirname,Form("Charge_v_Obj-E1_AoI%s",hsuf.c_str()),
                  200,-550,-350,s800->GetCorrTOF_OBJ(),
                  700,12000,19000,s800->GetIonChamber().Charge());

  obj.FillHistogram(dirname,Form("Charge_v_XToF-E1_AoI%s",hsuf.c_str()),
                  200,-1200,-1000,s800->GetCorrTOF_XFP(),
                  700,12000,19000,s800->GetIonChamber().Charge());
*/
//obj.FillHistogram(dirname,Form("Take2_Charge_v_Obj-E1%s",hsuf.c_str()),
//		    2000,-2000,2000,s800->GetTof().GetOBJ() - s800->GetTof().GetRF(),
//                    3000,0,30000,s800->GetIonChamber().Charge());

  
  //obj.FillHistogram(dirname,Form("dE_v_Obj-E1%s",hsuf.c_str()),
  //		    1400,-1500,-800,s800->GetMTofObjE1(),
  //                550,1200,2300,s800->GetIonChamber().GetdE(s800->GetXFP(),s800->GetYFP()));
//s800->GetOBJ_E1Raw()
//Histogram Charge vs XTof
//obj.FillHistogram(dirname,Form("Charge_v_CRDC1X%s",hsuf.c_str()),
//    1000,-2000,2000,s800->GetXFP(),
//                  1000,0,30000,s800->GetIonChamber().Charge());
//
//obj.FillHistogram(dirname,Form("CRDC1X_v_Obj-E1%s",hsuf.c_str()),
//                  1000,-2000,2000,s800->GetOBJ_E1Raw(),
//                  1250,-500,2000,s800->GetXFP());

  //obj.FillHistogram(dirname,Form("dE_v_CRDC1X%s",hsuf.c_str()),
  //		    300,-300,300,s800->GetXFP(),
  //                550,1200,2300,s800->GetIonChamber().GetdE(s800->GetXFP(),s800->GetYFP()));
  //Register and Trig Bit
  dirname = "S800" + dsuf;
  
  obj.FillHistogram(dirname,Form("Register%s",hsuf.c_str()),10,0,10,s800->GetReg());
  
  //unsigned int reg = s800->GetReg();
  unsigned int Sing = (s800->GetReg())&1;
  unsigned int Coin = (s800->GetReg())&2;

  if(Sing!=0) 
    {obj.FillHistogram(dirname,Form("TrigBit%s",hsuf.c_str()),10,0,10,Sing);}

  if(Coin!=0) 
    {obj.FillHistogram(dirname,Form("TrigBit%s",hsuf.c_str()),10,0,10,Coin);} 

  //CRDC Positions
  dirname = "CRDCPos" + dsuf;
  
  obj.FillHistogram(dirname,Form("CRDC1X%s",hsuf.c_str()),300,-1000,1000,s800->GetXFP());

  obj.FillHistogram(dirname,Form("CRDC1Y%s",hsuf.c_str()),300,-1000,5000,s800->GetYFP());

  obj.FillHistogram(dirname,Form("CRDC2X%s",hsuf.c_str()),300,-1000,1000,s800->GetXFP(1));

  obj.FillHistogram(dirname,Form("CRDC2Y%s",hsuf.c_str()),300,-1000,1000,s800->GetYFP(1));
  obj.FillHistogram(dirname,Form("CRDC1Y_v_CRDC1X%s",hsuf.c_str()),
  		    300,0,2000,s800->GetYFP(),
  		    600,-300,300,s800->GetXFP());

  obj.FillHistogram(dirname,Form("CRDC2Y_v_CRDC2X%s",hsuf.c_str()),
  		    300,-300,300,s800->GetYFP(1),
  		    600,-300,300,s800->GetXFP(1));
  obj.FillHistogram(dirname,Form("MaxPad%s",hsuf.c_str()),
                    100,0,500,s800->GetCrdc().GetMaxPad());
  obj.FillHistogram(dirname,Form("MaxPadSum%s",hsuf.c_str()),
                    100,0,500,s800->GetCrdc().GetMaxPadSum());
  //Focal Plane and Target
/*
  dirname = "FPandTA" + dsuf;
  
  obj.FillHistogram(dirname,Form("AFP%s",hsuf.c_str()),1000,-0.1,0.1,s800->GetAFP());

  obj.FillHistogram(dirname,Form("BFP%s",hsuf.c_str()),1000,-0.5,0.5,s800->GetBFP());

  obj.FillHistogram(dirname,Form("AFP_v_BFP%s",hsuf.c_str()),
  		    1000,-0.5,0.5,s800->GetBFP(),
  		    1000,-0.1,0.1,s800->GetAFP());

  obj.FillHistogram(dirname,Form("AFP_v_XFP%s",hsuf.c_str()),
  		    600,-300,300,s800->GetXFP(),
  		    1000,-0.1,0.1,s800->GetAFP());
  
  obj.FillHistogram(dirname,Form("ATA%s",hsuf.c_str()),1000,-0.1,0.1,s800->GetAta());

  obj.FillHistogram(dirname,Form("BTA%s",hsuf.c_str()),1000,-0.1,0.1,s800->GetBta());
  
  obj.FillHistogram(dirname,Form("ATA_v_BTA%s",hsuf.c_str()),
                    1000,-0.1,0.1,s800->GetBta(),
  		    1000,-0.1,0.1,s800->GetAta());

  obj.FillHistogram(dirname,Form("DTA%s",hsuf.c_str()),1000,-0.5,0.5,s800->GetDta());
  

  obj.FillHistogram(dirname,Form("YTA%s",hsuf.c_str()),1000,-50,50,s800->GetYta());
  
  obj.FillHistogram(dirname,Form("YFP_v_BFP%s",hsuf.c_str()),
  		    1000,-0.1,0.1,s800->GetBFP(),
  		    600,-300,300,s800->GetYFP());
  
  obj.FillHistogram(dirname,Form("Azita%s",hsuf.c_str()),400,-6.3,6.3,s800->Azita());
  //Focal Plane vs Time of Flight
  dirname = "FPvsToF" + dsuf;
  obj.FillHistogram(dirname,Form("XFP_v_Xfp-E1%s",hsuf.c_str()),
		    600,-2500,5000,s800->GetMTofXfpE1(),
                    600,-300,300,s800->GetXFP());

  obj.FillHistogram(dirname,Form("XFP_v_Xfp-E1_Uncorrected%s",hsuf.c_str()),
		    600,-2500,5000,s800->GetMTof().GetCorrelatedXfpE1(),
                    600,-300,300,s800->GetXFP());

  obj.FillHistogram(dirname,Form("AFP_v_Xfp-E1%s",hsuf.c_str()),
		    400,-2500,5000,s800->GetMTofXfpE1(),
		    1000,-0.05,0.05,s800->GetAFP());

  obj.FillHistogram(dirname,Form("AFP_v_Xfp-E1_Uncorrected%s",hsuf.c_str()),
		    400,-2500,5000,s800->GetMTof().GetCorrelatedXfpE1(),
                    1000,-0.05,0.05,s800->GetAFP());

  obj.FillHistogram(dirname,Form("XFP_v_Obj-E1%s",hsuf.c_str()),
		    600,-2500,2500,s800->GetMTofObjE1(),
                    600,-300,300,s800->GetXFP());

  obj.FillHistogram(dirname,Form("XFP_v_Obj-E1_Uncorrected%s",hsuf.c_str()),
		    600,-2500,2500,s800->GetMTof().GetCorrelatedObjE1(),
                    600,-300,300,s800->GetXFP());

  obj.FillHistogram(dirname,Form("AFP_v_Obj-E1%s",hsuf.c_str()),
		    400,-2500,2500,s800->GetMTofObjE1(),
		    1000,-0.05,0.05,s800->GetAFP());

  obj.FillHistogram(dirname,Form("AFP_v_Obj-E1_Uncorrected%s",hsuf.c_str()),
		    400,-2500,2500,s800->GetMTof().GetCorrelatedObjE1(),
                    1000,-0.05,0.05,s800->GetAFP());
*/
  return true;
}
/*
bool Gated_CS_Spectra(TRuntimeObjects &obj,TCutG *incoming, TCutG* outgoing, TCutG* xToF) {
  TS800 *s800  = obj.GetDetector<TS800>();

  if(!s800)
    {return false;}

  std::string dirname = "CrossSection";

  ///////////////////////////////////////////////////////////
  //Incoming  Gates on Xfp - E1 vs. Object - E1 Uncorrected//
  ///////////////////////////////////////////////////////////
  if(incoming) {  
    if(!incoming->IsInside(s800->MCorrelatedOBJ_E1(),s800->MCorrelatedXFP_E1()))
      {return false;}
    dirname += Form("_%s",incoming->GetName());
  }

  std::string suf = "";

  ///////////////////////////////////////////////////////////////////
  //Outgoing Gates on Ion Chamber Charge vs Object - E1 (Corrected)//
  ///////////////////////////////////////////////////////////////////
  if(outgoing) {
    if(!outgoing->IsInside(s800->GetMTOF_ObjE1(),s800->GetIonChamber().Charge()))
      {return false;}
    suf += Form("_%s",outgoing->GetName());
  }

  //////////////////////////////////////////////////////////
  //Correlation Gates on CRDC1X vs Object - E1 (Corrected)//
  //////////////////////////////////////////////////////////
  if(xToF) {
    if(!xToF->IsInside(s800->GetMTOF_ObjE1(),s800->GetXFP()))
      {return false;}
    suf += Form("_%s",xToF->GetName());
  }

  obj.FillHistogram(dirname,Form("Register%s",suf.c_str()),10,0,10,s800->GetReg());

  unsigned int reg = s800->GetReg();
  unsigned int Sing = reg&1;
  unsigned int Coin = reg&2;

  if(Sing!=0) 
    {obj.FillHistogram(dirname,Form("TrigBit%s",suf.c_str()),10,0,10,Sing);}

  if(Coin!=0) 
    {obj.FillHistogram(dirname,Form("TrigBit%s",suf.c_str()),10,0,10,Coin);} 

  obj.FillHistogram(dirname,Form("IC_Charge%s",suf.c_str()),
		    1600,20000,35000,s800->GetIonChamber().Charge());

  obj.FillHistogram(dirname,Form("AFP_v_Obj-Xfp_Uncorrected%s",suf.c_str()),
		    600,-5100,-4500,s800->MCorrelatedOBJ_E1()-s800->MCorrelatedXFP_E1(),
		    500,-0.05,0.05,s800->GetAFP());

  return true;

}
*/

bool HandleGretina(TRuntimeObjects &obj) {
  TGretina *gretina  = obj.GetDetector<TGretina>(); 
  //TBank29 *bank29 = obj.GetDetector<TBank29>();
   
  if(!gretina)
    {return false;}

  std::string dirname = "Gretina_Ungated";

  obj.FillHistogram(dirname,"Multiplicity",25,0,25,gretina->Size());

  for(unsigned int i=0;i<gretina->Size();i++) {
    TGretinaHit hit = gretina->GetGretinaHit(i);

    obj.FillHistogram(dirname,"GretinaCoreEnergy_NoThresh",2000,0,4000,hit.GetCoreEnergy());
    //if(bank29) {
    //  obj.FillHistogram(dirname,"CoreEn_v_TSDiff",125,-500,1500,bank29->Timestamp() - hit.GetTimestamp(),250,0,4000,hit.GetCoreEnergy());
    //}
    if(hit.GetCoreEnergy() > 100) {

      obj.FillHistogram(dirname,"GretinaCoreEnergy",2000,0,4000,hit.GetCoreEnergy());
      obj.FillHistogram(dirname,"GretinaSummarySpectrum",200,0,200,hit.GetCrystalId()
			                                ,2000,0,4000,hit.GetCoreEnergy());
                                                 
      obj.FillHistogram(dirname,"GretinaPositionSpectrum",226,0,3.2,hit.GetTheta()
                                                         ,452,0,6.3,hit.GetPhi());

      obj.FillHistogram(dirname,"HitTheta_v_DetMap",100,0,100,detMap[hit.GetCrystalId()]
			                              ,226,0,3.2,hit.GetTheta());

      obj.FillHistogram(dirname,"HitTheta_v_CrystalId",100,0,100,hit.GetCrystalId()
			                              ,226,0,3.2,hit.GetTheta());

      obj.FillHistogram(dirname,"HitPhi_v_DetMap",100,0,100,detMap[hit.GetCrystalId()]
			                         ,452,0,6.3,hit.GetPhi());

      obj.FillHistogram(dirname,"HitTheta_v_CrysThetaMap",100,0,100,crysThetaMap[hit.GetCrystalId()]
			                                 ,226,0,3.2,hit.GetTheta());

      obj.FillHistogram(dirname,"HitPhi_v_CrysThetaMap",100,0,100,crysThetaMap[hit.GetCrystalId()]
			                               ,452,0,6.3,hit.GetPhi());

      obj.FillHistogram(dirname,"HitTheta_v_QuadMap",100,0,100,quadMap[hit.GetHoleNumber()]
			                               ,226,0,3.2,hit.GetTheta());

      obj.FillHistogram(dirname,"HitPhi_v_QuadMap",100,0,100,quadMap[hit.GetHoleNumber()]
			                          ,452,0,6.3,hit.GetPhi());

    } //end if(hit.GetCoreEnergy() > 100)
  } //end gretina hit loop

  return true;
}

bool Gated_Gretina_Spectra(TRuntimeObjects &obj, TGretinaHit hit, TS800* s800, TBank29* bank29, int mult, std::string dsuf,
			   std::string hsuf, bool ab, double beta) {
    
  std::string pref = "";
  if(ab)
    {pref = "AB_";}
  
  std::string dirname = pref + "Gretina" + dsuf;


//Track is a complicated reconstruction of particle flight paths to get better accuracy on Doppler corrections. Right now I do not know how to set it up so it is commented out.
  TVector3 track = s800->Track();
  double gretina_ata_offset = GValue::Value("GRETINA_ATA_OFFSET");
  double gretina_bta_offset = GValue::Value("GRETINA_BTA_OFFSET");
  TVector3 adjusted_track = s800->Track(gretina_ata_offset,gretina_bta_offset);
  double yta = s800->GetYta();
/*
  for(double b : BinCenters(100,0.3,0.4)) {
   obj.FillHistogram(dirname,Form("Gamma(Beta)_BetaScan%s",hsuf.c_str()),
                     100,0.3,0.4,b,
                     500,0,4000,hit.GetDoppler(b));
   obj.FillHistogram(dirname,Form("Gamma(Beta&Corrections)_BetaScan%s",hsuf.c_str()),
                     100,0.3,0.4,b,
                     500,0,4000,hit.GetDopplerYta(b,yta,&track));
  }
  for(double b : BinCenters(100,0.3,0.4)) {
   obj.FillHistogram(dirname,Form("Gamma(Beta)_BetaScan%s",hsuf.c_str()),
                     100,0.3,0.4,b,
                     1000,0,4000,hit.GetDoppler(b,&track));
   if(quadMap[hit.GetHoleNumber()] <= 4) {
   obj.FillHistogram(dirname,Form("Forward_Gamma(Beta)_BetaScan%s",hsuf.c_str()),
                     100,0.3,0.4,b,
                     1000,0,4000,hit.GetDoppler(b,&track));
   }
   else {
   obj.FillHistogram(dirname,Form("90_Gamma(Beta)_BetaScan%s",hsuf.c_str()),
                     100,0.3,0.4,b,
                     1000,0,4000,hit.GetDoppler(b,&track));
   }
  }
*/
/*
  for(double b : BinCenters(150,0.25,0.4)) {
   obj.FillHistogram(dirname,Form("Gamma(Beta&Track&YTA)_BetaScan%s",hsuf.c_str()),
                     150,0.25,0.4,b,
                     1000,0,4000,hit.GetDopplerYta(b,yta,&track));
   if(quadMap[hit.GetHoleNumber()] <= 4) {
   obj.FillHistogram(dirname,Form("Forward_Gamma(Beta&Track&YTA)_BetaScan%s",hsuf.c_str()),
                     150,0.25,0.4,b,
                     1000,0,4000,hit.GetDopplerYta(b,yta,&track));
   }
   else {
   obj.FillHistogram(dirname,Form("90_Gamma(Beta&Track&YTA)_BetaScan%s",hsuf.c_str()),
                     150,0.25,0.4,b,
                     1000,0,4000,hit.GetDopplerYta(b,yta,&track));
   }
  }
*/
//Histogram the Doppler-corrected detection energies with 4 and 8 keV resolution.
//  obj.FillHistogram(dirname,Form("Gamma_Raw%s",hsuf.c_str()),500,0,4000,hit.GetCoreEnergy());
  obj.FillHistogram(dirname,Form("Gamma(Beta)%s",hsuf.c_str()),2000,0,4000,hit.GetDoppler(beta));
//obj.FillHistogram(dirname,Form("Double_Resolution_Gamma(Beta)%s",hsuf.c_str()),1000,0,4000,hit.GetDoppler(beta));
//obj.FillHistogram(dirname,Form("Half_Resolution_Gamma(Beta)%s",hsuf.c_str()),250,0,4000,hit.GetDoppler(beta));
  obj.FillHistogram(dirname,Form("Gamma(Beta&Track)%s",hsuf.c_str()),2000,0,4000,hit.GetDoppler(beta,&track));
//obj.FillHistogram(dirname,Form("Double_Resolution_Gamma(Beta&Track)%s",hsuf.c_str()),1000,0,4000,hit.GetDoppler(beta,&track));
//obj.FillHistogram(dirname,Form("Half_Resolution_Gamma(Beta&Track)%s",hsuf.c_str()),250,0,4000,hit.GetDoppler(beta,&track));
  obj.FillHistogram(dirname,Form("Gamma(Beta&Track&YTA)%s",hsuf.c_str()),2000,0,4000,hit.GetDopplerYta(beta,yta,&track));
  obj.FillHistogram(dirname,Form("Gamma(Beta&Track&YTA&DTA)%s",hsuf.c_str()),2000,0,4000,hit.GetDopplerYta(s800->AdjustedBeta(beta),yta,&track));
  obj.FillHistogram(dirname,Form("Gamma(Beta&AdjustedTrack&YTA&DTA)%s",hsuf.c_str()),2000,0,4000,hit.GetDopplerYta(s800->AdjustedBeta(beta),yta,&adjusted_track));
/*
  obj.FillHistogram(dirname,Form("Energy_vs_YTA%s",hsuf.c_str()),
                    20,-10,10,yta,
                    1000,0,4000,hit.GetDopplerYta(beta,yta,&track));
*/
/*
  if(mult == 1) {
    obj.FillHistogram(dirname,Form("Gamma_Singles%s",hsuf.c_str()),2000,0,4000,hit.GetDopplerYta(beta,yta,&track));
  }
  else if(mult == 2) {
    obj.FillHistogram(dirname,Form("Gamma_Mult2%s",hsuf.c_str()),2000,0,4000,hit.GetDopplerYta(beta,yta,&track));
  }
  else if(mult > 1) {
    obj.FillHistogram(dirname,Form("Gamma_Mult2Plus%s",hsuf.c_str()),2000,0,4000,hit.GetDopplerYta(beta,yta,&track));
  }
*/
//obj.FillHistogram(dirname,Form("Double_Resolution_Gamma(Beta&Track&YTA)%s",hsuf.c_str()),1000,0,4000,hit.GetDopplerYta(beta,yta,&track));
//obj.FillHistogram(dirname,Form("Half_Resolution_Gamma(Beta&Track&YTA)%s",hsuf.c_str()),250,0,4000,hit.GetDopplerYta(beta,yta,&track));
  obj.FillHistogram(dirname,Form("Gamma_v_Time%s",hsuf.c_str()),
		    250,-500,1000,bank29->Timestamp() - hit.GetTime(),
		    2000,0,16000,hit.GetDopplerYta(s800->AdjustedBeta(beta),yta,&adjusted_track));
  obj.FillHistogram(dirname,Form("Gamma_v_TSDiff%s",hsuf.c_str()),
		    250,-500,1000,bank29->Timestamp() - hit.Timestamp(),
		    2000,0,16000,hit.GetDopplerYta(s800->AdjustedBeta(beta),yta,&adjusted_track));
  obj.FillHistogram(dirname,Form("CoreEn_v_TSDiff%s",hsuf.c_str()),
		    250,-500,1000,bank29->Timestamp() - hit.Timestamp(),
		    2000,0,16000,hit.GetCoreEnergy());
  obj.FillHistogram(dirname,Form("CoreEn_v_Time%s",hsuf.c_str()),
		    250,-500,1000,bank29->Timestamp() - hit.GetTime(),
		    2000,0,16000,hit.GetCoreEnergy());

//  obj.FillHistogram(dirname,Form("GretinaCoreEnergy%s",hsuf.c_str()),500,0,4000,hit.GetCoreEnergy());



//This section of code is all for track-reconstructed corrections, so it is commented out until I figure out how to use that.
  //obj.FillHistogram(dirname,Form("Gamma(Beta&Track)%s",hsuf.c_str()),500,0,4000,hit.GetDoppler_dB(beta,&track));
  /*
  if(mult < 4) {
    obj.FillHistogram(dirname,Form("Gamma(Beta&Track)Mult%i%s",mult,hsuf.c_str()),
		      4000,0,4000,hit.GetDoppler(beta,&track));

    obj.FillHistogram(dirname,Form("Gamma(Beta&Track&YTA)Mult%i%s",mult,hsuf.c_str()),
		      4000,0,4000,hit.GetDoppler(beta,yta,&track));
  }

  else if(mult >= 4) {

    obj.FillHistogram(dirname,Form("Gamma(Beta&Track)Mult4+%s",hsuf.c_str()),
		      4000,0,4000,hit.GetDoppler(beta,&track));

    obj.FillHistogram(dirname,Form("Gamma(Beta&Track&YTA)Mult4+%s",hsuf.c_str()),
		      4000,0,4000,hit.GetDoppler(beta,yta,&track));
  }
  
  obj.FillHistogram(dirname,Form("Gamma(Beta&Track&YTA)%s",hsuf.c_str()),
		    4000,0,4000,hit.GetDoppler(beta,yta,&track));

  obj.FillHistogram(dirname,Form("Gamma(Beta&Track&YTA&DTA)%s",hsuf.c_str()),
                    4000,0,4000,hit.GetDoppler(s800->AdjustedBeta(beta),yta,&track));
  
  //obj.FillHistogram(dirname,Form("CoreEn_v_TDiff%s",hsuf.c_str()),
  //		    1200,-400,800,bank29->Timestamp() - hit.GetTime(),
  //		    2000,0,8000,hit.GetCoreEnergy());

  obj.FillHistogram(dirname,Form("CoreEn_v_TSDiff%s",hsuf.c_str()),
		    600,-400,800,bank29->Timestamp() - hit.Timestamp(),
		    2000,0,8000,hit.GetCoreEnergy());
  obj.FillHistogram(dirname,Form("Gamma(Beta&Adjusted_Track)_v_DetMap%s",hsuf.c_str()),
		    38,0,38,detMap[hit.GetCrystalId()],
		    1000,0,4000,hit.GetDopplerYta(s800->AdjustedBeta(beta),yta,&adjusted_track));
  obj.FillHistogram(dirname,Form("Gamma(Beta)_v_DetMap%s",hsuf.c_str()),
		    38,0,38,detMap[hit.GetCrystalId()],
		    1000,0,4000,hit.GetDopplerYta(s800->AdjustedBeta(beta),yta,&track));
*/
//obj.FillHistogram(dirname,Form("Hits/Crystal%s",hsuf.c_str()),
//                  38,0,38,detMap[hit.GetCrystalId()]);
//obj.FillHistogram(dirname,Form("Hits/CrystalID%s",hsuf.c_str()),
//                  100,0,100,hit.GetCrystalId());
/*
  obj.FillHistogram(dirname,Form("Gamma(Beta)_v_QuadMap%s",hsuf.c_str()),
		    10,0,10,quadMap[hit.GetHoleNumber()],
		    1000,0,4000,hit.GetDoppler(beta));

  obj.FillHistogram(dirname,Form("Gamma(Beta)_v_DetMap%s",hsuf.c_str()),
		    38,0,38,detMap[hit.GetCrystalId()],
		    1000,0,4000,hit.GetDoppler(beta));

  obj.FillHistogram(dirname,Form("CourseGamma(Beta)_v_DetMap%s",hsuf.c_str()),
		    38,0,38,detMap[hit.GetCrystalId()],
		    500,0,4000,hit.GetDoppler(beta));

  obj.FillHistogram(dirname,Form("CourseGamma(Beta)_v_QuadMap%s",hsuf.c_str()),
		    100,0,100,quadMap[hit.GetHoleNumber()],
		    500,0,4000,hit.GetDoppler(beta));
*/
/*
  if(quadMap[hit.GetHoleNumber()] <= 4) {
    obj.FillHistogram(dirname,Form("Gamma(Beta)_Forward_4keV%s",hsuf.c_str()),
      1000,0,4000,hit.GetDopplerYta(beta,yta,&track));
  }
  else {
    obj.FillHistogram(dirname,Form("Gamma(Beta)_90_4keV%s",hsuf.c_str()),
      1000,0,4000,hit.GetDopplerYta(beta,yta,&track));
  }
  if(quadMap[hit.GetHoleNumber()] <= 4) {
    obj.FillHistogram(dirname,Form("Gamma(Beta)_Forward_8keV%s",hsuf.c_str()),
      500,0,4000,hit.GetDopplerYta(beta,yta,&track));
  }
  else {
    obj.FillHistogram(dirname,Form("Gamma(Beta)_90_8keV%s",hsuf.c_str()),
      500,0,4000,hit.GetDopplerYta(beta,yta,&track));
  }
*/
  if(quadMap[hit.GetHoleNumber()] <= 4) {
    obj.FillHistogram(dirname,Form("Gamma(Beta)_Forward%s",hsuf.c_str()),
      2000,0,4000,hit.GetDopplerYta(s800->AdjustedBeta(beta),yta,&track));
  }
  else {
    obj.FillHistogram(dirname,Form("Gamma(Beta)_90%s",hsuf.c_str()),
      2000,0,4000,hit.GetDopplerYta(s800->AdjustedBeta(beta),yta,&track));
  }
/*  
  obj.FillHistogram(dirname,Form("Gamma(Beta&Track)_v_DetMap%s",hsuf.c_str()),
  		    38,0,38,detMap[hit.GetCrystalId()],
  		    4000,0,4000,hit.GetDoppler(beta,&track));

  obj.FillHistogram(dirname,Form("Gamma(Beta&Track&YTA)_v_DetMap%s",hsuf.c_str()),
  		    38,0,38,detMap[hit.GetCrystalId()],
  		    4000,0,4000,hit.GetDoppler(beta,yta,&track));

  obj.FillHistogram(dirname,Form("Gamma(Beta&Track&YTA&DTA)_v_DetMap%s",hsuf.c_str()),
  		    38,0,38,detMap[hit.GetCrystalId()],
  		    4000,0,4000,hit.GetDoppler(s800->AdjustedBeta(beta),yta,&track));
  */

  /*
  for(double b : BinCenters(80,0.38,0.42)) {
    obj.FillHistogram(dirname,Form("Gamma(Beta&Track)_BetaScan%s",hsuf.c_str()),
		      80,0.38,0.42,b,4000,0,4000,hit.GetDoppler(b,&track));
  }
  */
  /*
  for(double b : BinCenters(80,0.38,0.42)) {
    obj.FillHistogram(dirname,Form("Gamma(Beta&Track&YTA)_BetaScan%s",hsuf.c_str()),
		      80,0.38,0.42,b,4000,0,4000,hit.GetDoppler(b,&track));
  }
  */

  return true;
}

bool GammaCorrelations(TRuntimeObjects &obj, TGretinaHit hit, TS800* s800, std::string dsuf, std::string hsuf, double beta) {

  TVector3 track = s800->Track();
  double yta = s800->GetYta();
  double dta = s800->GetDta();

  double reac_phi = (track.Cross(TVector3(0.0,0.0,1.0))).Phi();
  if(reac_phi < 0)
    {reac_phi += TMath::TwoPi();}

  double det_phi = ((hit.GetPosition()).Cross(TVector3(0.0,0.0,1.0))).Phi();
  if(det_phi < 0)
    {det_phi += TMath::TwoPi();}

  double phi1 = reac_phi - det_phi;
  if(phi1 < 0)
    {phi1 += TMath::TwoPi();}

  double phi2 = TMath::TwoPi() - s800->Azita() - hit.GetPhi();
  if(phi2 < 0)
    {phi2 += TMath::TwoPi();}

  std::string dirname = "GammaCorrelations" + dsuf;

  obj.FillHistogram(dirname,Form("Phi2-Phi1%s",hsuf.c_str()),400,-6.3,6.3,phi2-phi1);

  obj.FillHistogram(dirname,Form("CoreEn_v_Phi1%s",hsuf.c_str()),
		    400,0,6.3,phi1,
		    4000,0,4000,hit.GetCoreEnergy());


  obj.FillHistogram(dirname,Form("Gamma(Beta)_v_Phi1%s",hsuf.c_str()),
		    400,0,6.3,phi1,
		    4000,0,4000,hit.GetDoppler(beta));

  obj.FillHistogram(dirname,Form("Gamma(Beta&Track)_v_Phi1%s",hsuf.c_str()),
		    400,0,6.3,phi1,
		    4000,0,4000,hit.GetDoppler(beta,&track));
  obj.FillHistogram(dirname,Form("CoreEn_v_Phi2%s",hsuf.c_str()),
		    400,0,6.3,phi2,
		    4000,0,4000,hit.GetCoreEnergy());

  obj.FillHistogram(dirname,Form("Gamma(Beta)_v_Phi2%s",hsuf.c_str()),
		    400,0,6.3,phi2,
		    4000,0,4000,hit.GetDoppler(beta));

  obj.FillHistogram(dirname,Form("Gamma(Beta&Track)_v_Phi2%s",hsuf.c_str()),
		    400,0,6.3,phi2,
		    4000,0,4000,hit.GetDoppler(beta,&track));

  double angle = (hit.GetPosition()).Angle(track)*TMath::RadToDeg();

  obj.FillHistogram(dirname,Form("CoreEn_v_Angle%s",hsuf.c_str()),
		    4000,0,4000,hit.GetCoreEnergy(),
		    180,0,180,angle);

  obj.FillHistogram(dirname,Form("Gamma(Beta)_v_Angle%s",hsuf.c_str()),
		    4000,0,4000,hit.GetDoppler(beta),
		    180,0,180,angle);

  obj.FillHistogram(dirname,Form("Gamma(Beta&Track)_v_Angle%s",hsuf.c_str()),
		    4000,0,4000,hit.GetDoppler(beta,&track),
		    180,0,180,angle);

  
/*
  dirname = "YTADetLevelCorrs" + dsuf;

  obj.FillHistogram(dirname,Form("Gamma(Beta&Track)_v_YTA_Det(T-P)%i%s",detMap[hit.GetCrystalId()],hsuf.c_str()),
		    60,-15,15,yta,
		    1000,0,4000,hit.GetDoppler(beta,&track));

  obj.FillHistogram(dirname,Form("Gamma(Beta&Track&YTA)_v_YTA_Det(T-P)%i%s",detMap[hit.GetCrystalId()],hsuf.c_str()),
		    60,-15,15,yta,
		    1000,0,4000,hit.GetDopplerYta(beta,yta,&track));
  dirname = "DTADetLevelCorrs" + dsuf;

  obj.FillHistogram(dirname,Form("Gamma(Beta&Track&YTA)_v_DTA_Det(T-P)%02i%s",detMap[hit.GetCrystalId()],hsuf.c_str()),
		    500,-0.5,0.5,dta,
		    1000,0,4000,hit.GetDopplerYta(beta,yta,&track));

  obj.FillHistogram(dirname,Form("Gamma(Beta&Track&YTA&DTA)_v_DTA_Det(T-P)%02i%s",detMap[hit.GetCrystalId()],hsuf.c_str()),
		    500,-0.5,0.5,dta,
		    1000,0,4000,hit.GetDopplerYta(s800->AdjustedBeta(beta),yta,&track));
*/  
  return true;
}

bool GammaGamma(TRuntimeObjects &obj, TGretinaHit hit1, TGretinaHit hit2, TS800* s800, std::string dsuf, std::string hsuf,
		bool ab, double beta) {

  std::string pref = "";
  if(ab)
    {pref = "AB_";}
  
  TVector3 track = s800->Track();
  double yta = s800->GetYta();
  double gretina_ata_offset = GValue::Value("GRETINA_ATA_OFFSET");
  double gretina_bta_offset = GValue::Value("GRETINA_BTA_OFFSET");
  TVector3 adjusted_track = s800->Track(gretina_ata_offset,gretina_bta_offset);


  std::string dirname = pref + "GammaGamma" + dsuf;
/*
  obj.FillHistogram(dirname,Form("8keVGammaGamma(Beta&Track&YTA)%s",hsuf.c_str()),
		    500,0,4000,hit1.GetDopplerYta(beta,yta,&track),
		    500,0,4000,hit2.GetDopplerYta(beta,yta,&track));
  obj.FillHistogram(dirname,Form("16keVGammaGamma(Beta&Track&YTA)%s",hsuf.c_str()),
		    250,0,4000,hit1.GetDopplerYta(beta,yta,&track),
		    250,0,4000,hit2.GetDopplerYta(beta,yta,&track));
*/
  obj.FillHistogram(dirname,Form("GammaGamma(Beta&Track&YTA)%s",hsuf.c_str()),
		    4000,0,4000,hit1.GetDopplerYta(s800->AdjustedBeta(beta),yta,&adjusted_track),
		    4000,0,4000,hit2.GetDopplerYta(s800->AdjustedBeta(beta),yta,&adjusted_track));

  
/*
  obj.FillHistogram(dirname,Form("GammaGamma(Beta&Track)%s",hsuf.c_str()),
		    2500,0,2500,hit1.GetDoppler(beta,&track),
		    2500,0,2500,hit2.GetDoppler(beta,&track));

  obj.FillHistogram(dirname,Form("GammaGamma(Beta&Track&YTA)%s",hsuf.c_str()),
		    2500,0,2500,hit1.GetDopplerYta(beta,yta,&track),
		    2500,0,2500,hit2.GetDopplerYta(beta,yta,&track));

  obj.FillHistogram(dirname,Form("GammaGamma(Beta&Track)_Exp%s",hsuf.c_str()),
		    2500,0,4000,hit1.GetDoppler(beta,&track),
		    2500,0,4000,hit2.GetDoppler(beta,&track));

  obj.FillHistogram(dirname,Form("GammaGamma(Beta&Track&YTA)_Exp%s",hsuf.c_str()),
		    2500,0,4000,hit1.GetDopplerYta(beta,yta,&track),
		    2500,0,4000,hit2.GetDopplerYta(beta,yta,&track));
*/
  return true;
}

bool GammaGammaGamma(TRuntimeObjects &obj, TGretinaHit hit1, TGretinaHit hit2, TGretinaHit hit3, TS800* s800, std::string dsuf, std::string hsuf,
		bool ab, double beta) {

  std::string pref = "";
  if(ab)
    {pref = "AB_";}
  
  TVector3 track = s800->Track();
  double yta = s800->GetYta();
  double gretina_ata_offset = GValue::Value("GRETINA_ATA_OFFSET");
  double gretina_bta_offset = GValue::Value("GRETINA_BTA_OFFSET");
  TVector3 adjusted_track = s800->Track(gretina_ata_offset,gretina_bta_offset);


  std::string dirname = pref + "GammaGammaGamma" + dsuf;

  if((hit1.GetDopplerYta(s800->AdjustedBeta(beta),yta,&adjusted_track) > 1510) && (hit1.GetDopplerYta(s800->AdjustedBeta(beta),yta,&adjusted_track) < 1530)) {
    obj.FillHistogram(dirname,Form("1520GammaGamma(Beta&Track&YTA)%s",hsuf.c_str()),
		    4000,0,4000,hit2.GetDopplerYta(s800->AdjustedBeta(beta),yta,&adjusted_track),
		    4000,0,4000,hit3.GetDopplerYta(s800->AdjustedBeta(beta),yta,&adjusted_track));
  }
  if((hit1.GetDopplerYta(s800->AdjustedBeta(beta),yta,&adjusted_track) > 1565) && (hit1.GetDopplerYta(s800->AdjustedBeta(beta),yta,&adjusted_track) < 1585)) {
    obj.FillHistogram(dirname,Form("1575GammaGamma(Beta&Track&YTA)%s",hsuf.c_str()),
		    4000,0,4000,hit2.GetDopplerYta(s800->AdjustedBeta(beta),yta,&adjusted_track),
		    4000,0,4000,hit3.GetDopplerYta(s800->AdjustedBeta(beta),yta,&adjusted_track));
  }
  if((hit1.GetDopplerYta(s800->AdjustedBeta(beta),yta,&adjusted_track) > 1110) && (hit1.GetDopplerYta(s800->AdjustedBeta(beta),yta,&adjusted_track) < 1130)) {
    obj.FillHistogram(dirname,Form("1120GammaGamma(Beta&Track&YTA)%s",hsuf.c_str()),
		    4000,0,4000,hit2.GetDopplerYta(s800->AdjustedBeta(beta),yta,&adjusted_track),
		    4000,0,4000,hit3.GetDopplerYta(s800->AdjustedBeta(beta),yta,&adjusted_track));
  }
  if((hit1.GetDopplerYta(s800->AdjustedBeta(beta),yta,&adjusted_track) > 375) && (hit1.GetDopplerYta(s800->AdjustedBeta(beta),yta,&adjusted_track) < 385)) {
    obj.FillHistogram(dirname,Form("380GammaGamma(Beta&Track&YTA)%s",hsuf.c_str()),
		    4000,0,4000,hit2.GetDopplerYta(s800->AdjustedBeta(beta),yta,&adjusted_track),
		    4000,0,4000,hit3.GetDopplerYta(s800->AdjustedBeta(beta),yta,&adjusted_track));
  }
  if((hit1.GetDopplerYta(s800->AdjustedBeta(beta),yta,&adjusted_track) > 300) && (hit1.GetDopplerYta(s800->AdjustedBeta(beta),yta,&adjusted_track) < 385)) {
    obj.FillHistogram(dirname,Form("380TailGammaGamma(Beta&Track&YTA)%s",hsuf.c_str()),
		    4000,0,4000,hit2.GetDopplerYta(s800->AdjustedBeta(beta),yta,&adjusted_track),
		    4000,0,4000,hit3.GetDopplerYta(s800->AdjustedBeta(beta),yta,&adjusted_track));
  }
  if((hit1.GetDopplerYta(s800->AdjustedBeta(beta),yta,&adjusted_track) > 350) && (hit1.GetDopplerYta(s800->AdjustedBeta(beta),yta,&adjusted_track) < 360)) {
    obj.FillHistogram(dirname,Form("355GammaGamma(Beta&Track&YTA)%s",hsuf.c_str()),
		    4000,0,4000,hit2.GetDopplerYta(s800->AdjustedBeta(beta),yta,&adjusted_track),
		    4000,0,4000,hit3.GetDopplerYta(s800->AdjustedBeta(beta),yta,&adjusted_track));
  }

  return true;
}

bool HandleGretina_Gated(TRuntimeObjects &obj,TCutG *incoming, TCutG* outgoing, TCutG* time_energy, TCutG* xToF) {
  TS800 *s800  = obj.GetDetector<TS800>();
  TGretina *gretina  = obj.GetDetector<TGretina>();
  TBank29 *bank29 = obj.GetDetector<TBank29>();

  //if(!s800 || !gretina || !bank29 || !outgoing || ! incoming || !time_energy || !xToF)
//This line ensures that we are only executing this if we have a gate for each field (currently incoming and outgoing)
  if(!s800 || !gretina || !bank29 || !outgoing || !incoming || !time_energy)
    {return false;}

  std::string d_suf = "";
  
  //////////////////////////////////////////////////////////
  //Incoming  Gate on Xfp - E1 vs. Object - E1 Uncorrected//
  //////////////////////////////////////////////////////////
//IsInside tests whether or not the event is inside the gate
  if(incoming) {  
    if(!incoming->IsInside(s800->GetMTof().GetCorrelatedObjE1(),s800->GetMTof().GetCorrelatedXfpE1()))
      {return false;}
    d_suf += Form("_%s",incoming->GetName());
  }

  std::string h_suf  = "";

  //////////////////////////////////////////////////////////////////
  //Outgoing Gate on Ion Chamber Charge vs Object - E1 (Corrected)//
  //////////////////////////////////////////////////////////////////
  if(outgoing) {
    //if(!outgoing->IsInside(s800->GetOBJRaw_TAC(),s800->GetIonChamber().Charge()))
    if(!outgoing->IsInside(s800->GetMTofObjE1(),s800->GetIonChamber().Charge()))
      {return false;}
    h_suf += Form("_%s",outgoing->GetName());
  }

  /////////////////////////////////////////////////////////
  //Correlation Gate on CRDC1X vs Object - E1 (Corrected)//
  /////////////////////////////////////////////////////////
  //if(xToF) {
  //  if(!xToF->IsInside(s800->GetMTofXfpE1(),s800->GetXFP()))
  //    {return false;}
  //  h_suf += Form("_%s",xToF->GetName());
  //}

  int size = 0;
  if(time_energy) {
    d_suf+=Form("_%s",time_energy->GetName());
  }

//There can be multiple detections associated with a single event, so we will loop over those here
  TGretina* good_gret = new TGretina(); 
  for(unsigned int i=0;i<gretina->Size();i++) {
    TGretinaHit hit = gretina->GetGretinaHit(i);

    //////////////////////////////////////////////////////////////////////////////
    //Time-Energy Gate on Gretina Core Enery vs Bank29 TimeStamp - Hit Timestamp//
    //////////////////////////////////////////////////////////////////////////////
    if(time_energy) {
      if(!time_energy->IsInside(bank29->Timestamp() - hit.GetTime(),hit.GetCoreEnergy()))
        {continue;}  
    }
    size++;
    good_gret->InsertHit(hit);
  } //end gretina hit loop

//Gated_Gretina_Spectra does histograms for gated gamma detection data, but we filter out energies less than 100keV here
  std::string dirname = "Gretina" + d_suf;
  obj.FillHistogram(dirname,Form("GoodTE_Mult%s",h_suf.c_str()),25,0,25,size);
  obj.FillHistogram(dirname,Form("Event_Mult%s",h_suf.c_str()),25,0,25,gretina->Size());
  for(unsigned int i=0;i<good_gret->Size();i++) {

    TGretinaHit hit = good_gret->GetGretinaHit(i);
//    if(hit.GetCoreEnergy() > 100) {
      double beta = GValue::Value(Form("%s_BETA",outgoing->GetName()));
      Gated_Gretina_Spectra(obj,hit,s800,bank29,size,d_suf,h_suf,false,beta);
      GammaCorrelations(obj,hit,s800,d_suf,h_suf,beta);
//    }
    for(unsigned int j=0;j<good_gret->Size();j++) {
      if(i==j) //want different hits
	{continue;}
      
      TGretinaHit hit2 = good_gret->GetGretinaHit(j);
//      if(hit.GetCoreEnergy() > 100 && hit2.GetCoreEnergy() > 100) {
        double beta = GValue::Value(Form("%s_BETA",outgoing->GetName()));
        GammaGamma(obj,hit,hit2,s800,d_suf,h_suf,false,beta);
//      }
    } //end second good gretina hit loop
  } //end first good gretina loop

//I don't know what this does yet
  for(int i=0;i<good_gret->AddbackSize();i++) {
    TGretinaHit hit = good_gret->GetAddbackHit(i);

//    if(hit.GetCoreEnergy() > 100) {
      double beta = GValue::Value(Form("%s_BETA",outgoing->GetName()));
      Gated_Gretina_Spectra(obj,hit,s800,bank29,good_gret->AddbackSize(),d_suf,h_suf,true,beta);
//    }
    for(int j=0;j<good_gret->AddbackSize();j++) {
      if(i==j) //want different hits
	{continue;}

      TGretinaHit hit2 = good_gret->GetAddbackHit(j);
//      if(hit.GetCoreEnergy() > 100 && hit2.GetCoreEnergy() > 100) {
        double beta = GValue::Value(Form("%s_BETA",outgoing->GetName()));
        GammaGamma(obj,hit,hit2,s800,d_suf,h_suf,true,beta);
//      }
       for(int k=0;k<good_gret->AddbackSize();k++) {
         if((k==i) || (k==j))
           {continue;}
         TGretinaHit hit3 = good_gret->GetAddbackHit(k);
          double beta = GValue::Value(Form("%s_BETA",outgoing->GetName()));
          GammaGammaGamma(obj,hit,hit2,hit3,s800,d_suf,h_suf,true,beta);
       }
         
      
    } //end second gretina addback hit loop
  } //end gretina addback hit loop

  dirname = "AB_Gretina" + d_suf;
  obj.FillHistogram(dirname,Form("GoodTE_ABMult%s",h_suf.c_str()),25,0,25,good_gret->AddbackSize());
  return true;
}

std::vector<GCutG*> incoming_cuts = {0};
std::vector<GCutG*> outgoing_cuts = {0};
std::vector<GCutG*> time_energy_cuts = {0};
std::vector<GCutG*> xToF_cuts = {0};

/*
GCutG* incoming_cut = 0;
GCutG* outgoing_cut = 0;
GCutG* time_energy_cut = 0;
GCutG* xToF_cut = 0;
*/

int gates_loaded=0;

// extern "C" is needed to prevent name mangling.
// The function signature must be exactly as shown here,
//   or else bad things will happen.
extern "C"
//Make Histograms is the 'main' function that actually gets executed and calls all of the other code in here. It must always take TRuntimeObjects as its argumnet.
void MakeHistograms(TRuntimeObjects& obj) {

  TList *list = &(obj.GetObjects());
  int numobj = list->GetSize();  
   
//Read in gates. Gates are created graphically in the GUI. Tags can be set to anything as long as you are consistant. Here I have tags for incoming (Tof v XTof), outgoing (Tof v Charge), time energy (Don't know, something with gammas), and x_tof (x_tof v charge)
  TList *gates = &(obj.GetGates());
  
  if(gates_loaded!=gates->GetSize()) {
    TIter iter(gates);
    while(TObject *obj = iter.Next()) {
      GCutG *gate = (GCutG*)obj;
      std::string tag = gate->GetTag();
      if(!tag.compare("incoming")) {
        incoming_cuts.push_back(gate);
        std::cout << "incoming: << " << gate->GetName() << std::endl;
      } else if(!tag.compare("outgoing")) {
        outgoing_cuts.push_back(gate);
        std::cout << "outgoing: << " << gate->GetName() << std::endl;
      } else if(!tag.compare("time_energy")) {
        time_energy_cuts.push_back(gate);
        std::cout << "time energy: << " << gate->GetName() << std::endl;
//      } else if(!tag.compare("xtof")) {
//        xToF_cuts.push_back(gate);
//        std::cout << "xToF: << " << gate->GetName() << std::endl;
      }
      gates_loaded++;
    }
  }
  
  /*
  if(gates_loaded!=gates->GetSize()) {
    TIter iter(gates);
    while(TObject *obj = iter.Next()) {
      GCutG *gate = (GCutG*)obj;
      std::string tag = gate->GetTag();
      if(!tag.compare("incoming")) {
        incoming_cut = gate;
        std::cout << "incoming: << " << gate->GetName() << std::endl;
      } else if(!tag.compare("outgoing")) {
        outgoing_cut = gate;
        std::cout << "outgoing: << " << gate->GetName() << std::endl;
      } else if(!tag.compare("time_energy")) {
        time_energy_cut = gate;
        std::cout << "time energy: << " << gate->GetName() << std::endl;
      } else if(!tag.compare("x_tof")) {
        xToF_cut = gate;
        std::cout << "xToF: << " << gate->GetName() << std::endl;
      }
      gates_loaded++;
    }
  }
  
  
  HandleGretina(obj);
  HandleGretina_Gated(obj,incoming_cut,outgoing_cut,time_energy_cut,xToF_cut);
  */
  
//Histogram overall results from Gretina, without filters.
  HandleGretina(obj);
  
//Loop over all combinations of incoming, outgoing, and xTof gates and histogram S800 results
  for(size_t i=0;i<incoming_cuts.size();i++) {
   for(size_t j=0;j<outgoing_cuts.size();j++) {
    for(size_t l=0;l<xToF_cuts.size();l++) {
      //Gated_CS_Spectra(obj,incoming_cuts.at(i),outgoing_cuts.at(j),xToF_cuts.at(l));
//    HandleTiming_Gated(obj,incoming_cuts.at(i),outgoing_cuts.at(j),xToF_cuts.at(l));
      HandleS800_Gated(obj,incoming_cuts.at(i),outgoing_cuts.at(j),xToF_cuts.at(l));
//Loop over all time energy gates (within the loops for the other three types) and histogram Gretina results
     HandleS800_Gated(obj,incoming_cuts.at(i),outgoing_cuts.at(j),xToF_cuts.at(l));
     for(size_t k=0;k<time_energy_cuts.size();k++) {
       HandleGretina_Gated(obj,incoming_cuts.at(i),outgoing_cuts.at(j),time_energy_cuts.at(k),xToF_cuts.at(l));
     } //end time_energy gate loop
    } //end x_tof gate loop
   } //end outgoing gate loop
  } //end incoming gate loop
  if(numobj!=list->GetSize())
    list->Sort();
}
