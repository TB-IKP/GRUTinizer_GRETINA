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

//------------------------------------------------------//
//	Bin stuff					//
//------------------------------------------------------//

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

//------------------------------------------------------//
//	Detector positions?				//
//------------------------------------------------------//

std::map<int,int> detMap = {{26,1},{30,2},{34,3},{38,4},{25,5},{29,6},{33,7},{37,8},{27,9},{31,10},{35,11},{39,12},
				{24,13},{28,14},{32,15},{36,16},{63,17},{71,18},{79,19},{59,20},{67,21},{58,22},{66,23},{60,24},
				{68,25},{76,26},{62,27},{70,28},{78,29},{56,30},{64,31},{57,32},{65,33},{61,34},
				{69,35},{77,36}};

std::map<int,int> quadMap = {{5,1},{6,2},{7,3},{8,4},{14,5},{16,6},{18,7},{13,8},{15,9}};

std::map<int,int> crysThetaMap = {{26,1},{30,1},{34,1},{38,1},{25,2},{29,2},{33,2},{37,2},{27,3},{31,3},{35,3},{39,3},
				{24,4},{28,4},{32,4},{36,4},{63,5},{71,5},{79,5},{59,6},{67,6},{58,7},{66,7},{60,8},
				{68,8},{76,8},{62,9},{70,9},{78,9},{56,10},{64,10},{57,11},{65,11},{61,12},{69,12},
				{77,12}};

//------------------------------------------------------//
//	Timing						//
//------------------------------------------------------//

bool HandleTiming_Gated(TRuntimeObjects &obj, TCutG *incoming, TCutG* outgoing, TCutG* xToF) {
	TS800 *s800  = obj.GetDetector<TS800>();

	if(!s800)
		{return false;}

	std::string dirname = "Timing";

	//////////////////////////////////////////////////////////
	//incoming gates on Xfp - E1 vs. Object - E1 Uncorrected//
	//////////////////////////////////////////////////////////
	if(incoming) {
		if(!incoming->IsInside(s800->GetMTof().GetCorrelatedObjE1(),s800->GetMTof().GetCorrelatedXfpE1()))
			{return false;}
		dirname += Form("_%s",incoming->GetName());
	}

	std::string suf = "";

	///////////////////////////////////////////////////////////////////
	//Outgoing Gates on Ion Chamber Charge vs Object - E1 (Corrected)//
	///////////////////////////////////////////////////////////////////
	if(outgoing) {
		if(!outgoing->IsInside(s800->GetMTofObjE1(),s800->GetIonChamber().Charge()))
			{return false;}
		suf += Form("_%s",outgoing->GetName());
	}

	//////////////////////////////////////////////////////////
	//Correlation Gates on CRDC1X vs Object - E1 (Corrected)//
	//////////////////////////////////////////////////////////
	if(xToF) {
		if(!xToF->IsInside(s800->GetMTofObjE1(),s800->GetXFP()))
			{return false;}
		suf += Form("_%s",xToF->GetName());
	}

	//Raw time histograms (Xfp, Obj, and E1)

	obj.FillHistogram(dirname,Form("E1_raw%s",suf.c_str()),
		4000,-1000,1000,s800->GetTof().GetRF());
	obj.FillHistogram(dirname,Form("Xfp_raw%s",suf.c_str()),
		4000,-1000,1000,s800->GetTof().GetXFP());
	obj.FillHistogram(dirname,Form("Obj_raw%s",suf.c_str()),
		4000,-1000,1000,s800->GetTof().GetOBJ());

	//Histogram Tof Xfp-E1

	obj.FillHistogram(dirname,Form("Xfp-E1%s",suf.c_str()),
		4000,-500000,500000,s800->GetMTof().GetCorrelatedXfpE1());

	//Histogram Tof Obj-E1

	obj.FillHistogram(dirname,Form("Obj-E1%s",suf.c_str()),
		4000,-10000,10000,s800->GetMTof().GetCorrelatedObjE1());

	return true;

}

//------------------------------------------------------//
//	S800						//
//------------------------------------------------------//

bool HandleS800_Gated(TRuntimeObjects &obj,TCutG *incoming, TCutG* outgoing, TCutG* xToF) {
	TS800 *s800  = obj.GetDetector<TS800>();

	if(!s800)
		{return false;}

	std::string dsuf = "";

	///////////////////////////////////////////////////////////
	//Incoming  Gates on Xfp - E1 vs. Object - E1 Uncorrected//
	///////////////////////////////////////////////////////////

	if(incoming) {
		if(!incoming->IsInside(s800->GetMTofObjE1(),s800->GetMTofXfpE1()))
			{return false;}
		dsuf += Form("_%s",incoming->GetName());
	}

	std::string hsuf = "";

	///////////////////////////////////////////////////////////////////
	//Outgoing Gates on Ion Chamber Charge vs Object - E1 (Corrected)//
	///////////////////////////////////////////////////////////////////

	if(outgoing) {
		if(!outgoing->IsInside(s800->GetMTofObjE1(),s800->GetIonChamber().Charge()))
			{return false;}
		hsuf += Form("_%s",outgoing->GetName());
	}

	std::string dirname = "";

	//Incoming Beam and Outgoing PIDs
	dirname = "S800" + dsuf;

	//Histogram Charge vs Tof

	obj.FillHistogram(dirname,Form("Charge_v_Obj-E1Raw%s",hsuf.c_str()),
		1000,-5000,5000,s800->GetMTof().GetCorrelatedObjE1(),
		1000,0,40000,s800->GetIonChamber().Charge());
	obj.FillHistogram(dirname,Form("Charge_v_Obj-E1%s",hsuf.c_str()),
		160,-2400,-1600,s800->GetMTofObjE1(),
		4000,0,40000,s800->GetIonChamber().Charge());

	//Histogram Focal Plane vs Tof

	obj.FillHistogram(dirname,Form("XFP_v_Obj-E1_Uncorrected%s",hsuf.c_str()),
		600,-5000,0,s800->GetMTof().GetCorrelatedObjE1(),
		1000,-500,500,s800->GetXFP());
	obj.FillHistogram(dirname,Form("XFP_v_Obj-E1%s",hsuf.c_str()),
		600,-5000,0,s800->GetMTofObjE1(),
		1000,-500,500,s800->GetXFP());
	obj.FillHistogram(dirname,Form("AFP_v_Obj-E1_Uncorrected%s",hsuf.c_str()),
		400,-5000,0,s800->GetMTof().GetCorrelatedObjE1(),
		1000,-0.05,0.05,s800->GetAFP());
	obj.FillHistogram(dirname,Form("AFP_v_Obj-E1%s",hsuf.c_str()),
		400,-5000,0,s800->GetMTofObjE1(),
		1000,-0.05,0.05,s800->GetAFP());

	return true;

}

//------------------------------------------------------//
//	GRETINA	ungated					//
//------------------------------------------------------//

bool HandleGretina(TRuntimeObjects &obj) {
	TGretina *gretina  = obj.GetDetector<TGretina>(); 

	if(!gretina)
		{return false;}

	std::string dirname = "Gretina_Ungated";

	obj.FillHistogram(dirname,"Multiplicity",25,0,25,gretina->Size());

	for(unsigned int i=0;i<gretina->Size();i++) {
		TGretinaHit hit = gretina->GetGretinaHit(i);

		obj.FillHistogram(dirname,"GretinaCoreEnergy_NoThresh",
				2000,0,4000,hit.GetCoreEnergy());

		if(hit.GetCoreEnergy() > 100) {

			obj.FillHistogram(dirname,"GretinaCoreEnergy",
				2000,0,4000,hit.GetCoreEnergy());
			obj.FillHistogram(dirname,"GretinaSummarySpectrum",
				200,0,200,hit.GetCrystalId(),
				2000,0,4000,hit.GetCoreEnergy());

		} //end if(hit.GetCoreEnergy() > 100)
	} //end gretina hit loop

	return true;
}

//------------------------------------------------------//
//	GRETINA	gated					//
//------------------------------------------------------//
/*
bool Gated_Gretina_Spectra(TRuntimeObjects &obj, TGretinaHit hit, TS800* s800, TBank29* bank29, int mult, 
		std::string dsuf,std::string hsuf, bool ab, double beta) {

	std::string pref = "";
	
	if(ab)
		{pref = "AB_";}

	std::string dirname = pref + "Gretina" + dsuf;

	obj.FillHistogram(dirname,Form("Gamma_Raw%s",hsuf.c_str()),500,0,4000,hit.GetCoreEnergy());

	return true;
}
*/
bool HandleGretina_Gated(TRuntimeObjects &obj,TCutG *incoming, TCutG* outgoing, TCutG* time_energy, TCutG* xToF) {
	TS800 *s800  = obj.GetDetector<TS800>();
	TGretina *gretina  = obj.GetDetector<TGretina>();

	//Make sure the event has data on S800 and Gretina
	if(!s800 || !gretina)
		{return false;}

	std::string d_suf = "";
	std::string h_suf = "";

	//////////////////////////////////////////////////////////////////
	//Outgoing Gate on Ion Chamber Charge vs Object - E1 (Corrected)//
	//////////////////////////////////////////////////////////////////

	if(outgoing) {
		if(!outgoing->IsInside(s800->GetMTofObjE1(),s800->GetIonChamber().Charge()))
			{return false;}
		d_suf += Form("_%s",outgoing->GetName());
	}

	std::string dirname = "Gretina_Gated" + d_suf;

	obj.FillHistogram(dirname,"Multiplicity",25,0,25,gretina->Size());

	for(unsigned int i=0;i<gretina->Size();i++) {
		TGretinaHit hit = gretina->GetGretinaHit(i);

		obj.FillHistogram(dirname,"GretinaCoreEnergy_NoThresh",
				2000,0,4000,hit.GetCoreEnergy());

		if(hit.GetCoreEnergy() > 100) {

			obj.FillHistogram(dirname,"GretinaCoreEnergy",
				2000,0,4000,hit.GetCoreEnergy());
			obj.FillHistogram(dirname,"GretinaSummarySpectrum",
				200,0,200,hit.GetCrystalId(),
				2000,0,4000,hit.GetCoreEnergy());

		} //end if(hit.GetCoreEnergy() > 100)
	} //end gretina hit loop

	return true;
}

//------------------------------------------------------//
//	???						//
//------------------------------------------------------//

std::vector<GCutG*> incoming_cuts = {0};
std::vector<GCutG*> outgoing_cuts = {0};
std::vector<GCutG*> time_energy_cuts = {0};
std::vector<GCutG*> xToF_cuts = {0};

int gates_loaded=0;

//extern "C" is needed to prevent name mangling.
//The function signature must be exactly as shown here,
//or else bad things will happen.
extern "C"
//Make Histograms is the 'main' function that actually gets executed and c
//alls all of the other code in here. It must always take TRuntimeObjects as its argumnet.
void MakeHistograms(TRuntimeObjects& obj) {

	TList *list = &(obj.GetObjects());
	int numobj = list->GetSize();  

	//Read in gates. Gates are created graphically in the GUI. 
	//Tags can be set to anything as long as you are consistant. 
	//Here I have tags for incoming (Tof v XTof), outgoing (Tof v Charge), 
	//time energy (Don't know, something with gammas), and x_tof (x_tof v charge)
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
			}

			gates_loaded++;
		}
	}

	//Histogram overall results from Gretina, without filters.
	HandleGretina(obj);

	//Loop over all combinations of incoming, outgoing, and xTof gates and histogram S800 results
	for(size_t i=0;i<incoming_cuts.size();i++) {

		for(size_t j=0;j<outgoing_cuts.size();j++) {

			for(size_t l=0;l<xToF_cuts.size();l++) {

				//Loop over all time energy gates (within the loops for the other three types) 
				//and histogram Gretina results
				HandleTiming_Gated(obj,incoming_cuts.at(i),outgoing_cuts.at(j),xToF_cuts.at(l));
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

