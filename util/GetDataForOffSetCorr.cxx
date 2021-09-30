#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include "TFile.h"
#include "TTree.h"

#include "GValue.h"
#include "TChannel.h"
#include "TGretinaHit.h"
#include "TGretina.h"
#include "TS800.h"
#include "TBank29.h"
#include "TInverseMap.h"
#include "GCutG.h"



//Loop through list_of_files and obtain core energy, position in lab frame, and
//yta for each event. This will be used to determine beam spot/offset
//correction. File names in list_of_files should correspond to 
void GetDataForOffsetCorr(std::vector<std::string> list_of_files,  
    std::string gval_file_prefix, 
    std::string main_gval_fn, std::string inv_map_fn, GCutG *in_cut, GCutG *out_cut, GCutG *prompt_cut, GCutG *x_tof_cut = 0){
  TS800 *s800 = 0;
  TGretina *gret = 0;
  TBank29 *bank29 = 0;

  std::ofstream out_file;
  out_file.open("output_for_shifts.dat");

  std::map<int,int> detMap = {{26,0}, {30,1}, {34,2}, {38,3}, {25,4}, {29,5}, {33,6}, {37,7}, {27,8}, {31,9}, 
                              {35,10}, {39,11}, {24,12}, {28,13}, {32,14}, {36,15}, {47,16}, {51,17}, {59,18}, 
                              {63,19}, {67,20}, {71,21}, {79,22}, {44,23}, {50,24}, {58,25}, {60,26}, {66,27}, 
                              {68,28}, {76,29}, {46,30}, {48,31}, {56,32}, {62,33}, {64,34}, {70,35}, {78,36}, 
                              {45,37}, {49,38}, {57,39}, {61,40}, {65,41}, {69,42}, {77,43}};

  if(GValue::ReadValFile(main_gval_fn.c_str()) <= 0){
    std::cout << "Failed to open GValue overall file.\n";
    return;
  }

  TInverseMap::Get(inv_map_fn.c_str());

  for (auto &file_name : list_of_files){
    TFile cur_file(file_name.c_str(), "read");
    //Have to get all corrections and cuts loaded as well!
    std::string val_fn(gval_file_prefix+file_name.replace(8,4,"val"));
    if (GValue::ReadValFile(val_fn.c_str()) <= 0){
      std::cout << "Failed to open run-by-run GValue file.\n";
      return;
    }
  

    TTree *tr = (TTree*)cur_file.Get("EventTree");
    tr->SetBranchAddress("TS800", &s800);
    tr->SetBranchAddress("TGretina", &gret);
    tr->SetBranchAddress("TBank29", &bank29);

    for (int i = 0; i < tr->GetEntries(); i++){
      tr->GetEntry(i);

      if (!s800 || !gret || !bank29){
        continue;
      }

      double crdc_1_x = s800->GetCrdc(0).GetDispersiveX();
      double crdc_1_y = s800->GetCrdc(0).GetNonDispersiveY();
      double corr_obj = s800->GetMTofObjE1();
      //double ic_de = s800->GetIonChamber().GetdE(crdc_1_x,crdc_1_y);
      double ic_de = s800->GetIonChamber().Charge();
      //Check if incoming, outgoing satisfied
      if(!s800->GetMTof().ObjSize() || !s800->GetMTof().XfpSize() || !s800->GetMTof().E1UpSize()){
        continue;
      } else if(!in_cut->IsInside(s800->GetMTof().GetCorrelatedObjE1(), 
            s800->GetMTof().GetCorrelatedXfpE1())){
        continue;
      } 

      if (out_cut && !out_cut->IsInside(corr_obj,ic_de)){
          continue; 
      }

      if (x_tof_cut && !x_tof_cut->IsInside(corr_obj, crdc_1_x)){
        continue;
      }

      TVector3 track = s800->Track();
      double yta = s800->GetYta();
      double ata = s800->GetAta();
      double bta = s800->GetBta();
      double adjusted_beta = s800->AdjustedBeta(GValue::Value("p38_BETA"));
      for (unsigned int j = 0; j < gret->Size(); j++){
        //Check if doppler corrected energy is in correct range (1220, 1300)
        const TGretinaHit &hit = gret->GetGretinaHit(j);

        double energy = hit.GetDoppler(GValue::Value("p38_BETA"));
        double time = bank29->Timestamp()-hit.GetTime();
        if (prompt_cut->IsInside(time, energy)){
          double energy_track_yta_dta = hit.GetDopplerYta(adjusted_beta, yta, 0.0,0.0,0.0, &track);
          if (energy_track_yta_dta > 514 && energy_track_yta_dta < 526){//corrected energy in correct range 
            int num = detMap[hit.GetCrystalId()];
            TVector3 pos = hit.GetPosition();
            double core_energy = hit.GetCoreEnergy();

            //dump to file
            out_file << num << " "  << pos.X() << " " << pos.Y() << " " << pos.Z() << " "; 
            out_file << core_energy << " " << yta << " " << adjusted_beta << " "; 
            out_file << ata << " " << bta << "\n"; 
          }
        }
      }
      s800->Clear();
      gret->Clear();
      bank29->Clear();
    }//loop over tree entries

    delete tr;
  }//loop over files 
  out_file.close();
}
void runAr46(){
  std::vector<std::string> file_names = {
    "run0034.root",
    "run0035.root",
    "run0036.root",
    "run0037.root",
    "run0038.root",
    "run0039.root",
    "run0040.root",
    "run0044.root",
    "run0045.root",
    "run0046.root",
    "run0047.root",
    "run0048.root",
    "run0049.root",
    "run0050.root",
    "run0051.root",
    "run0052.root",
    "run0053.root",
    "run0054.root",
    "run0057.root",
    "run0058.root",
    "run0059.root",
    "run0060.root",
    "run0061.root",
    "run0062.root",
    "run0063.root",
    "run0064.root",
    "run0065.root",
    "run0066.root",
    "run0067.root",
    "run0068.root",
    "run0069.root",
    "run0070.root",
    "run0077.root",
    "run0078.root",
    "run0079.root",
    "run0081.root",
    "run0082.root",
    "run0083.root",
    "run0087.root",
    "run0088.root",
    "run0089.root",
    "run0090.root",
    "run0091.root",
    "run0092.root",
    "run0093.root",
    "run0094.root",
    "run0095.root",
    "run0096.root",
    "run0097.root",
    "run0098.root",
    "run0099.root",
    "run0100.root",
    "run0101.root",
    "run0102.root",
    "run0103.root",
    "run0104.root",
    "run0105.root",
    "run0106.root",
    "run0107.root",
    "run0108.root",
    "run0109.root",
    "run0110.root",
    "run0111.root",
    "run0135.root",
    "run0136.root",
    "run0137.root",
    "run0138.root",
    "run0139.root",
    "run0140.root",
    "run0141.root",
    "run0142.root",
    "run0143.root",
  };


  std::string gval_file_prefix("/user/hillm/analysis/e19002/config/");
  std::string main_gval_fn("/user/hillm/analysis/e19002/config/values.val");
  std::string inv_map_fn("/user/hillm/analysis/e19002/inverseMaps/in46Ar_out39P.inv");
  GCutG *in_cut = new GCutG("incoming",9);
  in_cut->SetPoint(0,-2852.96,1219.3);
  in_cut->SetPoint(1,-2677.15,1495.26);
  in_cut->SetPoint(2,-2595.15,1591.08);
  in_cut->SetPoint(3,-2517.58,1656.24);
  in_cut->SetPoint(4,-2262.46,1857.46);
  in_cut->SetPoint(5,-2398.14,1644.74);
  in_cut->SetPoint(6,-2478.4,1554.67);
  in_cut->SetPoint(7,-2646.57,1384.11);
  in_cut->SetPoint(8,-2852.96,1219.3);
  GCutG *out_cut = new GCutG("p38",9);
  out_cut->SetPoint(0,-2551.72,23679.1);
  out_cut->SetPoint(1,-2563.71,24731.6);
  out_cut->SetPoint(2,-2555.07,26757.4);
  out_cut->SetPoint(3,-2529.9,26902.1);
  out_cut->SetPoint(4,-2515.75,26665.3);
  out_cut->SetPoint(5,-2508.08,25290.6);
  out_cut->SetPoint(6,-2508.08,23679.1);
  out_cut->SetPoint(7,-2534.46,23416);
  out_cut->SetPoint(8,-2551.72,23679.1);
  GCutG *prompt_cut = new GCutG("timeEnergy",9);
  prompt_cut->SetPoint(0,47.4586,-118.467);
  prompt_cut->SetPoint(1,48.7724,819.432);
  prompt_cut->SetPoint(2,68.5113,6885.02);
  prompt_cut->SetPoint(3,69.2623,15909.4);
  prompt_cut->SetPoint(4,103.098,15874.6);
  prompt_cut->SetPoint(5,101.584,8609.76);
  prompt_cut->SetPoint(6,100.09,-362.369);
  prompt_cut->SetPoint(7,79.0376,-118.467);
  prompt_cut->SetPoint(8,47.4586,-118.467);

  GetDataForOffsetCorr(file_names, gval_file_prefix, 
      main_gval_fn, inv_map_fn, in_cut, out_cut, prompt_cut);
}

int main(){
//std::cout << "Running Cu71.." << std::endl;
//runCu71();
//std::cout << "Running Cu71 First Set.." << std::endl;
//runCu71_firstset();
  runAr46();
//std::cout << "Running Ni71.." << std::endl;
//runNi71();
}
