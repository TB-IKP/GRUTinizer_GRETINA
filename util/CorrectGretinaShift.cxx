#include <string>
#include <iostream>
#include <fstream>
#include <limits>


#include "TVector3.h"
#include "TMath.h"
//Run GetDataForOffsetCorr to get proper input file for this.
//The format is:
//  DirkOrderDetNum GretPosX GretPostY GretPosZ CoreEnergy YTA DTAAdjustedBeta ATA BTA 

struct Event {
  //detector number ordered from lowest theta to highest, 
  //for each theta from lowest to highest phi
  int det_num;
  TVector3 gret_pos;//mm
  double core_energy;//kev
  //yta from s800 inverse map
  double yta;//mm
  //beta corrected for DTA
  double beta;
  double ata;
  double bta;
  double gamma; 

  Event(int in_det_num, TVector3 in_gret_pos, double in_energy, double in_yta, 
        double in_beta, double in_ata, double in_bta):
    det_num(in_det_num),
    gret_pos(in_gret_pos),
    core_energy(in_energy),
    yta(in_yta),
    beta(in_beta),
    ata(in_ata),
    bta(in_bta) {
      gamma = 1./(sqrt(1.-pow(this->beta,2.)));
    }

  //Pass x,y,z,ata,bta offsets
  double GetDoppler(double, double, double, double, double) const;
  //pass ata,bta offsets
  TVector3 GetTrack(double, double) const;

};

TVector3 Event::GetTrack(double ata_offset, double bta_offset) const{
  double sata = TMath::Sin(ata + ata_offset);
  double sbta = TMath::Sin(bta + bta_offset);

  TVector3 track(sata, -sbta, sqrt(1-sata*sata-sbta*sbta));
  return track;
}

double Event::GetDoppler(double x_offset, double y_offset, double z_offset, 
                         double ata_offset, double bta_offset) const {

  //Target offsets determine new reference point in lab frame
  TVector3 cur_gret_pos(gret_pos);
  cur_gret_pos.SetX(gret_pos.X() - x_offset);
  cur_gret_pos.SetY(gret_pos.Y() - (y_offset - yta));
  //testing without yta offset
  //cur_gret_pos.SetY(gret_pos.Y() - y_offset);
  cur_gret_pos.SetZ(gret_pos.Z() - z_offset);

  return core_energy*gamma *(1 - beta*TMath::Cos(cur_gret_pos.Angle(GetTrack(ata_offset,bta_offset))));
}


//Readin input file with described format and fill vector with all events
int LoadInputFile(std::string input_fn, std::vector<Event> &all_events){

  std::ifstream input_file(input_fn.c_str());
  if (!input_file.is_open()){
    std::cout << "Failed to open file: " << input_fn << "\n";
    return -1;
  }
  int det_num;
  double gret_pos_x;
  double gret_pos_y;
  double gret_pos_z;
  double energy;
  double yta;
  double beta;
  double ata;
  double bta;
  std::cout << "Loading input file.. \n";
  while (input_file >> det_num >> gret_pos_x >> gret_pos_y >> gret_pos_z 
                    >> energy >> yta >> beta >> ata >> bta){
  
    std::cout <<  det_num << " " << gret_pos_x << " " << gret_pos_y << " " << gret_pos_z 
              << " " << energy << " " << yta << " " << beta << " " << ata << " " << bta << " " << "\n";
    TVector3 gret_pos(gret_pos_x, gret_pos_y, gret_pos_z);
    Event evt(det_num,gret_pos,energy,yta,beta,ata,bta);
    all_events.push_back(evt);
  }
  input_file.close();
  return 0;
}

struct double_iota {
  double_iota(double inc, double init_value = 0.0) : _value(init_value), _inc(inc){}

  operator double() const {return _value; }
  double_iota& operator++() { _value += _inc; return *this; }
  double _value;
  double _inc;
};

void CorrectGretinaShift(std::string input_fn){
  //These are the variables which will determine the range to try shifts, e.g.
  //ATA shift can range from (-ATA_MAG, ATA_MAG) with step size
  //determined by 2*ATA_MAG / N_STEPS 
  //(b/c/ ATA_MAG - (-ATA_MAG) = 2*ATA_MAG) 
  const double ATA_MAG = 0.008; //radians
  const double BTA_MAG = 0.004; //radians
  const double X_MAG = 5.0; //mm
  const double Y_MAG = 1.0;
  const double Z_MAG = 5.0;
  const double N_STEPS = 30; 

  const double ATA_STEP = 2*ATA_MAG / N_STEPS;
  const double BTA_STEP = 2*BTA_MAG / N_STEPS;
  const double X_STEP = 2*X_MAG / N_STEPS;
  const double Y_STEP = 2*Y_MAG / N_STEPS;
  const double Z_STEP = 2*Z_MAG / N_STEPS;


  std::vector<Event> all_events;

  if (LoadInputFile(input_fn, all_events) == -1){
    return;
  }
//  return;

  std::cout << "Loaded: " << all_events.size() << " events\n";
  std::vector<double> ata_vals(std::size_t(((ATA_MAG + ATA_STEP - std::numeric_limits<double>::epsilon()) + ATA_MAG) / ATA_STEP));
  std::vector<double> bta_vals(std::size_t(((BTA_MAG + BTA_STEP - std::numeric_limits<double>::epsilon()) + BTA_MAG) / BTA_STEP));
  std::vector<double> x_vals(std::size_t(((X_MAG + X_STEP - std::numeric_limits<double>::epsilon()) + X_MAG) / X_STEP));
  std::vector<double> y_vals(std::size_t(((Y_MAG + Y_STEP - std::numeric_limits<double>::epsilon()) + Y_MAG) / Y_STEP));
  std::vector<double> z_vals(std::size_t(((Z_MAG + Z_STEP - std::numeric_limits<double>::epsilon()) + Z_MAG) / Z_STEP));
  //Fill vectors with possible values
  std::iota(ata_vals.begin(), ata_vals.end(), double_iota(ATA_STEP,-ATA_MAG));
  std::iota(bta_vals.begin(), bta_vals.end(), double_iota(BTA_STEP,-BTA_MAG));
  std::iota(x_vals.begin(), x_vals.end(), double_iota(X_STEP,-X_MAG));
  std::iota(y_vals.begin(), y_vals.end(), double_iota(Y_STEP,-Y_MAG));
  std::iota(z_vals.begin(), z_vals.end(), double_iota(Z_STEP,-Z_MAG));

  const double DESIRED_ENERGY = 1296;

  struct BestOffsets{
    double x_off;
    double y_off;
    double z_off;
    double ata_off;
    double bta_off;
  };

  BestOffsets best_offsets;
  double best_energy = 0;
  double best_chi_sq = std::numeric_limits<double>::max();
  double crystal_vals[44] = {0};
  double crystal_occurrences[44] {0};
  double xoff = -0.3333;
  double ata_off = 0.002133;
//for (auto &ata_off : ata_vals){
    for (auto &bta_off : bta_vals){
//      for (auto &xoff : x_vals){
        for (auto &yoff : y_vals){
//          for (auto &zoff : z_vals){
            double curr_chi_sq = 0;
            double mean_energy = 0;
            memset(crystal_vals, 0, sizeof(crystal_vals));
            memset(crystal_occurrences, 0, sizeof(crystal_occurrences));
            for (auto &event : all_events){
              crystal_vals[event.det_num] += event.GetDoppler(xoff,yoff,0,ata_off,bta_off);
              //mean_energy += event.GetDoppler(xoff,yoff,0,ata_off,bta_off);
//              std::cout << "Crystal_vals["<<event.det_num<<"] = " << crystal_vals[event.det_num] << "\n";
              crystal_occurrences[event.det_num]++;
            }//loop over events

            //mean_energy = mean_energy / (double)all_events.size();
            for (int i = 0; i < 44; i++){
              if (i == 43) {
                continue;
              }
              crystal_vals[i] /= crystal_occurrences[i];
              mean_energy += crystal_vals[i];
//            std::cout << "Crystal " << i << " has average energy: " << crystal_vals[i] << " and occurrences " << crystal_occurrences[i] << "\n";
            }
            mean_energy /= 43;
            for (int i =0; i < 44; i++) {
              if (i==43) {
                continue;
              }
              curr_chi_sq += pow(crystal_vals[i] - mean_energy, 2);
            }
//          std::cout << "Mean Energy: " << mean_energy << "\n";

            //curr_chi_sq += pow(mean_energy - DESIRED_ENERGY, 2)*10;
            if(curr_chi_sq < best_chi_sq){
              best_offsets.x_off = xoff;
              best_offsets.y_off = yoff;
              best_offsets.z_off = 0;
              best_offsets.ata_off = ata_off;
              best_offsets.bta_off = bta_off;
              best_chi_sq = curr_chi_sq;
              best_energy = mean_energy;
            }
//          std::cout << "==========================================================\n";
        //    return;
//          }//loop over z
        }//loop over y
//      }//loop over x
    }//loop over bta
//  }//loop over ata

  std::cout << "Best Energy = " << best_energy << "\n";
  std::cout << "Best Values (chi-2 = " << best_chi_sq << ") are:\n";
  std::cout << "ATA Offset: " << best_offsets.ata_off << "\n";
  std::cout << "BTA Offset: " << best_offsets.bta_off << "\n";
  std::cout << "x Offset: "   << best_offsets.x_off << "\n";
  std::cout << "y Offset: "   << best_offsets.y_off << "\n";
  std::cout << "z Offset: "   << best_offsets.z_off << "\n";
}//end of function

int main(int argc, char **argv){

  std::cout << "Running... " << std::endl;
  if (argc < 2){
    std::cout << "USAGE: CorrectGretinaOffset INPUT_FN\n";
  }
  CorrectGretinaShift(argv[1]);
}
