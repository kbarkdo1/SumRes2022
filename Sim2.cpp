// Welcome to the Machine

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <unistd.h>
#include <cmath>

using namespace std;

/*
Overall problems:
Firing is calculated without checking inhibition first. A neuron that is inhibited
won't see that inhibition if it fires from just the input and ICs.

Redundancy: input spiking is calculated twice-as a spiking child
and again for the total updating)(and discarded the first time)

Switching conditionals for file writing:
Metric is above 1/3 for 10 consecutive cases
Leave when metric is below 1/3 for 5 cases

*/


void init_global(float* connect, float* image, float* voltages, int* fired, int neur_num, int excitatory);
void chase_spikes(vector<int> to_check, float* voltages, float* delta, int* fired, float* connect, int neur_num, float s, float* threshold);

int main() {
  // initialize all variables
  // make neuron connections
  float forcings[20];
  for(int i=0; i<20; i++) {
    forcings[i] = 0.2*(i+1);
  }

  float lambs[50];
  float phis[50];
  for(int i=0; i< 50; i++) {
    lambs[i] = float(i)/100;
    phis[i] = float(i)/100;
  }

  ofstream f_vals;
  ofstream l_vals;
  ofstream p_vals;
  f_vals.open("forcings.txt");
  l_vals.open("lambdas.txt");
  p_vals.open("phis.txt");
  for(int i=0; i < 20; i++) {
    f_vals << forcings[i] << endl;
  }
  for(int i=0; i< 50; i++) {
    l_vals << lambs[i] << endl;
    p_vals << phis[i] << endl;
  }
  f_vals.close();
  l_vals.close();
  p_vals.close();

  int tn =30;//  and 928    //TRACKING NEURON
  int tn2=2;
  int round_counter=0;
  float fir_ex=0;
  float fir_in=0;
  int RC2=0;
  float temp_ex=0;
  float temp_in=0;
  
  int excitatory = 1000;
  int inhibitory = 1000;
  int neur_num = excitatory + inhibitory;
  float cycles = 40.00;
  float step = 0.01;
  float* connect = new float [neur_num*neur_num];
  float* image = new float [neur_num];
  float* voltages = new float [neur_num];
  int* fired = new int [neur_num];
  init_global(connect, image, voltages, fired, neur_num, excitatory);
  int fileCount =0;

  ofstream for_fir;
  for_fir.open("for_fir.txt");
  float* k1 = new float [neur_num];
  float* k2 = new float [neur_num];
  float* delta = new float [neur_num];
  
  vector<int> num_firing;
  vector<float> avg_volts;
  float threshold[neur_num];
  vector<float> neur_thresh; // rolling sum for average neuron threshold per neuron
  for(int i=0; i< neur_num; i++) {
    num_firing.push_back(0);
    threshold[i] = 1;
    avg_volts.push_back(0);
    neur_thresh.push_back(0);
  }
  

  for(int z=4; z < 5; z++) {
    for(int phi_index=30; phi_index < 31; phi_index+=10) {
      for(int lamb_index=5; lamb_index < 6; lamb_index+=5) {
        fileCount++;
        float f = forcings[z];
        // float phi = phis[phi_index];
        // float lambda  = lambs[lamb_index];
        float phi = 0.005;
        float lambda = 0.00625;
        float s = neur_num;

        string avg_voltages_str = "S2_avg_volages_"+to_string(fileCount)+".txt";
        string avg_neur_thresh_str= "S2_avg_neur_thresholds_"+to_string(fileCount)+".txt";
        string f_rate_str= "S2_f_rate_"+to_string(fileCount)+".txt";
        string avg_fir_str= "S2_avg_fir_"+to_string(fileCount)+".txt";
        ofstream avg_fir;
        ofstream act_ex;
        ofstream act_in;
        avg_fir.open(avg_fir_str);
        act_ex.open("act_ex"+to_string(fileCount)+".txt");
        act_in.open("act_in"+to_string(fileCount)+".txt");
        ofstream tn_ex;
        ofstream tn_in;
        tn_ex.open("tn_ex"+to_string(fileCount)+".txt");
        tn_in.open("tn_in"+to_string(fileCount)+".txt");

        ofstream avg_thresh;
        ofstream tn_thresh;
        ofstream tn_voltage;
        avg_thresh.open("avg_thresh"+to_string(fileCount)+".txt");
        tn_thresh.open("tn_thresh"+to_string(fileCount)+".txt");
        tn_voltage.open("tn_voltage"+to_string(fileCount)+".txt");

        float thresh_acc = 0;

        for(int i=0; i< neur_num; i++) {
          num_firing[i]=0;
          threshold[i] = 1;
          // ex_sum[i]=0;
          // in_sum[i]=0;
          avg_volts[i]=0;
          neur_thresh[i]=0;
        }

        for(float t=0; t<cycles; t = t+step) {
          // printf("Time = %2.3f \n", t);
          // printf("ex_sum, in_sum, %d: %2.4f %2.4f  ", tn, ex_sum[tn], in_sum[tn]);
          // printf("thresh %d: %2.4f \n", tn, threshold[tn]);
          tn_thresh << threshold[tn] << endl;
          tn_voltage << voltages[tn] << endl;


          // for every neuron

          vector<int> mybabies; // list of neurons that need firing

          for(int n=0; n < neur_num; n++) {
            thresh_acc += threshold[n]; // ADDED FOR Threshold Plot
            avg_volts[n] += voltages[n]; // rolling sum for eventual average volts
            neur_thresh[n] += threshold[n];
            
            // calculate voltage steps
            float volt_diff = -voltages[n]; // calculate resting voltage difference
            float img_in = f * image[n];// calculate image input
            //   k1[n] = volt_diff + img_in;
            k1[n] = step * (volt_diff+img_in);
            
            float volt_diff_2 = -(voltages[n]+k1[n]); // calculate k2 voltage difference
            float img_in_2 = f * image[n];  // calculate image input
            k2[n] = step * (volt_diff_2 + img_in_2);

            delta[n] = 0.5 * (k1[n] + k2[n]);

            if (voltages[n] + delta[n] >= threshold[n]) { // checks for spike
              fired[n] = 1; // spiked list
              int from = n; // syntax formality
              for(int target=0; target < neur_num; target++) {
                if (connect[target * neur_num + from]) {
                  mybabies.push_back(target);
                }
              }
            }
          }
          avg_thresh << thresh_acc / neur_num << endl;
          thresh_acc = 0;

          chase_spikes(mybabies, voltages, delta, fired, connect, neur_num, s, threshold);
          
          int fired_now = 0;
          for(int i=0; i< neur_num; i++) { // voltage updating with delta
            if (fired[i]==0) { // updates if not fired
              float neighbor_in = 0;
              int target = i;

              for(int from = 0; from < (neur_num); from++) {
                neighbor_in = neighbor_in + connect[target * neur_num + from]
                  * fired[from]; // sums neighbors if neighbors fire.
              }
              float s_N = s / float(neur_num); // S/N from DE
              neighbor_in = neighbor_in * s_N; // neighbor adjustment
              
              if (i==tn && neighbor_in !=0 ) {
                // printf(" %d receieved %2.4f fron neighs \n", tn, neighbor_in);
              }
              /*
              if (i==tn2 && neighbor_in !=0 ) {
                printf(" %d receieved %2.4f fron neighs \n", tn2, neighbor_in);
              }
              */

              if (voltages[i]+delta[i]+neighbor_in >= threshold[i]) { // failed overrun double check
                printf("unfired overrun, %d, %2.6f, %2.6f\n", i, voltages[i]+delta[i], voltages[i]+delta[i]+neighbor_in);
              } else {
                voltages[i] = voltages[i]+delta[i]+neighbor_in; // forward step
                float r1 = (-1 * lambda * (threshold[i] - 1)) * step;
                float r2 = (-1 * lambda * ((threshold[i] + r1) - 1)) * step;
                float rdelt = (r1 + r2) * 0.5;
                threshold[i] += rdelt;
              }
            } else if (fired[i]==1) { // resets spikes
              if (i == tn) {
                printf("Time at tn firing: %2.4f \n", t);
              }
              num_firing[i] = num_firing[i] + 1;
              voltages[i]=0;
              threshold[i] += phi;
            }
          }
          
          for(int i=0; i<excitatory; i++) {
            fir_ex+= fired[i];
          }
          for(int i=excitatory; i < neur_num; i++) {
            fir_in += fired[i];
          }
          
          round_counter++;
          if (round_counter==40) {
            act_ex << float(fir_ex)/float(excitatory) << endl;
            act_in << float(fir_in)/float(inhibitory) << endl;
            round_counter=0;
            fir_ex=0;
            fir_in=0;
          }
          
          fired_now=0;
          
          for(int i=0; i<neur_num; i++) {
            float ex_input_add=0;
            float in_input_add=0;
            for(int j=0; j<excitatory; j++) {
              ex_input_add += connect[i*neur_num + j] * fired[j];
              // printf(" ex: i: %d, con: %2.4f, fired: %d, ex_input_add: %2.4f, ex_sum: %2.4f", i, connect[i*neur_num+j], fired[j], ex_input_add, ex_sum[i]);
            }
            // ex_sum[i] = ex_sum[i] + ex_input_add;
            
            for(int j=excitatory; j<neur_num; j++) {
              in_input_add += connect[i*neur_num + j] * fired[j];
            }
            // in_sum[i] = in_sum[i] + in_input_add;

            if (i == tn) {
              RC2++;
              temp_ex += ex_input_add;
              temp_in += in_input_add;
              if (RC2 == 50) {
                tn_ex << temp_ex << endl;
                tn_in << temp_in << endl;
                RC2 = 0;
                temp_ex = 0;
                temp_in = 0;
              }
            }
          }
          

          for(int i=0; i< neur_num; i++) {
            fired[i]=0;
          }

        }
        float n_avg=0;
        float E_avg = 0;
        float I_avg = 0;
        ofstream f_rate;
        f_rate.open(f_rate_str);
        // printf("\n %2.4f %2.4f \n", cycles, (cycles)); 
        for(int i=0; i<neur_num; i+=1) {
          n_avg=float(num_firing[i]);
          n_avg = n_avg / cycles;

          f_rate << n_avg << endl;
          if (i<excitatory) {
            E_avg += n_avg;
          } else {
            I_avg += n_avg;
          }
          
          n_avg=0;
          // printf("total_avg: %d %2.4f  ", i, total_avg);
        }
        f_rate.close();
        // printf("%2.6f %d \n", total_avg, neur_num);
        E_avg = E_avg / (excitatory);
        I_avg = I_avg / (inhibitory);
        for_fir << E_avg << " " << I_avg << endl;

        ofstream avg_voltages;
        avg_voltages.open (avg_voltages_str);
        ofstream neur_thresholds;
        neur_thresholds.open (avg_neur_thresh_str);

        for(int i=0; i<neur_num; i++) {
          
          avg_voltages << avg_volts[i]/(cycles/step) << endl;
          neur_thresholds << neur_thresh[i]/(cycles/step) << endl;
        }
        neur_thresholds.close();
        avg_voltages.close();
        act_ex.close();
        act_in.close();

        /*
        ofstream ex_sum_out;
        ex_sum_out.open ("ex_sum.txt");
        ofstream in_sum_out;
        in_sum_out.open ("in_sum.txt");
        for(int i=0; i<neur_num; i++) {
          ex_sum_out << ex_sum[i]/cycles + image[i] << endl;
          in_sum_out << in_sum[i]/cycles << endl;
        }
        ex_sum_out.close();
        in_sum_out.close();
        */
        
      }
    }

  // Outputs:

  }
  for_fir.close();
  delete[] k1;
  delete[] k2;
  delete[] delta;

    // deallocation, baby!
  /*
  delete[] connect;
  delete[] image;
  delete[] voltages;
  delete[] fired;
  delete[] vec_arr;
  */
  return 0;
};

void init_global(float* connect, float* image, float* voltages, int* fired, int neur_num, int excitatory) {

  // float check_ratio=0;
  // 0<=n < 1000 = Excitatory
  // // 1000 <= n < 2000 = Inhibitory
  
  float j_ee = 1;
  float j_ei = -2;
  float j_ii = -1.8;
  float j_ie = 1;
  float K = 40;
  int inhibitory = neur_num - excitatory;
  for(int from = 0; from < excitatory; from++) {
    for(int target = 0; target < excitatory; target++) {
      int seed = rand() % excitatory/K;
      if (seed == 1) {
        connect[target*neur_num + from] = float(j_ee)/sqrt(K);
        // printf(" %d", connect[target*neur_num + from]);
        // check_ratio++;
      } else {
        connect[target*neur_num + from] = 0;
      }
    }
  }
  
  for(int from = excitatory; from < neur_num; from++) {
    for(int target = 0; target < excitatory; target++) {
      int seed = rand() % inhibitory/K;
      if (seed == 1) {
        connect[target*neur_num + from] = float(j_ei)/sqrt(K);
        // printf(" %d", connect[target*neur_num + from]);
        // check_ratio++;
      } else {
        connect[target*neur_num + from] = 0;
      }
    }
  }

  for(int from = excitatory; from < neur_num; from++) {
    for(int target = excitatory; target < neur_num; target++) {
      int seed = rand() % inhibitory/K;
      if (seed == 1) {
        connect[target*neur_num + from] = float(j_ii)/sqrt(K);
        // printf(" %d", connect[target*neur_num + from]);
        // check_ratio++;
      } else {
        connect[target*neur_num + from] = 0;
      }
    }
  }

  for(int from = 0; from < excitatory; from++) {
    for(int target = excitatory; target < neur_num; target++) {
      int seed = rand() % excitatory/K;
      if (seed == 1) {
        connect[target*neur_num + from] = float(j_ie)/sqrt(K);
        // printf(" %d", connect[target*neur_num + from]);
        // check_ratio++;
      } else {
        connect[target*neur_num + from] = 0;
      }
    }
  }

  // printf("check ratio: %3.6f", float(check_ratio)/(neur_num*neur_num));

  for(int i=0; i < neur_num; i++) {   // set horiztonal to zero, no reccurence
    connect[i*neur_num + i] = 0;
  }
  ofstream A;
  A.open ("connectivity.txt" );
  for(int i=0; i< neur_num; i++) {
    for(int j=0; j< neur_num; j++) {
      A << connect[i*neur_num + j] << ' ';
    }
    A << endl;
  }
  A.close();
  /*
  // file read in.
  float f_e = 1.25;
  float f_i = 1;
  // float m_0 = float(1)/float(750);
  float m_0 = float(1/152.000);
  ifstream library_reader;
  library_reader.clear();
  library_reader.open("im_stripes.txt");
  double temp=0;
  vector<float> IM;
  while(!library_reader.eof()) {     //create/fill image vector, IM (c1)
    library_reader >> temp;
    IM.push_back(float(temp));     //vector<double> IM;  //pixel matrix put into vector representation
    // printf("%2.6f\n", float(temp));
  }
  library_reader.close();
  int img_vec_len = IM.size() - 1;
  printf("   IM.size = %zu   ", IM.size());
  float img_vec[img_vec_len];
  for(int i=0; i<img_vec_len; i++) {
    img_vec[i] = IM[i];
    // printf(" img_vec %d: %4.4f \n", i, img_vec[i][0]);
  }
  float* B = new float [neur_num*img_vec_len];
  printf("IM length: %d\n", int(IM.size()));
  printf("B length: %d\n", neur_num*img_vec_len);
  
  ofstream Be;
  Be.open ("B.txt");
  
  for(int i=0; i< neur_num; i++) {
    for (int j=0; j<img_vec_len/2; j++) {
      int path = rand() % 2000;
      if (path == 1) {
        // all this stuff below is for dynamic connection strengths
        
        B[i*img_vec_len + j] = 1;
        Be << "1 ";
        // if (i==149 && j==5) {printf("\n at 149,5 \n");}
      } else {
        B[i*img_vec_len + j] = 0;
        Be << "0 ";
      }

    }

    for (int j=img_vec_len/2; j<img_vec_len; j++) {
      int path = rand() % 2000;
      if (path == 1) {
        // all this stuff below is for dynamic connection strengths
        
        B[i*img_vec_len + j] = 1;
        Be << "1 ";
        // if (i==149 && j==5) {printf("\n at 149,5 \n");}
      } else {
        B[i*img_vec_len + j] = 0;
        Be << "0 ";
      }

    }
    
    Be << endl;
  }
  Be.close();
  
  for(int i = 0; i < excitatory; i++) {
    image[i]=0; 
    // printf(" image %d, %2.4f  \n", i, image[i]);
    for(int k = 0; k < 10000; k++) {
      image[i] += B[i*img_vec_len + k] * IM[k];
    }
    image[i] = image[i] * f_e * m_0;
    // printf(" image %d, %2.4f  \n", i, image[i]);
  }

  for(int i = excitatory; i < neur_num; i++) {
    image[i]=0; 
    // printf(" image %d, %2.4f  \n", i, image[i]);
    for(int k = 0; k < 10000; k++) {
      image[i] += B[i*img_vec_len + k] * IM[k];
    }
    image[i] = image[i] * f_i * m_0;
    // printf(" image %d, %2.4f  \n", i, image[i]);
  }
  */
  for(int i = 0; i < neur_num/2; i++) {  // make image input
    image[i] = 1*1*sqrt(K);
  }
  for(int i=neur_num/2; i<neur_num; i++) {
    image[i]=0.8*1*sqrt(K);
  }
  
  ofstream fBp;
  fBp.open ("fBp.txt" );
  for(int i=0; i<neur_num; i++) {
    fBp << image[i] << endl;
  }
          
  fBp.close();
  //read in image C++

  for(int i = 0; i < neur_num; i++) {   // make neuron voltages
    int seed = rand() % 1000;
    float volt = seed / 1000.0000;
    voltages[i] = volt;
  }

  for(int i =0; i< neur_num; i++) { // set fired to zero
    fired[i] = 0;
  }


}

void chase_spikes(vector<int> to_check, float* voltages, float* delta, int* fired, float* connect, int neur_num, float s, float* threshold) {
  // printf("%s\n", "in chase spikes \n");
  vector<int> mybabies;
  int num_fi = 0;
  if (to_check.size()< 1) { // base case
    return;
  }

  for(int i=0; i<to_check.size(); i++) {
    if (fired[to_check[i]]) { // if its already fired
      num_fi++;
      continue;
    }

    int here = to_check[i];
    float neighbor_in = 0;
    int target = here;

    for(int from = 0; from < (neur_num); from++) { // sums up fired inputs
      neighbor_in = neighbor_in + connect[target * neur_num + from]
        * fired[from];
    }
    /*
    if (to_check[i] == 915) {
      for(int from = 0; from < (neur_num); from++) { // sums up fired inputs
        if (connect[target * neur_num + from] * fired[from]) {
          printf("from: %d\n", from);
        }
      }
    }
    */
    float s_N = s / float(neur_num); // constant calculation
    neighbor_in = neighbor_in * s_N; // constant calculation

    if (voltages[here] + delta[here] + neighbor_in >= threshold[here]) { // check for new firing
      if (fired[here]!=1) {
        fired[here] = 1; // sets as fired
        int from = here;
        // push children for analysis
        for(int target=0; target < neur_num; target++) { // adds anyone who receives output
          if (connect[target * neur_num + from]) {
            mybabies.push_back(target);
          }
        }

      } else {
        printf("%s\n", "refiring in chase_spikes");
      }
    }
  }

  chase_spikes(mybabies, voltages, delta, fired, connect, neur_num, s, threshold); // recursion
  return;
}