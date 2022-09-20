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

*/


void init_global(float* connect, float* image, float* voltages, int* fired, int neur_num, int ex1, int ex2, int in1, int in2);
void chase_spikes(vector<int> to_check, float* voltages, float* delta, int* fired, float* connect, int neur_num, float s, float* threshold, float* neighbor_in);

int main() {
  // initialize all variables
  // make neuron connections
  srand(40);
  float forcings[20];
  forcings[0] = 1.000;
  forcings[1] = .9;
  forcings[2] = .8;
  forcings[3] = .7;
  forcings[4] = .6;
  forcings[5] = .5;
  forcings[6] = .4;
  forcings[7] = .3;
  forcings[8] = .2;
  forcings[9] = .1;
  forcings[10] = 1;

  float lambs[50];
  float phis[50];
  for(int i=0; i< 50; i++) {
    lambs[i] = float(i)/100;
    phis[i] = float(i)/100;
  }

  /*
  ofstream f_vals;
  ofstream l_vals;
  ofstream p_vals;
  f_vals.open("forcings.txt");
  l_vals.open("lambdas.txt");
  p_vals.open("phis.txt");
  f_vals << f << endl;
  l_vals << lambda << endl;
  p_vals << phi << endl;
  f_vals.close();
  l_vals.close();
  p_vals.close();
  */
  int file_count = 0;
  float m_2 = 1;
  for (float h =1; h < 1.01; h += 0.025) {
  // for (int h = 40; h < 50; h += 10) {
    float m_1 = h;
    // srand(h);
    
    file_count++;
    /*
    ofstream events;
    events.open("events.txt");
    ofstream times;
    times.open("times.txt");
    */

    ofstream pos_times;
    pos_times.open("pos_times_" + to_string(file_count) +".txt");
    ofstream neg_times;
    neg_times.open("neg_times_" + to_string(file_count) +".txt");
    // printf("1st force = %2.2f \n", h);
    
    int tn =13;//  and 928    //TRACKING NEURON
    int tn2=2;
    int round_counter=0;
    float fir_ex_1=0;
    float fir_in_1 = 0;
    float fir_ex_2=0;
    float fir_in_2 = 0;
    
    int ex1 = 1000;
    int in1 = 1000;
    int ex2 = 1000;
    int in2 = 1000;
    int neur_num = ex1 + ex2 + in1 + in2;
    float cycles = 500.00;
    float step = 0.01;

    float* connect = new float [neur_num*neur_num];
    float* image = new float [neur_num];
    float* voltages = new float [neur_num];
    int* fired = new int [neur_num];
    init_global(connect, image, voltages, fired, neur_num, ex1, ex2, in1, in2);
    

    float neighbor_in[neur_num];
    int num_firing[neur_num];
    int nf_temp[neur_num];
    // int nfn_temp[neur_num];
    float avg_volts_pos[neur_num];
    float avp_temp[neur_num];
    float avg_volts_neg[neur_num];
    float avn_temp[neur_num];
    float ex_sum[neur_num];
    float in_sum[neur_num]; 
    float threshold[neur_num];
    float neur_thresh_pos[neur_num]; // rolling sum for average neuron threshold per neuron
    float ntp_temp[neur_num];
    float neur_thresh_neg[neur_num];
    float ntn_temp[neur_num];
    int pos_collection=0; //number of time steps that we have collected positive data for
    int neg_collection=0;
    float metric = 0;
    float ae1 = 0;
    float ae2 = 0;
    float met = 0;
    int met_counter_pos = 0;
    int met_counter_neg = 0;
    int last_round = 0;
    int write_to = 0;
    int nf_count = 0;
    int pf_count = 0;
    float lt = 0;
    float atp1 = 0;
    float atn1 = 0;
    float atp2 = 0;
    float atn2 = 0;

    for(int z=0; z<1; z++) {
      for(int phi_index=0; phi_index < 1; phi_index++) {
        for(int lamb_index=0; lamb_index < 1; lamb_index++) {
          // file_count++;
          float f = forcings[z];
          // float phi = phis[phi_index];
          float phi = 0.3;
          // float phi = 0;
          // float lambda  = lambs[lamb_index];
          float lambda = 0.005;
          // float lambda = 0;
          float s = neur_num;

          printf("f = %2.4f, phi = %2.4f, lambda = %2.4f \n", f, phi, lambda);

          ofstream avg_fir;
          ofstream act_ex_1;
          ofstream act_ex_2;
          ofstream act_in_1;
          ofstream act_in_2;
          ofstream metric;
          string avg_fir_str = "avg_fir_"+to_string(file_count)+".txt";
          string metric_str = "metric_"+to_string(file_count)+".txt";
          string act_ex_1_str = "act_ex_1_"+to_string(file_count)+".txt";
          string act_ex_2_str = "act_ex_2_"+to_string(file_count)+".txt";
          avg_fir.open(avg_fir_str);
          act_ex_1.open(act_ex_1_str);
          act_ex_2.open(act_ex_2_str);
          act_in_1.open("act_in_1_"+to_string(file_count)+".txt");
          act_in_2.open("act_in_2_"+to_string(file_count)+".txt");
          metric.open(metric_str);

          ofstream avg_voltages;
          ofstream neur_thresholds;
          ofstream f_rate;
          ofstream pos1_at;
          ofstream neg1_at;
          ofstream pos2_at;
          ofstream neg2_at;
          pos1_at.open("atp1_"+to_string(file_count)+".txt");
          neg1_at.open("ntp1_"+to_string(file_count)+".txt");
          pos2_at.open("atp2_"+to_string(file_count)+".txt");
          neg2_at.open("ntp2_"+to_string(file_count)+".txt");
          ofstream pos1_fr;
          ofstream neg1_fr;
          ofstream pos2_fr;
          ofstream neg2_fr;
          pos1_fr.open("frp1_"+to_string(file_count)+".txt");
          neg1_fr.open("frn1_"+to_string(file_count)+".txt");
          pos2_fr.open("frp2_"+to_string(file_count)+".txt");
          neg2_fr.open("frn2_"+to_string(file_count)+".txt");
          

          

          /*
          string avg_voltages_str = "avg_volages_"+to_string(0)+".txt";
          string avg_neur_thresh_str = "avg_neur_thresholds_"+to_string(0)+".txt";
          string f_rate_str = "f_rate_"+to_string(0)+".txt";
          
          avg_voltages.open(avg_voltages_str);
          neur_thresholds.open(avg_neur_thresh_str);
          f_rate.open(f_rate_str);
          */
          ofstream t_tr;
          t_tr.open("thresh_trace.txt");
          for(int i=0; i< neur_num; i++) {
            neighbor_in[i] = 0;
            num_firing[i]=0;
            nf_temp[i]=0;
            // nfn_temp[i]=0;
            avg_volts_pos[i]=0;
            avp_temp[i]=0;
            avg_volts_neg[i]=0;
            avn_temp[i]=0;
            ex_sum[i]=0;
            in_sum[i]=0; 
            threshold[i]=1;
            neur_thresh_pos[i]=0; // rolling sum for average neuron threshold per neuron
            ntp_temp[i]=0;
            neur_thresh_neg[i]=0;
            ntn_temp[i]=0;
          }

          for(int i=ex1; i< ex1+in1; i++) {
            threshold[i] = 0.8;
          }
          for(int i=ex1+in1+ex2; i < neur_num; i++) {
            threshold[i] = 0.8;
          }
          float* k1 = new float [neur_num];
          float* k2 = new float [neur_num];
          float* delta = new float [neur_num];
          
          for (float t=0; t<cycles; t = t+step) {
            
            //printf("Time = %2.3f \n", t);
            // printf("ex_sum, in_sum, %d: %2.4f %2.4f  ", tn, ex_sum[tn], in_sum[tn]);
            // printf("thresh %d: %2.4f \n", tn, threshold[tn]);
            // tn_thresh << threshold[tn] << endl;
            // for every neuron

            vector<int> mybabies; // list of neurons that need firing
            t_tr << threshold[tn] << endl;
            for(int n=0; n < neur_num; n++) {
              /*
              if (n == 200) {
                printf("%2.6f \n", voltages[n]);
              }
              */
              // avg_volts[n] += voltages[n]; // rolling sum for eventual average volts
              // neur_thresh[n] += threshold[n];
              // calculate voltage steps
              if (n < ex1) {
                atp1 += threshold[n];
              } else if (n < ex1+in1) {
                atn1 += threshold[n];
              } else if (n < ex1+in1+ex2) {
                atp2 += threshold[n];
              } else {
                atn2 += threshold[n];
              }

              float volt_diff = -voltages[n]; // calculate resting voltage difference
              float img_in;
              if (n < ex1+in1) {
                img_in = m_1 * image[n];// calculate image input
              } else {
                img_in = m_2 * image[n];
              }
              
              //   k1[n] = volt_diff + img_in;
              k1[n] = step * (volt_diff+img_in);
              
              float volt_diff_2 = -(voltages[n]+k1[n]); // calculate k2 voltage difference
              float img_in_2;
              if (n < ex1+in1) {
                img_in_2 = m_1 * image[n];// calculate image input
              } else {
                img_in_2 = m_2 * image[n];
              }
              k2[n] = step * (volt_diff_2 + img_in_2);

              delta[n] = 0.5 * (k1[n] + k2[n]);

              if (voltages[n] + delta[n] >= threshold[n]) { // checks for spike
                fired[n] = 1; // spiked list
                int from = n; // syntax formality
                for(int target=0; target < neur_num; target++) {
                  if (connect[target * neur_num + from]) {
                    mybabies.push_back(target);
                    neighbor_in[target] += connect[target * neur_num + from] * fired[from];
                  }
                }
              }
            }

            pos1_at << atp1 / (ex1) << endl;
            neg1_at << atn1 / (in1) << endl;
            pos2_at << atp2 / (ex2) << endl;
            neg2_at << atn2 / (in2) << endl;
            atp1 = 0;
            atn1 = 0;
            atp2 = 0;
            atn2 = 0;
            // printf("at chase spikes \n");
            chase_spikes(mybabies, voltages, delta, fired, connect, neur_num, s, threshold, neighbor_in);
            // printf("after chase spikes \n");
            int fired_now = 0;
            for(int i=0; i< neur_num; i++) { // voltage updating with delta
              if (i < ex1) {
                avp_temp[i] += voltages[i];
                ntp_temp[i] += threshold[i];
              } else if (i < ex1+in1) {
                avp_temp[i] += voltages[i];
                ntp_temp[i] += threshold[i];
              } else if (i < ex1+in1+ex2) {
                avn_temp[i] += voltages[i];
                ntn_temp[i] += threshold[i];
              } else {
                avn_temp[i] += voltages[i];
                ntn_temp[i] += threshold[i];
              }
              if (fired[i]==0) { // updates if not fired
                /*
                float neighbor_in = 0;
                int target = i;

                for(int from = 0; from < (neur_num); from++) {
                  neighbor_in = neighbor_in + connect[target * neur_num + from]
                    * fired[from]; // sums neighbors if neighbors fire.
                }
                float s_N = s / float(neur_num); // S/N from DE
                neighbor_in = neighbor_in * s_N; // neighbor adjustment
                
                if (i==tn && neighbor_in !=0 ) {
                  printf(" %d receieved %2.4f fron neighs \n", tn, neighbor_in);
                }
                */
                if (voltages[i]+delta[i]+neighbor_in[i] >= threshold[i]) { // failed overrun double check
                  printf("unfired overrun, %d, %2.6f, %2.6f\n", i, voltages[i]+delta[i], voltages[i]+delta[i]+neighbor_in[i]);
                } else {
                  voltages[i] = voltages[i]+delta[i]+neighbor_in[i]; // forward step
                  if (i < ex1 || (i >= ex1+in1 && i < ex1+in1+ex2)) {
                    float r1 = (-1 * lambda * (threshold[i] - 1)) * step;
                    float r2 = (-1 * lambda * ((threshold[i] + r1) - 1)) * step;
                    float rdelt = (r1 + r2) * 0.5;
                    threshold[i] += rdelt;
                  } else {
                    float r1 = (-1 * lambda * (threshold[i] - 0.8)) * step;
                    float r2 = (-1 * lambda * ((threshold[i] + r1) - 0.8)) * step;
                    float rdelt = (r1 + r2) * 0.5;
                    threshold[i] += rdelt;
                  }
                }
              } else if (fired[i]==1) { // resets spikes
                // num_firing_pos[i] = num_firing[i] + 1;
                // events << i << endl;
                // times << t << endl;
                if (i == tn) {
                  printf("%4.2f\n", t);
                }
                voltages[i]=0;
                threshold[i] += phi;
                nf_temp[i] += nf_temp[i] + 1;
                if (i < ex1) {
                  fir_ex_1 += fired[i];
                } else if (i < ex1+in1) {
                  fir_in_1 += fired[i];
                } else if (i < ex1+in1+ex2) {
                  fir_ex_2 += fired[i];
                } else {
                  fir_in_2 += fired[i];
                }
              }
            }
            
            round_counter++;
            if (round_counter==50) {
              // printf("%2.4f \n", lt);
              act_ex_1 << float(fir_ex_1)/float(ex1) << endl;
              ae1 = float(fir_ex_1)/float(ex1);
              act_ex_2 << float(fir_ex_2)/float(ex2) << endl;
              ae2 = float(fir_ex_2)/float(ex2);
              act_in_1 << float(fir_in_1)/float(in1) << endl;
              act_in_2 << float(fir_in_2)/float(in2) << endl;
              metric << float((float(fir_ex_1)/float(ex1))-(float(fir_ex_2)/float(ex2))) / float((float(fir_ex_1)/float(ex1))+(float(fir_ex_2)/float(ex2))) << endl;
              met = float((float(fir_ex_1)/float(ex1))-(float(fir_ex_2)/float(ex2))) / float((float(fir_ex_1)/float(ex1))+(float(fir_ex_2)/float(ex2)));
              round_counter=0;
              fir_ex_1=0;
              fir_ex_2=0;
              fir_in_1=0;
              fir_in_2=0;

              if (write_to == 1) {
                for (int i = 0; i < neur_num; i++) {
                  num_firing[i] += nf_temp[i];
                  avg_volts_pos[i] += avp_temp[i];
                  neur_thresh_pos[i] += ntp_temp[i];
                }    
              }
              if (write_to == -1) {
                for (int i = 0; i < neur_num; i++) {
                  num_firing[i] += nf_temp[i];
                  avg_volts_neg[i] += avn_temp[i];
                  neur_thresh_neg[i] += ntn_temp[i];
                } 
              }
              for (int i = 0; i < neur_num; i++) {
                nf_temp[i] = 0;
                avp_temp[i] = 0;
                ntp_temp[i] = 0;
                avn_temp[i] = 0;
                ntn_temp[i] = 0;
              }
              if (met > 0.4) {
                if (last_round == 1) {
                  met_counter_pos += 1;
                } else {
                  last_round = 1;
                  if (write_to != 1) {
                    met_counter_pos = 1;
                  }
                }
              }
              if (met < -0.4) {
                if (last_round == -1) {
                  met_counter_neg += 1;
                } else {
                  last_round = -1;
                  if (write_to != -1) {
                    met_counter_neg = 1;
                  }
                }
              }
              if (met_counter_pos == 5 && write_to != 1) { // only open and straight write to them
                // Close negavtive and OPEN POSITIVE FILES
                nf_count -= 1;
                string avg_voltages_str = "avg_volages_"+to_string(nf_count)+".txt";
                string avg_neur_thresh_str = "avg_neur_thresholds_"+to_string(nf_count)+".txt";
                string f_rate_str = "f_rate_"+to_string(nf_count)+".txt";
          
                avg_voltages.open(avg_voltages_str);
                neur_thresholds.open(avg_neur_thresh_str);
                f_rate.open(f_rate_str);
                float n_avg = 0;
                float total_avg = 0;
                for(int i=0; i<neur_num; i++) {
                  avg_voltages << avg_volts_neg[i]/((t - lt)/step) << endl;
                  neur_thresholds << neur_thresh_neg[i]/((t - lt)/step) << endl;
                  n_avg=float(num_firing[i]);

                  // printf("n_avg_sum = %3.4f  ", n_avg);
                  n_avg = n_avg / (t - lt);
                  // printf("n_avg = %3.4f , t-lt = %3.4f ", n_avg, t - lt);
                  f_rate << n_avg << endl;
                  total_avg+=n_avg;
                  n_avg=0;

                  avg_volts_neg[i] = 0;
                  neur_thresh_neg[i] = 0;
                  num_firing[i] = 0;
                }
                neur_thresholds.close();
                avg_voltages.close();
                f_rate.close();

                neg_times << (t - lt) << endl;

                // printf("Time: %2.4f, wrote_to == %d, file = %d \n", t, write_to, nf_count);
                
                write_to = 1;
                met_counter_neg = 0;
                lt = t;
              }
              if (met_counter_neg == 5 && write_to != -1) {
                // close positive and open NEGATIVE FILES
                pf_count += 1; 
                string avg_voltages_str = "avg_volages_"+to_string(pf_count)+".txt";
                string avg_neur_thresh_str = "avg_neur_thresholds_"+to_string(pf_count)+".txt";
                string f_rate_str = "f_rate_"+to_string(pf_count)+".txt";
          
                avg_voltages.open(avg_voltages_str);
                neur_thresholds.open(avg_neur_thresh_str);
                f_rate.open(f_rate_str);
                float n_avg = 0;
                float total_avg = 0;
                for(int i=0; i<neur_num; i++) {
                  avg_voltages << avg_volts_pos[i]/((t - lt)/step) << endl;
                  neur_thresholds << neur_thresh_pos[i]/((t - lt)/step) << endl;
                  n_avg=float(num_firing[i]);
                  n_avg = n_avg / (t - lt);

                  f_rate << n_avg << endl;
                  total_avg+=n_avg;
                  n_avg=0;
                  avg_volts_pos[i] = 0;
                  neur_thresh_pos[i] = 0;
                  num_firing[i] = 0;
                }


                neur_thresholds.close();
                avg_voltages.close();
                f_rate.close();

                pos_times << (t - lt) << endl;

                // printf("Time: %2.4f, wrote_to == %d, file = %d \n", t, write_to, pf_count);

                write_to = -1;
                met_counter_pos = 0;
                lt = t;
              }
            }
            if (t >= cycles-step) {
              printf("Final file creation \n");
              if (write_to == 1) {
                pf_count += 1; 
                string avg_voltages_str = "avg_volages_"+to_string(pf_count)+".txt";
                string avg_neur_thresh_str = "avg_neur_thresholds_"+to_string(pf_count)+".txt";
                string f_rate_str = "f_rate_"+to_string(pf_count)+".txt";
          
                avg_voltages.open(avg_voltages_str);
                neur_thresholds.open(avg_neur_thresh_str);
                f_rate.open(f_rate_str);
                float n_avg = 0;
                float total_avg = 0;
                for(int i=0; i<neur_num; i++) {
                  avg_voltages << avg_volts_pos[i]/((t - lt)/step) << endl;
                  neur_thresholds << neur_thresh_pos[i]/((t - lt)/step) << endl;
                  n_avg=float(num_firing[i]);
                  n_avg = n_avg / (t - lt);

                  f_rate << n_avg << endl;
                  total_avg+=n_avg;
                  n_avg=0;
                  avg_volts_pos[i] = 0;
                  neur_thresh_pos[i] = 0;
                  num_firing[i] = 0;
                }
                neur_thresholds.close();
                avg_voltages.close();
                f_rate.close();
              } else {
                nf_count -= 1;
                string avg_voltages_str = "avg_volages_"+to_string(nf_count)+".txt";
                string avg_neur_thresh_str = "avg_neur_thresholds_"+to_string(nf_count)+".txt";
                string f_rate_str = "f_rate_"+to_string(nf_count)+".txt";
          
                avg_voltages.open(avg_voltages_str);
                neur_thresholds.open(avg_neur_thresh_str);
                f_rate.open(f_rate_str);
                float n_avg = 0;
                float total_avg = 0;
                for(int i=0; i<neur_num; i++) {
                  avg_voltages << avg_volts_neg[i]/((t - lt)/step) << endl;
                  neur_thresholds << neur_thresh_neg[i]/((t - lt)/step) << endl;
                  n_avg=float(num_firing[i]);
                  n_avg = n_avg / (t - lt);

                  f_rate << n_avg << endl;
                  total_avg+=n_avg;
                  n_avg=0;

                  avg_volts_neg[i] = 0;
                  neur_thresh_neg[i] = 0;
                  num_firing[i] = 0;
                }
                neur_thresholds.close();
                avg_voltages.close();
                f_rate.close();
              }
            }
            fired_now=0;
            
            for(int i=0; i< neur_num; i++) {
              fired[i]=0;
              neighbor_in[i] = 0;
              // nf_temp[i] = 0;
              // avp_temp[i] = 0;
              // ntp_temp[i] = 0;
              // avn_temp[i] = 0;
              // ntn_temp[i] = 0;
            }
            
          }

          delete[] k1;
          delete[] k2;
          delete[] delta;
          
          
          metric.close();
          act_ex_1.close();
          act_ex_2.close();
        }
      }

    // Outputs:

    }

    pos_times.close();
    neg_times.close();
    delete[] connect;
    delete[] image;
    delete[] voltages;
    delete[] fired;
  }
  

    // deallocation, baby!
  /*
  
*/
  return 0;
};

void init_global(float* connect, float* image, float* voltages, int* fired, int neur_num, int ex1, int ex2, int in1, int in2) {

  // float check_ratio=0;
  // 0<=n < 1000 = Excitatory
  // // 1000 <= n < 2000 = Inhibitory
  
  float j_ee = 1;
  float j_ei = -2;
  float j_ii = -1.8;
  float j_ie = 1;
  float K = 40;
  int neur1 = ex1 + in1;
  int neur2 = ex2 + in2;

  for(int from = 0; from < ex1; from++) {
    for(int target = 0; target < ex1; target++) {
      int seed = rand() % ex1/K;
      if (seed == 1) {
        connect[target*neur_num + from] = float(j_ee)/sqrt(K);
      } else {
        connect[target*neur_num + from] = 0;
      }
    }
  }
  
  for(int from = ex1; from < neur1; from++) {
    for(int target = 0; target < ex1; target++) {
      int seed = rand() % in1/K;
      if (seed == 1) {
        connect[target*neur_num+ from] = float(j_ei)/sqrt(K);
      } else {
        connect[target*neur_num + from] = 0;
      }
    }
  }

  for(int from = ex1; from < neur1; from++) {
    for(int target = ex1; target < neur1; target++) {
      int seed = rand() % in1/K;
      if (seed == 1) {
        connect[target*neur_num + from] = float(j_ii)/sqrt(K);
      } else {
        connect[target*neur_num + from] = 0;
      }
    }
  }

  for(int from = 0; from < ex1; from++) {
    for(int target = ex1; target < neur1; target++) {
      int seed = rand() % ex1/K;
      if (seed == 1) {
        connect[target*neur_num + from] = float(j_ie)/sqrt(K);
      } else {
        connect[target*neur_num + from] = 0;
      }
    }
  }
  // Second A Generation
  for(int from = neur1; from < neur1+ex2; from++) {
    for(int target = neur1; target < neur1+ex2; target++) {
      int seed = rand() % ex2/K;
      if (seed == 1) {
        connect[target*neur_num + from] = float(j_ee)/sqrt(K);
      } else {
        connect[target*neur_num + from] = 0;
      }
    }
  }
  for(int from = neur1+ex2; from < neur_num; from++) {
    for(int target = neur1; target < neur1+ex2; target++) {
      int seed = rand() % in2/K;
      if (seed == 1) {
        connect[target*neur_num + from] = float(j_ei)/sqrt(K);
      } else {
        connect[target*neur_num + from] = 0;
      }
    }
  }
  for(int from = neur1+ex2; from < neur_num; from++) {
    for(int target = neur1+ex2; target < neur_num; target++) {
      int seed = rand() % in2/K;
      if (seed == 1) {
        connect[target*neur_num + from] = float(j_ii)/sqrt(K);
      } else {
        connect[target*neur_num + from] = 0;
      }
    }
  }
  for(int from = neur1; from < neur1+ex2; from++) {
    for(int target = neur1+ex2; target < neur_num; target++) {
      int seed = rand() % ex2/K;
      if (seed == 1) {
        connect[target*neur_num + from] = float(j_ie)/sqrt(K);
      } else {
        connect[target*neur_num + from] = 0;
      }
    }
  }

  for(int from = 0; from < neur1; from++) {
    for(int target = neur1; target < neur_num; target++) {
      connect[target*neur_num + from] = 0;
    }
  }
  for(int from = neur1; from < neur_num; from++) {
    for(int target = 0; target < neur1; target++) {
      connect[target*neur_num + from] = 0;
    }
  }
  for(int from = 0; from < ex1; from++) {
      for(int target = neur1+ex1; target < neur_num; target++) {
        int seed = rand() % ex1/(K);
        if (seed == 1) {
          // connect[target*neur_num + from] = float(j_ie)/sqrt(K);
          connect[target*neur_num + from] = 0;
        } else {
          connect[target*neur_num + from] = 0;
        }
      }
  }
  for(int from = neur1; from < neur1+ex1; from++) {
      for(int target = ex1; target < neur1; target++) {
        int seed = rand() % ex2/(K);
        if (seed == 1) {
          // connect[target*neur_num + from] = float(j_ie)/sqrt(K);
          connect[target*neur_num + from] = 0;
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

  // file read in.
  float f_e = 1.25;
  float f_i = 1;
  //float m_0 = float(1)/float(750);
  // float m_0 = float(1/10.000);
  float m_0 = 0.2;
  float img1_average = float(1/152.000);
  float img2_average = float(1/152.000);

  ifstream image_input_1;
  image_input_1.clear();
  image_input_1.open("im_stripes.txt");
  double temp=0;
  vector<float> IM1;
  while(!image_input_1.eof()) {     //create/fill image vector, IM (c1)
    image_input_1 >> temp;
    IM1.push_back(float(temp));     //vector<double> IM;  //pixel matrix put into vector representation
    // printf("%2.6f\n", float(temp));
  }
  image_input_1.close();

  int img_vec_len = IM1.size() - 1;
  printf("   IM.size = %zu   ", IM1.size());
  float* B1 = new float [neur1*img_vec_len];
  printf("B length: %d\n", neur1*img_vec_len);
  ofstream Bone;
  Bone.open ("B1.txt");
  for(int i=0; i< neur1; i++) {
    for (int j=0; j<img_vec_len; j++) {
      int path = rand() % 2000;
      if (path == 1) {
        B1[i*img_vec_len + j] = 1;
        Bone << "1 ";
      } else {
        B1[i*img_vec_len + j] = 0;
        Bone << "0 ";
      }
    }
    Bone << endl;
  }
  Bone.close();
  
  for(int i = 0; i < ex1; i++) {
    image[i]=0; 
    for(int k = 0; k < img_vec_len; k++) {
      image[i] += B1[i*img_vec_len + k] * IM1[k];
    }
    image[i] = image[i] * f_e * m_0 * img1_average;
  }
  
  for(int i = ex1; i < neur1; i++) {
    image[i]=0; 
    // printf(" image %d, %2.4f  \n", i, image[i]);
    for(int k = 0; k < img_vec_len; k++) {
      image[i] += B1[i*img_vec_len + k] * IM1[k];
    }
    image[i] = image[i] * f_i * m_0 * img1_average;
    // printf(" image %d, %2.4f  \n", i, image[i]);
  }

  ifstream image_input_2;
  image_input_2.clear();
  image_input_2.open("Im_side_stripes.txt");
  temp=0;
  vector<float> IM2;
  while(!image_input_2.eof()) {     //create/fill image vector, IM (c1)
    image_input_2 >> temp;
    IM2.push_back(float(temp));     //vector<double> IM;  //pixel matrix put into vector representation
    // printf("%2.6f\n", float(temp));
  }
  image_input_2.close();

  img_vec_len = IM2.size() - 1;
  printf("   IM.size = %zu   ", IM2.size());
  float* B2 = new float [neur2*img_vec_len];
  printf("B length: %d\n", neur2*img_vec_len);
  ofstream Btwo;
  Btwo.open ("B2.txt");
  
  for(int i=0; i < neur2; i++) {
    // printf("i= %d neur1 = %d neur_num = %d\n", i, neur1, neur_num);
    for (int j=0; j<img_vec_len; j++) {
      int path = rand() % 2000;
      if (path == 1) {
        B2[i*img_vec_len + j] = 1;
        Btwo << "1 ";
      } else {
        B2[i*img_vec_len + j] = 0;
        Btwo << "0 ";
      }
    }
    Btwo << endl;
  }
  Btwo.close();
  printf("init complete \n");
  for(int i = 0; i < ex2; i++) {
    image[i+neur1]=0; 
    for(int k = 0; k < img_vec_len; k++) {
      image[i+neur1] += B2[i*img_vec_len + k] * IM2[k];
    }
    image[i+neur1] = image[i+neur1] * f_e * m_0 * img2_average;
  }

  for(int i = ex2; i < neur2; i++) {
    image[i+neur1]=0; 
    // printf(" image %d, %2.4f  \n", i, image[i]);
    for(int k = 0; k < img_vec_len; k++) {
      image[i+neur1] += B2[i*img_vec_len + k] * IM2[k];
    }
    image[i+neur1] = image[i+neur1] * f_i * m_0 * img2_average;
    // printf(" image %d, %2.4f  \n", i, image[i]);
  }

  
  ofstream fBp;
  fBp.open ("fBp.txt" );
  for(int i=0; i<neur_num; i++) {
    fBp << image[i] << endl;
  }        
  fBp.close();
  //read in image C++
  
  /*
  for(int i = 0; i < neur_num/2; i++) {  // make image input
    image[i] = 1*0.5*sqrt(K);
  }
  for(int i=neur_num/2; i<neur_num; i++) {
    image[i]=0.8*0.5*sqrt(K);
  }
  */

  for(int i = 0; i < neur_num; i++) {   // make neuron voltages
    int seed = rand() % 1000;
    float volt = seed / 1000.0000;
    voltages[i] = volt;
  }

  for(int i =0; i< neur_num; i++) { // set fired to zero
    fired[i] = 0;
  }


}

void chase_spikes(vector<int> to_check, float* voltages, float* delta, int* fired, float* connect, int neur_num, float s, float* threshold, float* neighbor_in) {
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
    /*
    float neighbor_in = 0;
    int target = here;

    for(int from = 0; from < (neur_num); from++) { // sums up fired inputs
      neighbor_in = neighbor_in + connect[target * neur_num + from]
        * fired[from];
    }
    
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
    */
    // neighbors[target] += connect[target * neur_num + from] * fired[from];
    if (voltages[here] + delta[here] + neighbor_in[here] >= threshold[here]) { // check for new firing
      if (fired[here]!=1) {
        fired[here] = 1; // sets as fired
        int from = here;
        // push children for analysis
        for(int target=0; target < neur_num; target++) { // adds anyone who receives output
          if (connect[target * neur_num + from]) {
            mybabies.push_back(target);
            neighbor_in[target] += connect[target * neur_num + from] * fired[from];
          }
        }

      } else {
        printf("%s\n", "refiring in chase_spikes");
      }
    }
  }

  chase_spikes(mybabies, voltages, delta, fired, connect, neur_num, s, threshold, neighbor_in); // recursion
  return;
}