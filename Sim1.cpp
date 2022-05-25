// Welcome to the Machine

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <unistd.h>

using namespace std;

/*
Overall problems:
Firing is calculated without checking inhibition first. A neuron that is inhibited
won't see that inhibition if it fires from just the input and ICs.

Redundancy: input spiking is calculated twice-as a spiking child
and again for the total updating)(and discarded the first time)

*/


void init_global(float* connect, float* image, float* voltages, int* fired, int neur_num);
void chase_spikes(vector<int> to_check, float* voltages, float* delta, int* fired, float* connect, int neur_num, float s);
float make_outputs(int neur_num, float cycles, vector<int> num_firing); // returns avg firing rate

int main() {
  // initialize all variables
  // make neuron connections
  float forcings[20];
  int tn =926;//  and 928    //TRACKING NEURON
  int tn2=928;
  float start = float(4)/float(1500);
  // float start = 1;
  float force_inc = 0.0;
  ofstream f_vals;
  f_vals.open("forcings.txt");
  for(int i=0; i < 20; i++) {
    forcings[i] = start+force_inc;
    start+=force_inc;
    f_vals << forcings[i] << endl;
    printf("%2.6f\n", forcings[i]);
  }

  ofstream avg_fir;
  avg_fir.open("avg_firing.txt");
  for(int z=0; z<1; z++) {
    ofstream raster;
    raster.open ("spikes.txt");
    ofstream times;
    times.open ("times.txt");
    ofstream volt_file;
    volt_file.open ("tn_volts.txt" );
    ofstream volt_file2;
    ofstream volt_times;
    volt_file2.open ("tn2_volts.txt");
    volt_times.open ("volt_times.txt");
    
    int neur_num = 1000;
    float cycles = 100.00;


    float* connect = new float [neur_num*neur_num];
    float* image = new float [neur_num];
    float* voltages = new float [neur_num];
    int* fired = new int [neur_num];
    init_global(connect, image, voltages, fired, neur_num);

    printf("%2.4f, %2.4f, \n", voltages[tn], image[tn]);
    printf("%2.4f, %2.4f, \n", voltages[tn2], image[tn2]);
    /*
    for(int from = 0; from < (neur_num); from++) {
      if (connect[0 * neur_num + from]) {
        printf("%d \n", from);
      }// sums neighbors if neighbors fire.
    }
    */
    vector<int> num_firing;
    for(int i=0; i< neur_num; i++) {
      num_firing.push_back(0);
    }
    // vector<int>* vec_arr = new vector<int> [neur_num];
    // vector<float> n0_volts;
    // vector<float> n915_volts;
    // iterate through all times steps
    // iterate for 10 by 0.01

    float step = 0.001;
    float f = forcings[z];
    // float f = float(1)/(150*10); // 150 = average image strength * 10 = avg of 10 connections per neurons.
    printf("fL %2.6f", f);
    float s = 1;
    // printf("all values \n");
    //float B = 1;
    for(float t=0; t<cycles; t = t+step) {
      printf("Time = %2.3f", t);
      // for every neuron
      float* k1 = new float [neur_num];
      float* k2 = new float [neur_num];
      float* delta = new float [neur_num];

      vector<int> mybabies; // list of neurons that need firing

      for(int n=0; n < neur_num; n++) {
        if (n==tn) {
          volt_file << voltages[tn] << endl; // voltage trace
        }
        if (n==tn2) {
          volt_file2 << voltages[tn2] << endl; // voltage trace
        }
        // calculate voltage steps
        float volt_diff = -voltages[n]; // calculate resting voltage difference
        float img_in = f * image[n];// calculate image input
        //   k1[n] = volt_diff + img_in;
        k1[n] = step * (volt_diff+img_in);
        if (n==tn) {
          // printf("%s: %2.6f %2.6f %2.6f %2.6f %2.6f\n", "915", k1[n], voltages[n], image[n], img_in, f); // voltage trace
        }
      }

      for(int n=0; n < neur_num; n++) {
        float volt_diff_2 = -(voltages[n]+k1[n]); // calculate k2 voltage difference
        float img_in_2 = f * image[n];  // calculate image input
        //  k2[n] = volt_diff_2 + img_in_2;
        k2[n] = step * (volt_diff_2 + img_in_2);

      }

      for(int n=0; n < neur_num; n++) {
        delta[n] = 0.5 * (k1[n] + k2[n]);
        if (n==tn) {
          //printf("%s: %2.6f\n", "915", delta[n]);; // voltage trace
        }

        if (voltages[n] + delta[n] >= 1) { // checks for spike
          fired[n] = 1; // spiked list
          int from = n; // syntax formality
          for(int target=0; target < neur_num; target++) {
            if (connect[target * neur_num + from]) {
              mybabies.push_back(target);
            }
          }
        }
      }

      chase_spikes(mybabies, voltages, delta, fired, connect, neur_num, s);

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
          /*
           * if (i==tn && neighbor_in !=0 ) {
            printf(" %d receieved %2.4f fron neighs \n", tn, neighbor_in);
          }
          if (i==tn2 && neighbor_in !=0 ) {
            printf(" %d receieved %2.4f fron neighs \n", tn2, neighbor_in);
          }
          */

          if (voltages[i]+delta[i]+neighbor_in >= 1) { // failed overrun double check
            printf("unfired overrun, %d, %2.6f, %2.6f\n", i, voltages[i]+delta[i], voltages[i]+delta[i]+neighbor_in);
          } else {
            voltages[i] = voltages[i]+delta[i]+neighbor_in; // forward step
          }
        } else if (fired[i]==1) { // resets spikes
          raster << i << endl;
          times << t << endl;
          num_firing[i] = num_firing[i] + 1;
          voltages[i]=0;
        }
      }
      for(int i=0; i< neur_num; i++) {
        fired[i]=0;
      }
      delete[] k1;
      delete[] k2;
      delete[] delta;
    }

  // Outputs:

    avg_fir << make_outputs(neur_num, cycles, num_firing) << endl;
    ofstream A;
    A.open ("connectivity.txt" );
    for(int i=0; i< 1000; i++) {
      for(int j=0; j< 1000; j++) {
        A << connect[i*1000 + j] << ' ';
      }
      A << endl;
    }
    A.close();
    
    ofstream fBp;
    fBp.open ("fBp.txt" );
    for(int i=0; i<1000; i++) {
      fBp << image[i] << endl;
    }
    fBp.close();

    raster.close();
    times.close();
    volt_file.close();
    volt_file2.close();
    volt_times.close();





  }

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

void init_global(float* connect, float* image, float* voltages, int* fired, int neur_num) {

  // float check_ratio=0;
  for(int from = 0; from < (neur_num); from++) {
    for(int target = 0; target < (neur_num); target++) {
      int seed = rand() % 100 ;
      if (seed == 1) {
        connect[target*neur_num + from] = 1;
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
  
  // file read in.
  ifstream library_reader;
  library_reader.clear();
  library_reader.open("im_peppers.txt");
  double temp=0;
  vector<float> IM;
  while(!library_reader.eof()) {     //create/fill image vector, IM (c1)
    library_reader >> temp;
    IM.push_back(float(temp));     //vector<double> IM;  //pixel matrix put into vector representation
    printf("%2.6f\n", float(temp));
  }
  library_reader.close();
  int img_vec_len = IM.size() - 1;
  printf("   IM.size = %zu   ", IM.size());
  float img_vec[img_vec_len][1];
  for(int i=0; i<img_vec_len; i++) {
    img_vec[i][0] = IM[i];
    // printf(" img_vec %d: %4.4f \n", i, img_vec[i][0]);
  }
  float* B = new float [neur_num*img_vec_len];
  printf("IM length: %d\n", int(IM.size()));
  printf("B length: %d\n", neur_num*img_vec_len);
  
  ofstream Be;
  Be.open ("B.txt");
  int counter=0;
  for(int i=0; i< neur_num; i++) {
    for (int j=0; j<img_vec_len; j++) {
      int path = rand() % 1000;
      if (path == 1) {
        // all this stuff below is for dynamic connection strengths
        /* 
        int seed = rand() % 40;
        seed = 120-seed;
        float volt = seed / (150);
        // printf("index: %d ", i*img_vec_len + j);
        B[i*img_vec_len + j] = volt;
        */
        B[i*img_vec_len + j] = 1;
        Be << "1 ";
        counter += 1;
      } else {
        B[i*img_vec_len + j] = 0;
        Be << "0 ";
      }

    }
    printf(" counter: row %d: %d \n", i, counter); 
    counter=0;
    Be << endl;
  }
  Be.close();

  for(int i = 0; i < neur_num; i++) {
    image[i]=0; 
    printf(" image %d, %2.4f  \n", i, image[i]);
    for(int k = 0; k < 10000; k++) {
      image[i] += B[i*neur_num + k] * img_vec[k][0];
    }
    printf(" image %d, %2.4f  \n", i, image[i]);
 }
  
  //read in image C++
  /*
  
  for(int i = 0; i < neur_num; i++) {  // make image input
    srand(i);
    int seed = rand() % 2000;
    seed = 5000-seed;
    float volt = seed / 1000.0000;
    image[i] = volt;
    if (i ==599 || i==609 || i==617 || i==622 || i==624 || i==628 || i==637 || i==641 || i==616 || i==620 || i==629 ) {
      printf("\n %d image: %2.4f", i, image[i]);
    }

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

void chase_spikes(vector<int> to_check, float* voltages, float* delta, int* fired, float* connect, int neur_num, float s) {
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

    if (voltages[here] + delta[here] + neighbor_in >= 1) { // check for new firing
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

  chase_spikes(mybabies, voltages, delta, fired, connect, neur_num, s); // recursion
  return;
}

float make_outputs(int neur_num, float cycles, vector<int> num_firing) {
  printf("%s\n", "making ouputs");
  float n_avg=0;
  float total_avg = 0;
  ofstream f_rate;
  f_rate.open("f_rate.txt");
  // printf("\n %2.4f %2.4f \n", cycles, (cycles)); 
  for(int i=0; i<neur_num; i+=1) {
    n_avg=float(num_firing[i]);
    if(i==213) {
      printf("213: %2.4f, %2.4f", n_avg, n_avg/cycles);
    }
    if(i==215) {
      printf("215: %2.4f, %2.4f", n_avg, n_avg/cycles);
    }
    n_avg = n_avg / cycles;
    if(i==213) { printf(" 213: %2.4f", n_avg); }
    if(i==215) { printf(" 215: %2.4f", n_avg); }
    if(i==23) { printf(" 23: %2.4f", n_avg); }
    if(i==29) { printf(" 29: %2.4f", n_avg); }
    if(i==39) { printf(" 39: %2.4f", n_avg); }

    f_rate << n_avg << endl;
    total_avg+=n_avg;
    n_avg=0;
    // printf("total_avg: %d %2.4f  ", i, total_avg);
  }
  f_rate.close();
  // printf("%2.6f %d \n", total_avg, neur_num);
  total_avg = total_avg / (neur_num);
  return total_avg;
}
