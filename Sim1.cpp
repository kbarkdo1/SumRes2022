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
void make_outputs(int neur_num, float cycles, vector<int>* vec_arr, vector<float> n0_volts, vector<float> n915_volts);

int main() {
  // initialize all variables
  // make neuron connections
  int neur_num = 1000;
  float cycles = 10.00;


  float* connect = new float [neur_num*neur_num];
  float* image = new float [neur_num];
  float* voltages = new float [neur_num];
  int* fired = new int [neur_num];
  init_global(connect, image, voltages, fired, neur_num);

  printf("%2.4f, %2.4f, \n", voltages[0], image[0]);
  printf("%2.4f, %2.4f, \n", voltages[915], image[915]);
  vector<int>* vec_arr = new vector<int> [neur_num];
  vector<float> n0_volts;
  vector<float> n915_volts;
  // iterate through all times steps
  // iterate for 10 by 0.01

  float step = 0.01;
  float f = 1;
  float s = 1;
  printf("all values \n");
  //float B = 1;
  for(float t=0; t<cycles; t = t+step) {
    //printf("Time = %2.2f", t);
    // for every neuron
    float* k1 = new float [neur_num];
    float* k2 = new float [neur_num];
    float* delta = new float [neur_num];

    vector<int> mybabies; // list of neurons that need firing

    for(int n=0; n < neur_num; n++) {
      if (n==0) {
        n0_volts.push_back(voltages[0]); // voltage trace
      }
      if (n==915) {
        n915_volts.push_back(voltages[915]); // voltage trace
      }
      // calculate voltage steps
      float volt_diff = -voltages[n]; // calculate resting voltage difference
      float img_in = image[n];// calculate image input
      k1[n] = volt_diff + img_in;
    }

    for(int n=0; n < neur_num; n++) {
      float volt_diff_2 = -(voltages[n]+k1[n]); // calculate k2 voltage difference
      float img_in_2 = image[n];  // calculate image input
      k2[n] = volt_diff_2 + img_in_2;
    }

    for(int n=0; n < neur_num; n++) {
      delta[n] = 0.5 * (k1[n] + k2[n]);

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

    for(int i=0; i< neur_num; i++) { // add to vector arrays for ouput
      vec_arr[i].push_back(fired[i]);
    }

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
        if (voltages[i]+delta[i]+neighbor_in >= 1) { // failed overrun double check
          printf("unfired overrun, %d, %2.6f, %2.6f\n", i, voltages[i]+delta[i], voltages[i]+delta[i]+neighbor_in);
        } else {
          voltages[i] = voltages[i]+delta[i]+neighbor_in; // forward step
        }
      } else if (fired[i]==1) { // resets spikes
        voltages[i]=0;
        fired[i]=0;
      }
    }
    delete[] k1;
    delete[] k2;
    delete[] delta;
  }

  // Outputs:

  make_outputs(neur_num, cycles, vec_arr, n0_volts, n915_volts);

  // deallocation, baby!
  delete[] connect;
  delete[] image;
  delete[] voltages;
  delete[] fired;
  delete[] vec_arr;
  return 0;
};

void init_global(float* connect, float* image, float* voltages, int* fired, int neur_num) {

  // float check_ratio=0;
  for(int from = 0; from < (neur_num); from++) {
    for(int target = 0; target < (neur_num); target++) {
      int seed = rand() % 10;
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

  for(int i = 0; i < neur_num; i++) {  // make image input
    int seed = rand() % 400;
    seed = 1200-seed;
    float volt = seed / 1000.0000;
    image[i] = volt;
  }

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

void make_outputs(int neur_num, float cycles, vector<int>* vec_arr, vector<float> n0_volts, vector<float> n915_volts) {

  // Raster Plot for MatLab:
  ofstream raster;
  raster.open ("spikes.txt");
  ofstream times;
  times.open ("times.txt");

  for(int i=0; i<neur_num; i+=1) {
    for(int j=0; j < cycles/0.01; j++) {
      if (vec_arr[i][j]) {
        raster << i << endl;
        times << float(j)*0.01 << endl;
        // printf("cycle: %2.6f\n", j);
      }
    }
  }
  raster.close();
  times.close();

  // Boring Voltage Trace:
  ofstream volt_file;
  volt_file.open ("outputvecs/n0_volts.txt" );
  for(int j=0; j< cycles/0.01; j++) {
    volt_file << n0_volts[j] << endl;
  }
  volt_file.close();
  // Cool Voltage Trace
  ofstream volt_file2;
  volt_file2.open ("outputvecs/n915_volts.txt");
  for(int j=0; j< cycles/0.01; j++) {
    volt_file2 << n915_volts[j] << endl;
  }
  volt_file2.close();


  return;
}
