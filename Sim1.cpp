// Welcome to the Machine

#include <iostream>
#include <string>
#include <vector>
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <unistd.h>

using namespace std;

int main() {
  // initialize all variables
  // make neuron connections
  int* connect = new int [1000*1000];

  float check_ratio=0;
  for(int from = 0; from < (1000); from++) {
    for(int target = 0; target < (1000); target++) {
      int seed = rand() % 10;
      if (seed == 1) {
        connect[target*1000 + from] = 1;
        // printf(" %d", connect[target*1000 + from]);
        check_ratio++;
      }
    }
  }
  // check_ratio = check_ratio/(1000*1000);
  // printf("   %2.6f   ", check_ratio);

  // make image input
  float* image = new float [1000];

  for(int i = 0; i < 1000; i++) {
    int seed = rand() % 400;
    seed = 1200-seed;
    float volt = seed / 1000.0000;
    image[i] = volt;
    // printf(" %2.4f", image[i]);
  }
  // make neuron voltages
  float* voltages = new float [1000];

  for(int i = 0; i < 1000; i++) {
    int seed = rand() % 1000;
    float volt = seed / 1000.0000;
    voltages[i] = volt;
    // printf(" %2.4f", voltages[i]);
  }




  // iterate through all times steps
  // iterate for 10 by 0.01
  float step = 0.01;
  for(float t=0; t<10; t = t+step) {

    // for every neuron
    float* k1 = new float [1000];
    float* k2 = new float [1000];
    float* delta = new float [1000];
    int fired[1000];
    vector<int> mybabies; // list of neurons that need firing

    for(int n=0; n < 1000; n++) {
      // calculate voltage steps
      // calculate resting voltage difference
      float volt_diff = -voltages[n];
      printf("voltage_diff = %2.4f\n", volt_diff);
      // calculate image input
      float img_in = image[n];
      printf("img_in = %2.4f\n", img_in);
      // calculate neighboring inputs

      k1[n] = volt_diff + img_in;

      float volt_diff_2 = -(voltages[n]+k1[n]);
      printf("voltage_diff = %2.4f\n", volt_diff);
      // calculate image input
      float img_in_2 = image[n];
      printf("img_in = %2.4f\n", img_in);
      // calculate neighboring inputs

      k2[n] = volt_diff_2 + img_in_2;

      delta[n] = 0.5 * (k1[n] + k2[n]);

      /*
      int target = n;
      for(int from = 0; from < (1000); from++) {
        neighbor_in = neighbor_in + connect[target * 1000 + from]
          * voltages[from];
        //printf("neighbor_in = %2.4f\n", neighbor_in);
      }
      int s = 1;
      int nn = 1000;
      float s_N = s / float(nn);
      neighbor_in = neighbor_in * s_N;
      printf("neighbor_in = %2.10f\n", neighbor_in);
      // sleep(5);
      */
      // check overrun - note k1 is just the forward step

      if (voltages[n] + delta[n] >= 1) {
        voltages[n] = 0;
        fired[n] = 1;
        int from = n;
        for(int target=0; target < 1000; target++) {
          if (connect[target * 1000 + from]) {
            mybabies.push_back(target);
          }
        }
      }

      // iterating through mybabies - note k1 is just the forward step

    }

    for(int i=0; i< mybabies.size(); i++) {
      // recalculating action potential
      if (!fired[mybabies[i]]) { // i f its already fired
        continue;
      }
      int here = mybabies[i];
      float neighbor_in = 0;
      int target = here;
      for(int from = 0; from < (1000); from++) {
        // note this is not how to calcualte neighbors if only action
        // potential neighbors count
        neighbor_in = neighbor_in + connect[target * 1000 + from]
          * fired[from];
      }
      int s = 1;
      int nn = 1000;
      float s_N = s / float(nn);
      neighbor_in = neighbor_in * s_N;
      delta[here] = delta[here] + neighbor_in; // note no RK2 on this
      if (voltages[here] + delta[here] >= 1) {
        voltages[here] = 0;
        fired[here] = 1;
        int from = here;
        for(int target=0; target < 1000; target++) {
          if (connect[target * 1000 + from]) {
            mybabies.push_back(target);
          }
        }
      }
    }

    /*

    // calculate the second steps
    
    for rec
    int s = 1;
    int nn = 1000;
    float s_N = s / float(nn);
    neighbor_in = neighbor_in * s_N;
    k1[rec] = k1[rec] + neighbor_in;
    if k1[rec]>=1 {

    }
    */
    // take average voltage step w/ heun

    // check overrun to 1
    // if 1
      // record as spiked (hashtable?)
      // note connections
      // for each connection
        // update with new input
        // check overrun to 1
        // if 1
          // this seems like a great place for recursion, carefully
          // or a list
      // save data
      // save which neurons fired
      // save selective single neuron voltages
  }
  return 0;
};
