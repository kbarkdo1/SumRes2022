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



void chase_spikes(vector<int> to_check, float* voltages, float* delta, int* fired, int* connect, int neur_num);

int main() {
  // initialize all variables
  // make neuron connections
  int neur_num = 1000;
  float cycles = 10.00;
  int* connect = new int [neur_num*neur_num];

  float check_ratio=0;
  for(int from = 0; from < (neur_num); from++) {
    for(int target = 0; target < (neur_num); target++) {
      int seed = rand() % 10;
      if (seed == 1) {
        connect[target*neur_num + from] = 1;
        // printf(" %d", connect[target*neur_num + from]);
        check_ratio++;
      }
    }
  }
  // set horiztonal to zero, no reccurence
  for(int i=0; i < neur_num; i++) {
    connect[i*neur_num + i] = 0;
  }
  // check_ratio = check_ratio/(neur_num*neur_num);
  // printf("   %2.6f   ", check_ratio);

  // make image input
  float* image = new float [neur_num];

  for(int i = 0; i < neur_num; i++) {
    int seed = rand() % 400;
    seed = 1200-seed;
    float volt = seed / 1000.0000;
    image[i] = volt;
    // printf(" %2.4f", image[i]);
  }
  // make neuron voltages
  float* voltages = new float [neur_num];

  for(int i = 0; i < neur_num; i++) {
    int seed = rand() % 1000;
    float volt = seed / 1000.0000;
    voltages[i] = volt;
    // printf(" %2.4f", voltages[i]);
  }

  int* fired = new int [neur_num];
  for(int i =0; i< neur_num; i++) {
    fired[i] = 0;
  }


  printf("%2.4f, %2.4f, \n", voltages[0], image[0]);
  printf("%2.4f, %2.4f, \n", voltages[915], image[915]);
  vector<int>* vec_arr = new vector<int> [neur_num];
  vector<float> n0_volts;
  vector<float> n915_volts;
  // iterate through all times steps
  // iterate for 10 by 0.01

  float step = 0.01;
  float f = 1;
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
        n0_volts.push_back(voltages[0]);
      }
      if (n==915) {
        n915_volts.push_back(voltages[915]);
      }
      // calculate voltage steps
      // calculate resting voltage difference
      float volt_diff = -voltages[n];
      //printf("voltage_diff = %2.4f\n", volt_diff);
      // calculate image input
      float img_in = image[n];
      //printf("img_in = %2.4f\n", img_in);
      // calculate neighboring inputs

      k1[n] = volt_diff + img_in;
    }

    for(int n=0; n < neur_num; n++) {
      float volt_diff_2 = -(voltages[n]+k1[n]);
      //printf("voltage_diff = %2.4f\n", volt_diff);
      // calculate image input
      float img_in_2 = image[n];
      // float img_in_2 = f * B[n] * image[n]; // for future details
      //printf("img_in = %2.4f\n", img_in);
      // calculate neighboring inputs

      k2[n] = volt_diff_2 + img_in_2;
    }

    for(int n=0; n < neur_num; n++) {
      delta[n] = 0.5 * (k1[n] + k2[n]);

      /*
      int target = n;
      for(int from = 0; from < (neur_num); from++) {
        neighbor_in = neighbor_in + connect[target * neur_num + from]
          * voltages[from];
        //printf("neighbor_in = %2.4f\n", neighbor_in);
      }
      int s = 1;
      int nn = neur_num;
      float s_N = s / float(nn);
      neighbor_in = neighbor_in * s_N;
      printf("neighbor_in = %2.10f\n", neighbor_in);
      // sleep(5);
      */
      // check overrun - note k1 is just the forward step

      if (voltages[n] + delta[n] >= 1) {
        //printf("We have a firing!");
        //voltages[n] = 0;
        fired[n] = 1;
        int from = n;
        for(int target=0; target < neur_num; target++) {
          if (connect[target * neur_num + from]) {
            mybabies.push_back(target);
          }
        }
      }

      // iterating through mybabies - note k1 is just the forward step

    }
    chase_spikes(mybabies, voltages, delta, fired, connect, neur_num);



    for(int i=0; i< neur_num; i++) {
      vec_arr[i].push_back(fired[i]);
      /*if (i==915) {
        printf("915: %d %d %2.6f %2.6f \n", fired[i], vec_arr[i][vec_arr[i].size()-1], voltages[i], t );
      }*/
      /*if (i==915) {
        printf("915: %d %d   ", fired[i], vec_arr[i][vec_arr[i].size()-1]);
      }*/
    }
    for(int i=0; i< neur_num; i++) {

      if (fired[i]==0) {
        float neighbor_in = 0;
        int target = i;
        for(int from = 0; from < (neur_num); from++) {
          // note this is not how to calcualte neighbors if only action
          // potential neighbors count
          neighbor_in = neighbor_in + connect[target * neur_num + from]
            * fired[from];
        }
        int s = 1;
        int nn = neur_num;
        float s_N = s / float(nn);
        neighbor_in = neighbor_in * s_N;
        if (voltages[i]+delta[i]+neighbor_in >= 1) {
          printf("unfired overrun, %d, %2.6f, %2.6f\n", i, voltages[i]+delta[i], voltages[i]+delta[i]+neighbor_in);
        } else {
          voltages[i] = voltages[i]+delta[i]+neighbor_in;
        }
      } else if (fired[i]==1) {

        voltages[i]=0;
        fired[i]=0;
        //if (voltages[i] != 0) {
          //printf("error: fired neuron has nonzero voltage");
        //}
      }
    }
    /*

    // calculate the second steps

    for rec
    int s = 1;
    int nn = neur_num;
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
  /*
  //for n output files
  for(int i=0; i<neur_num; i+=1) {
    ofstream file;
    float k = float(i)/neur_num;
    file.open ("outputvecs/"+ to_string(k) + ".txt");
    for(int j=0; j< 10/0.01; j++) {
      file << vec_arr[i][j] << endl;
    }
    file.close();

  }
  */
  ofstream n915;
  n915.open ("n915_fire_track.txt");
  for(int j=0; j< cycles/0.01; j++) {
    printf("%d\n", vec_arr[915][j]);
    n915 << vec_arr[915][j];
  }
  n915.close();

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


  ofstream volt_file;
  volt_file.open ("outputvecs/n0_volts.txt" );
  for(int j=0; j< cycles/0.01; j++) {
    volt_file << n0_volts[j] << endl;
  }
  volt_file.close();
  ofstream volt_file2;
  volt_file2.open ("outputvecs/n915_volts.txt");
  for(int j=0; j< cycles/0.01; j++) {
    volt_file2 << n915_volts[j] << endl;
  }
  volt_file2.close();
  return 0;
};


void chase_spikes(vector<int> to_check, float* voltages, float* delta, int* fired, int* connect, int neur_num) {
  // printf("chasing_spikes with %zu \n", to_check.size());
  // printf("voltage 0: %2.4f ", voltages[915]);
  /*int connect_sum=0;
  for(int i=0;i<1000*1000;i++) {
    connect_sum += connect[i];
  }
  printf("connect sum: %d ", connect_sum);*/
  vector<int> mybabies;
  int num_fi = 0;
  if (to_check.size()< 1) {
    return;
  }
  /*
  for (int i=0; i<1000; i++) {
    printf("%d %2.4f  ", i, voltages[i]);
  }
  */
  for(int i=0; i<to_check.size(); i++) {
    // printf("child: %d %2.6f %2.6f \n", to_check[i], voltages[to_check[i]], voltages[to_check[i]] + delta[to_check[i]]);
    // recalculating action potential
    if (fired[to_check[i]]) { // i f its already fired
      num_fi++;
      continue;
    }
    int here = to_check[i];
    float neighbor_in = 0;
    int target = here;

    for(int from = 0; from < (neur_num); from++) {
      // note this is not how to calcualte neighbors if only action
      // potential neighbors count
      neighbor_in = neighbor_in + connect[target * neur_num + from]
        * fired[from];
        /*
        if (neighbor_in != 0) {
          printf("  %2.4f,%2.4f,%2.4f,%2.4f ", voltages[i], delta[i], neighbor_in, voltages[here] + delta[here] + neighbor_in);
        }
        */
    }

    int s = 1;
    int nn = neur_num;
    float s_N = s / float(nn);
    neighbor_in = neighbor_in * s_N;

    // printf("neighbors_in: %2.4f", neighbor_in);
    /* wait this doesn't work. Right now if it has five neurons firing,
    it might count up three the first time, then all five the second time,
    and approximate a total of 8 inputs, because it recounted the three the
    second time around. if 5 neurons point into it, and 4 fire, it will be
    in mybabies 4 times, maybe this is the solution */
    /*
    for(int from = 0; from < (neur_num); from++) {
      // note this is not how to calcualte neighbors if only action
      // potential neighbors count
      neighbor_in = neighbor_in + connect[target * neur_num + from]
        * fired[from];
    }
    int s = 1;
    int nn = neur_num;
    float s_N = s / float(nn);
    neighbor_in = neighbor_in * s_N;
    */
    // attempt 2

    // delta[here] = delta[here] + s_N; // one more tiny connection
    // delta[here] = delta[here] + neighbor_in; // note no RK2 on this
    // printf("future voltage: %d %2.4f, ", here, voltages[here] + delta[here] + neighbor_in);

    // printf("child: %d %2.4f %2.4f ", here, voltages[here], (1-voltages[here]-delta[here]) - neighbor_in);
    /*
    if (((1-voltages[here]-delta[here]) - neighbor_in) <= 0) {
      printf("CHECKCHEKC ");
    }
    */
    //printf("Child %d: %2.6f, %2.6f, %2.6f\n", here, voltages[here] + delta[here], neighbor_in, voltages[here] + delta[here] + neighbor_in);
    /*if (here==893) {
      printf("expo 893: %d %2.6f %2.6f \n", here, voltages[here] + delta[here], voltages[here] + delta[here] + neighbor_in );
    }*/
    if (voltages[here] + delta[here] + neighbor_in >= 1) {
      if (fired[here]!=1) {
        //  printf("<-%d spiked \n", here);
        //printf("child spike");
        // voltages[here] = 0;
        fired[here] = 1;
        int from = here;
        // push children for analysis
        for(int target=0; target < neur_num; target++) {
          if (connect[target * neur_num + from]) {
            mybabies.push_back(target);
          }
        }

      }

    }
  }
  // printf("num_fi: %d ", num_fi);
  chase_spikes(mybabies, voltages, delta, fired, connect, neur_num);
  return;
}
