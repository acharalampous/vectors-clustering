#include <iostream>
#include <cstdlib>
#include <fstream>
#include <unistd.h>

#include "dataset.h"
#include "clusters.h"
#include "cl_algorithms.h"

using namespace std;

int main(void){
    int i = 0;
    ifstream input("twitter_dataset_small_v2.csv");

    int metric = 2;
    int k = 40;
    int max_updates = 5;
    int L = 4;
    int hf = 5;
    int hc_probes = 3;
    int hc_M = 100;
    int init = 1;
    int assign = 1;
    int upd = 1;


    cl_management<double> cl_manage(metric, k, L, 1, max_updates, hf, hc_probes, hc_M, init, assign, upd);

    cl_manage.fill_dataset(input);

    cl_manage.tick();
    cl_manage.init_clusters();
    
    cout << "Done initializing clusters" << endl;
   
    cl_manage.assign_clusters();
    cout << "Done assigning to clusters" << endl;

    while(1){
        int made_changes = cl_manage.update_clusters();
        cl_manage.assign_clusters();
        cout << "Iteration #" << i++ << ":" << endl;
        if(made_changes == 0 || i >= max_updates)
            break;
    }
    cl_manage.tock();
    cl_manage.silhouette();
    cl_manage.print_to_file();
    cout << "-----------------------------------" << endl;
}