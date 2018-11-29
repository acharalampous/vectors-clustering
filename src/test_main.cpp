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

    cl_management<double> cl_manage(1, 4, 1, 1, 1);
    cl_manage.fill_dataset(input);

    cl_manage.init_clusters();
    cout << "Done initializing clusters" << endl;
   
    cl_manage.assign_clusters();
    cout << "Done assigning to clusters" << endl;

    while(1){
        int res = cl_manage.update_clusters();
        cl_manage.assign_clusters();
        cout << "Iteration #" << i++ << ":" << endl;
        if(res == 0)
            break;
        
    }

    cl_manage.evaluation();
    cout << "-----------------------------------" << endl;
}