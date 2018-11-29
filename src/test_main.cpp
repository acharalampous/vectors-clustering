#include <iostream>
#include <cstdlib>
#include <fstream>
#include <unistd.h>

#include "dataset.h"
#include "clusters.h"
#include "cl_algorithms.h"

using namespace std;

int main(void){

    ifstream input("twitter_dataset_small_v2.csv");

    cl_management<double> cl_manage(1, 10, 2, 1, 1);
    cl_manage.fill_dataset(input);

    cl_manage.init_clusters();
    cout << "Done initializing clusters" << endl;
   
    cl_manage.assign_clusters();
    cout << "Done assigning to clusters" << endl;

    cl_manage.print();

    cout << "About to update clusters" << endl;
    cl_manage.update_clusters();
    
    cout << "Done updating clusters" << endl;


    cout << "----------------Wait, turn around----------------" << endl;

    cl_manage.assign_clusters();

    cl_manage.evaluation();
    cout << "-----------------------------------" << endl;
}