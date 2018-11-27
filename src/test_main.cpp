#include <iostream>
#include <cstdlib>
#include <fstream>
#include <malloc.h>

#include "dataset.h"
#include "clusters.h"
#include "cl_algorithms.h"

using namespace std;

int main(void){

    ifstream input("twitter_dataset_small_v2.csv");

    cl_management<double> cl_manage(1, 10, 2);
    cl_manage.fill_dataset(input);

    cl_manage.init_clusters();

    cout << "Wait, turn around" << endl;

    cl_manage.print();
}