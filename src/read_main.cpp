#include <iostream>
#include <cstdlib>
#include <fstream>
#include <malloc.h>

#include "dataset.h"

using namespace std;

int main(void){

    ifstream new_file("twitter_dataset_small_v2.csv");
    string line;
    string id;

    dataset<double> new_dataset;



    while(getline(new_file, line)){
        new_dataset.add_vector(line);
    }

    
}