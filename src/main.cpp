#include <iostream>
#include <cstdlib>
#include <fstream>
#include <unistd.h>

#include "dataset.h"
#include "clusters.h"
#include "cl_algorithms.h"

using namespace std;

int main(int argc, char* argv[]){
    exe_args parameters;

    int result = get_parameters(argc, argv, parameters);
    if(result == -2){
        printValidParameters();
        return -1;
    }

    ofstream output;
    result = validate_parameters(parameters, output);
    if(result == -1){
        return -1;
    }
    else if(result == -2){
        printValidParameters();
        return -1;
    }
    else if(result == -3){
        printValidConfig();
        return -2;
    }

    int choice = 1; 
    int init;
    int assign;
    int upd;
    if(parameters.all_combinations == 0){

        choice = read_combination(init, assign, upd);
    }

    if(choice == 0){
        cl_management<double>* cl_manage = new cl_management<double>(parameters, init, assign, upd);   
        cl_manage->clustering(parameters, output, init, assign, upd);    
        delete cl_manage;
    }
    else if(choice == 1){
        for(int init = 1; init <= 2; init++){
            for(int assign = 1; assign <= 3; assign++){
                for(int upd = 1; upd <= 2; upd++){
                    cl_management<double>* cl_manage = new cl_management<double>(parameters, init, assign, upd);   
                    cl_manage->clustering(parameters, output, init, assign, upd);    
                    delete cl_manage;
                }
            }
        }
    }
    else if(choice == 2){
        cout << "Choosed to exit! Abort." << endl;
        return 1;
    }


    cout << "DONE" << endl;
    return 0;
}