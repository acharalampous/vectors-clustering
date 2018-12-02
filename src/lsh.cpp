/*******************************/
/* lsh.cpp */

/* Name:    Andreas Charalampous
 * A.M :    1115201500195
 * e-mail:  sdi1500195@di.uoa.gr
 */
/********************************/
#include <iostream>
#include <unordered_set>

#include "lsh.h"
#include "metrics.h"

using namespace std;

template class LSH<int>;
template class LSH<double>;

/*  Implementation of all functions of the class
 *  that is used for LSH. Definitions found in
 *  lsh.h.
 */

/* L: num of tables for each metric, k: num of hash functions, */
/* n: dataset size */ 
template <class T>
LSH<T>::LSH(int metrics, int L, int k, int n){
    this->L = L;
    if(metrics == 1){ // euclidean is defined
        for(int i = 0; i < L; i++)
            eu_tables.push_back(new euclidean<T>(k, n));
    }

    else if(metrics == 2){ // cosine is defined
        for(int i = 0; i < L; i++)
            cs_tables.push_back(new csimilarity<T>(k));
    }
}

template <class T>
LSH<T>::~LSH(){
    for(unsigned int i = 0; i < eu_tables.size(); i++){
        delete eu_tables[i];
    }

    eu_tables.clear();

    for(unsigned int i = 0; i < cs_tables.size(); i++){
        delete cs_tables[i];
    }

    cs_tables.clear();
}

template <class T>
void LSH<T>::add_vector(vector_item<T>* new_vector){
    for(unsigned int i = 0; i < eu_tables.size(); i++){
        eu_tables[i]->add_vector(new_vector);
    }

    for(unsigned int i = 0; i < cs_tables.size(); i++){
        cs_tables[i]->add_vector(new_vector);
    }
}

// template <class T>
// void LSH<T>::assign_clusters(cl_management<T>& cl_manage){
//     vector<cluster<T>*> clusters = cl.get_clusters();
//     int num_of_centroids = clusters.size(); // get number of centroids to be checked
//     vector<unordered_set<int>> checked_set; // set that holds items that were checked

//     vector<vector<vector_check<T>*>> vectors_to_check; // avoid finding vectors and recalculating same distances 
//     vector<cluster_info*> vectors_info = cl.get_vectors_info();

//     int vectors_left = vectors_info.size();

//     for(int i = 0; i < num_of_centroids; i++){
//         checked_set.push_back(unordered_set<string>);
//         vectors_to_check.push_back(vector<vector_item<T>*>);

//         vectors_left -= clusters[i]->get_centroid_type();
//         buckets_left.push_back(this->L);
//     }

//     double r = get_starting_r(cl, cl.get_dist_func());

//     for(int i = 0; i < num_of_centroids; i++){
//         for(int j = 0; j < eu_tables.size(); j++){
//             int res = eu_tables[j]->first_assign(clusters, i, r, checked_set[i], vectors_to_check[i], vectors_info, vectors_left);
//         }
//     }

// }

template <class T>
void LSH<T>::findANN(vector_item<T>& query, float radius, float& min_dist, string& ANN_name, ofstream& output){
    unordered_set<string> checked_set; // set that holds items that were checked
    for(unsigned int i = 0; i < eu_tables.size(); i++){
        eu_tables[i]->findANN(query, radius, min_dist, ANN_name, output, checked_set);
    }

    for(unsigned int i = 0; i < cs_tables.size(); i++){
        cs_tables[i]->findANN(query, radius, min_dist, ANN_name, output, checked_set);
    }
}

template <class T>
long int LSH<T>::get_total_size(){
    long int total_size = 0;

    total_size += sizeof(*this);
    for(unsigned int i = 0; i < eu_tables.size(); i++){
        total_size += eu_tables[i]->get_size();
    }

    for(unsigned int i = 0; i < cs_tables.size(); i++){
        total_size += cs_tables[i]->get_size();
    }

    return total_size;
}