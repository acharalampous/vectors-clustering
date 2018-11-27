/*******************************/
/* cl_algorithms.cpp */

/* Name:    Andreas Charalampous
 * A.M :    1115201500195
 * e-mail:  sdi1500195@di.uoa.gr
 */
/********************************/
#include <iostream>
#include <random>
#include <array>

#include "cl_algorithms.h"
using namespace std;

template class cl_init_random<double>;
template class cl_init_kmeans<double>;




////////////////////
// CL_INIT_RANDOM //
////////////////////
template <class T>
void cl_init_random<T>::init_clusters(cl_management<T>& cl_manage){
    random_device rd;
    mt19937 gen(rd()); // choose random vector as centroids
        
    dataset<T>* all_vectors = cl_manage.get_dataset();
    vector<cluster_info*>& vectors_info = cl_manage.get_vectors_info();
    vector<cluster<T>*>& clusters = cl_manage.get_clusters();


    int k = cl_manage.get_k(); // get number of total clusters

    int num_of_vectors = all_vectors->get_counter(); // get num of vectors

    uniform_int_distribution<int> rand_Z(0, num_of_vectors - 1);

    for(int i = 0; i < k; i++){
        int sel = rand_Z(gen); // select random vector
        
        if(vectors_info[sel]->get_is_centroid() == 1){ // vector already a centroid
            i--;
            continue;
        }
    
        /* Set point as centroid */
        clusters[i]->set_centroid(all_vectors->get_item(sel));
        clusters[i]->set_centroid_type(1);
        
        vectors_info[sel]->set_centroid();
        vectors_info[sel]->set_cluster(i);
    }
}



////////////////////
// CL_INIT_KMEANS //
////////////////////
template <class T>
void cl_init_kmeans<T>::init_clusters(cl_management<T>& cl_manage){
    random_device rd;
    mt19937 gen(rd()); // choose random vector as centroids
        
    dataset<T>* all_vectors = cl_manage.get_dataset();
    vector<cluster_info*>& vectors_info = cl_manage.get_vectors_info();
    vector<cluster<T>*>& clusters = cl_manage.get_clusters();
    
    dist_func dist_function;
    dist_function = cl_manage.get_dist_func();

    int k = cl_manage.get_k(); // get number of total clusters
    int num_of_vectors = all_vectors->get_counter(); // get num of vectors

    /* The first centroid is chosen randomly */
    uniform_int_distribution<int> rand_Z(0, num_of_vectors - 1);
    int sel = rand_Z(gen); // select random first centroid

    /* Set point as centroid */
    clusters[0]->set_centroid(all_vectors->get_item(sel));
    clusters[0]->set_centroid_type(1);
    
    vectors_info[sel]->set_centroid();
    vectors_info[sel]->set_cluster(0);

    /* Array that maps the minimum distance of each vector to all centroids */
    vector<T> min_distances;
    for(int i = 0; i < num_of_vectors; i++) // Initialize array
        min_distances.push_back(0.0);
    
    T max_distance = 0.0; // distance of the furthest vector

    vector<dist_mapping> probs_array; // array that holds probabilities

    for(int i = 1; i < k; i++){
        probs_array.clear();
        max_distance = 0.0;
        vector_item<T>* last_centroid = clusters[i - 1]->get_centroid();

        /* First iteration, so wont check for min_distance(minimize checks) */
        if(i == 1){
            
            /* Find minimum distance to centroids. Must check only with the last centroid */
            for(int j = 0; j < num_of_vectors; j++){

                /* Check if not centroid */
                if(vectors_info[j]->get_is_centroid() == 1){
                    continue;
                }

                vector_item<T>* curr_item = all_vectors->get_item(j);

                /* Set minimum distance */
                min_distances[j] = dist_function(*curr_item, *last_centroid);
                if(min_distances[j] > max_distance)
                    max_distance = min_distances[j];
            
                probs_array.push_back(dist_mapping(curr_item->get_index(), min_distances[j]));

            }
        } // end if i == 1
        else{
            /* Find minimum distance to centroids. Must check only with the last centroid */
            for(int j = 0; j < num_of_vectors; j++){
                /* Check if not centroid */
                if(vectors_info[j]->get_is_centroid() == 1){
                    continue;
                }

                vector_item<T>* curr_item = all_vectors->get_item(j);

                double dist = dist_function(*curr_item, *last_centroid);
                if(dist < min_distances[j]){
                    min_distances[j] = dist;

                    if(min_distances[j] > max_distance)
                        max_distance = min_distances[j];
                }

                probs_array.push_back(dist_mapping(curr_item->get_index(), min_distances[j]));
            }
        } // end else i != 1
        
        /* Normalize distances */
        probs_array[0].distance = ((double)pow(probs_array[0].distance, 2)) / max_distance;
        
        for(int j = 1; j < num_of_vectors - i; j++){
            probs_array[j].distance = (probs_array[j - 1].distance + ((double)pow(probs_array[j].distance, 2))) / max_distance;
        }
    }
}

template <class T>
int cl_init_kmeans<T>::get_next_centroid(vector<dist_mapping>& probs_array){
    random_device rd;
    mt19937 gen(rd());
    
    int start = 0;
    int end = probs_array.size() - 1;

    double min_val = 0.0;
    double max_val = probs_array[end].distance;

    uniform_real_distribution<> rand_R(min_val, max_val);

    double val = rand_R(gen);

    /* perform binary search for given val */
    while(1){
        int pivot = (end + start) / 2;
        
        if(val == probs_array[pivot].distance)
            return pivot;
        else if(val >= probs_array[pivot + 1].distance){
            start = pivot;
            continue;
        }
        else if(val < probs_array[pivot].distance){
            end = pivot;
            continue;
        }
        else if(val > probs_array[pivot].distance && val < probs_array[pivot + 1].distance)
            return pivot + 1;

        cout << "[BINARY SEARCH]I SHOULD NOT BE HERE" << endl;
    }
}