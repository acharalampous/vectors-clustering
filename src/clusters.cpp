/*******************************/
/* clusters.cpp */

/* Name:    Andreas Charalampous
 * A.M :    1115201500195
 * e-mail:  sdi1500195@di.uoa.gr
 */
/********************************/
#include <iostream>
#include <fstream>

#include "utils.h"
#include "clusters.h"

using namespace std;

template class cluster<double>;
template class cl_management<double>;


/*  Implementation of all functions of cluster_info
 *  Definitions found in clusters.h.
 */


//////////////////////
//** CLUSTER_INFO **//
//////////////////////
cluster_info::cluster_info(int index){ 
    is_centroid = 0;
    cluster_num = -1;
    this->index = index;
    distance = -1.0;
}

/* Returns is_centroid */
int cluster_info::get_is_centroid(){
    return this->is_centroid;
}

/* Returns cluster_num */
int cluster_info::get_cluster_num(){
    return this->cluster_num;
}

/* Returns index of vector in dataset */
int cluster_info::get_index(){
    return this->index;
}

/* Returns the distance of vector to centroid */
double cluster_info::get_distance(){
    return this->distance;
}

/* Set is_centroid as 1 */
void cluster_info::set_centroid(){
    this->is_centroid = 1;
}

/* Set cluster_num as the given int */
void cluster_info::set_cluster(int new_cluster){
    this->cluster_num = new_cluster;
}

/* Set index as the given int */
void cluster_info::set_index(int new_index){
    this->index = 1;
}

/* Set distance of vector to centroid */
void cluster_info::set_distance(double new_dist){
    this->distance = new_dist;
}


/* Unassign vector from cluster, setting cluster_num to -1 */
void cluster_info::reset_info(){
    this->is_centroid = 0;
    this->cluster_num = -1;
}



/////////////
// CLUSTER //
/////////////
template <class T>
cluster<T>::cluster(int cluster_num){
    this->cluster_num = cluster_num;
    this->centroid_type = -1;
    this->centroid = NULL;
}

template <class T>
void cluster<T>::add_vector(vector_item<T>* new_vector){
    this->vectors.push_back(new_vector);

}


template <class T>
void cluster<T>::evaluation(vector<cluster<T>*>& clusters, dist_func& dist){
    int num_of_vectors = vectors.size() + 1;
    double* s_values;
    double s_of_cluster = 0.0;

    /* Must not perform silhouette for centroid if it is not in dataset */
    if(this->get_centroid_type() == 0){
        num_of_vectors--;
        s_values = new double[num_of_vectors];

        /* Silhouettes Values */
        for(int i = 0; i < num_of_vectors; i++)
            s_values[i] = 0.0;
    }
    else{
        vector_item<T>* old_centroid = centroid;
        s_values = new double[num_of_vectors];
        for(int i = 0; i < num_of_vectors - 1; i++){
            
            /* Get current vector from cluster */
            vector_item<T>* curr_vector = vectors[i];

            /* Calculate distance and add to distances */
            double distance = dist(*old_centroid, *curr_vector);
            s_values[0] += distance;
            s_values[i] = distance;
        }

        s_values[0] /= num_of_vectors;

        int sec_best = get_second_best(*old_centroid, cluster_num, clusters, dist);
        double b_value = calculate_b(*old_centroid, clusters[sec_best], dist);
    
        if(s_values[0] >= b_value){
            s_values[0] = (b_value - s_values[0]) / s_values[0];
        }
        else{
            s_values[0] = (b_value - s_values[0]) / b_value;
        }

        s_of_cluster += s_values[0];
    }


    
    for(int i = 0; i < num_of_vectors; i++){
        /* Get current vector from cluster */
        vector_item<T>* curr_vector = vectors[i];

        /* Calculate the distance with the remaining vectors */
        for(int j = i + 1; j < num_of_vectors - 1; j++){

            vector_item<T>* next_vector = vectors[j];

            /* Calculate distance and add to distances */
            double distance = dist(*curr_vector, *next_vector);
            s_values[i] += distance;
            s_values[j] += distance;  

        }

        s_values[i] /= double(num_of_vectors);

        int sec_best = get_second_best(*curr_vector, cluster_num, clusters, dist);
        double b_value = calculate_b(*curr_vector, clusters[sec_best], dist);
    
        if(s_values[i] >= b_value){
            s_values[i] = (b_value - s_values[i]) / s_values[i];
        }
        else{
            s_values[i] = (b_value - s_values[i]) / b_value;
        }

        s_of_cluster += s_values[i];
    }

    s_of_cluster /= (double)num_of_vectors;

    cout << "s(" << cluster_num << ") = " << s_of_cluster << endl;
}


/* Returns number(id) of cluster */
template <class T>
int cluster<T>::get_cluster_num(){
    return this->cluster_num;
}

/* Returns type of centroid */
template <class T>
int cluster<T>::get_centroid_type(){
    return this->centroid_type;
}

/* Returns pointer to centroid */
template <class T>
vector_item<T>* cluster<T>::get_centroid(){
    return this->centroid;
}

/* Returns pointer to vector in given index of vector */
template <class T>
vector_item<T>* cluster<T>::get_vector(int index){
    return vectors[index ];
}

template <class T>
vector<vector_item<T>*>& cluster<T>::get_vectors(){
    return vectors;
}

/* Returns number of vectors in cluster */
template <class T>
int cluster<T>::get_size(){
    return this->vectors.size();
}

template <class T>
void cluster<T>::set_centroid(vector_item<T>* new_centroid){
    this->centroid = new_centroid;
}

template <class T>
void cluster<T>::set_centroid_type(int type){
    this->centroid_type = type;
}

template <class T>
void cluster<T>::print(){
    cout << "Cluster #" << this->get_cluster_num() << " || Centroid: " << get_centroid()->get_id() << endl;

    for(unsigned int i = 0; i < vectors.size(); i++){
        cout << "\t" << i << ". " << vectors[i]->get_id() << endl;
    }
}



///////////////////
// CL_MANAGEMENT //
///////////////////
template <class T>
cl_management<T>::cl_management(int metric, int k, int init_alg, int assign_alg, int update_alg){
    all_vectors = NULL;
    init_algorithm = NULL;

    this->metric = metric;
    this->k = k;

    for(int i = 0; i < k; i++)
        clusters.push_back(new cluster<T>(i));

    if(init_alg == 1)
        init_algorithm = new cl_init_random<T>;
    else if(init_alg == 2)
        init_algorithm = new cl_init_kmeans<T>;

    if(assign_alg == 1){
        assign_algorithm = new cl_assign_lloyd<T>;
    }

    if(update_alg == 1){
        update_algorithm = new cl_update_kmeans<T>;
    }
    else if(update_alg == 2){
        update_algorithm = new cl_update_pam<T>;
    }

    if(metric == 1)
        dist_function = &eucl_distance<T>;
    else if(metric == 2)
        dist_function = &cs_distance<T>;
}

template <class T>
void cl_management<T>::fill_dataset(ifstream& input){
    string line;
    
    if(this->all_vectors != NULL){
        cout << "Error. Dataset already created" << endl;
        exit(-1);
    }

    if(this->vectors_info.size() != 0){
        cout << "Error. Vectors info already created" << endl;
        exit(-1);
    }


    /* Create new dataset */
    this->all_vectors = new dataset<T>;

    int i = 0;
    /* Scan all vectors in input file, line by line */
    while(getline(input, line)){
        this->all_vectors->add_vector(line);
        this->vectors_info.push_back(new cluster_info(i++));
    }
}

template <class T>
void cl_management<T>::init_clusters(){
    this->init_algorithm->init_clusters(*this);
}

template <class T>
void cl_management<T>::assign_clusters(){
    this->assign_algorithm->assign_clusters(*this);
}

template <class T>
void cl_management<T>::update_clusters(){
    this->update_algorithm->update_clusters(*this);
}



template <class T>
dataset<T>* cl_management<T>::get_dataset(){
    return this->all_vectors;
}

template <class T>
vector<cluster_info*>& cl_management<T>::get_vectors_info(){
    return this->vectors_info;
}

template <class T>
vector<cluster<T>*>& cl_management<T>::get_clusters(){
    return this->clusters;
}


/* Returns number of k(total clusters) */
template <class T>
int cl_management<T>::get_k(){
    return this->k;
}


template <class T>
void cl_management<T>::print(){
    for(unsigned int i = 0; i < clusters.size(); i++)
        clusters[i]->print();
}

template <class T>
dist_func cl_management<T>::get_dist_func(){
    return this->dist_function;
}

template <class T>
void cl_management<T>::evaluation(){
    for(int i = 0; i < k; i++){
        cluster<T>* curr_cluster = clusters[i];
        curr_cluster->evaluation(clusters, dist_function);
    }
}
