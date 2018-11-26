/*******************************/
/* clusters.cpp */

/* Name:    Andreas Charalampous
 * A.M :    1115201500195
 * e-mail:  sdi1500195@di.uoa.gr
 */
/********************************/
#include <iostream>

#include "clusters.h"

using namespace std;

template class cluster<double>;


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

/* Unassign vector from cluster, setting cluster_num to -1 */
void cluster_info::reset_info(){
    this->cluster_num = -1;
}



/////////////
// CLUSTER //
/////////////
template <class T>
cluster<T>::cluster(int cluster_num){
    this->cluster_num = cluster_num;
    this->real_centroid = -1;
    this->centroid = NULL;
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

/* Returns pointer to centroid */
template <class T>
vector_item<T>* cluster<T>::get_centroid(){
    return this->centroid;
}

/* Returns number of vector in cluster */
template <class T>
int cluster<T>::get_size(){
    return this->vectors.get_size();
}