/*******************************/
/* clusters.h */

/* Name:    Andreas Charalampous
 * A.M :    1115201500195
 * e-mail:  sdi1500195@di.uoa.gr
 */
/********************************/
#pragma once
#include <fstream>
#include <vector>

#include "dataset.h"
#include "cl_algorithms.h"

template <class T> class cl_init_algorithm;

typedef double (*dist_func)(vector_item<double>&, vector_item<double>&);

/* Implementation of all necessary structs and methods that 
 * will be used in order to implement all the clustering 
 * algorithms 
 */


/* Class that holds info about the cluster that a vector_item belongs to */
class cluster_info{
    private:
        int is_centroid; // 1 if vector is the centroid of the cluster that it belongs to, else 0
        int cluster_num; // number of cluster that vector belongs to
        int index; // index of vector in dataset(that is mapped to this info)
    public:
        /* Con-De Structor */

        /* Given the index of vector, create a new mapped info about it */
        cluster_info(int);

        /* Accessors */
        int get_is_centroid(); 
        int get_cluster_num();
        int get_index();

        /* Mutators */
        void set_centroid(); // sets is_centroid as 1
        void set_cluster(int); // sets cluster_num to the given int 
        void set_index(int); // sets index to the given int
        void reset_info(); // sets cluster_num as -1
};


/* Class that represents the cluster. Holds all necessary info of cluster, */
/* like pointer to centroid and all points in it                           */
template <class T>
class cluster{
    private:
        int cluster_num; // number of cluster
        int centroid_type; // 1: if centroid in dataset, 0: not in dataset, -1 no centroid
        vector_item<T>* centroid; // pointer to vector in dataset that is centroid
        std::vector<vector_item<T>*> vectors; // holds pointer to all vectors that are in cluster


    public:
        /* Con-De Structor */

        /* Given the number of cluster, create a new empty cluster */
        cluster(int);

        /* Accessors */
        int get_cluster_num();
        int get_centroid_type();
        vector_item<T>* get_centroid(); 
        int get_size();

        /* Mutators */
        void set_centroid(vector_item<T>*); // sets is_centroid as 1
        void set_centroid_type(int); // sets centroid type 
        // void set_index(int); // sets index to the given int
        // void reset_info(); // sets cluster_num as -1
        
        void print();
};

/* Class that has all necessary structs collected */
template <class T>
class cl_management{
    private:
        dataset<T>* all_vectors; // dataset with all vectors
        std::vector<cluster_info*> vectors_info; // cluster info for each vector
        std::vector<cluster<T>*> clusters; // all clusters
        cl_init_algorithm<T>* init_algorithm; // algorithm for initiliazing clusters
        dist_func dist_function; // pointer to distance function, eucliean or cosine

        int metric; // metric to be used, 1: euclidean, 2: cosine
        int k; // number of total clusters

    public:
        /* Con-De Structor */

        /* Given the metric, algorithms and number of clusters, initiliaze clusters controlling class */
        cl_management(int, int, int);

        /* Given a file stream, get all vectors and assign in dataset */
        void fill_dataset(std::ifstream&);
        void init_clusters();

        /*Accessors */
        dataset<T>* get_dataset();
        std::vector<cluster_info*>& get_vectors_info();
        std::vector<cluster<T>*>& get_clusters();
        int get_k(); // returns number of k
        dist_func get_dist_func();

        void print();

};