/*******************************/
/* cl_algorithms.h */

/* Name:    Andreas Charalampous
 * A.M :    1115201500195
 * e-mail:  sdi1500195@di.uoa.gr
 */
/********************************/
#pragma once

#include "clusters.h"

template<class T> class cl_management;

/* Implementation of all algorithms that will be used in clustering */

/* Initializing algorithms */


/* Abstract type of initializing algorithm */
template <class T>
class cl_init_algorithm{
    protected:
        
    public:
        /* All init algorithms must have this method */
        virtual void init_clusters(cl_management<T>&) = 0;
};

/* Algorithm that selects randomly each centroid of cluster */
template <class T>
class cl_init_random: public cl_init_algorithm<T>{
    private:

    public:
        void init_clusters(cl_management<T>&); 

};

/* Implementation of k-means++ algorithm */
template <class T>
class cl_init_kmeans: public cl_init_algorithm<T>{
    private:
        struct dist_mapping{
            int index;
            double distance;

            dist_mapping(int i, double d) : index(i), distance(d){ }
        };
    public:
        void init_clusters(cl_management<T>&);
        
        /* Given the probabilities array, perform distribution and return */
        /* the array cell that was chosen, after doing binary search      */
        int get_next_centroid(std::vector<dist_mapping>&);
};