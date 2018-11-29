/*******************************/
/* cl_algorithms.h */

/* Name:    Andreas Charalampous
 * A.M :    1115201500195
 * e-mail:  sdi1500195@di.uoa.gr
 */
/********************************/
#pragma once

#include "dataset.h"
#include "clusters.h"

template<class T> class cl_management;
template<class T> class cluster;

typedef double (*dist_func)(vector_item<double>&, vector_item<double>&);

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



/* Assignment algorithms */
template <class T>
class cl_assign_algorithm{
    protected:

    public:
        virtual void assign_clusters(cl_management<T>&) = 0;

};

template <class T>
class cl_assign_lloyd : public cl_assign_algorithm<T>{
    private:

    public:
        void assign_clusters(cl_management<T>&);
};


/* Update algorithms */
template <class T>
class cl_update_algorithm{
    protected:

    public:
        virtual void update_clusters(cl_management<T>&) = 0;
};

template <class T>
class cl_update_kmeans : public cl_update_algorithm<T>{
    private:

    public:
        void update_clusters(cl_management<T>&);
        
        vector_item<T>* get_new_centroid(cluster<T>&);
};

template <class T>
class cl_update_pam : public cl_update_algorithm<T>{
    private:

    public:
        void update_clusters(cl_management<T>&);

        /* Calculates the average distance for every vector in given cluster */
        /* and returns the vector with minimum average distance              */
        vector_item<T>* get_new_centroid(cluster<T>&, dist_func&);
};