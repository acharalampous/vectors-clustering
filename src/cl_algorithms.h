/*******************************/
/* cl_algorithms.h */

/* Name:    Andreas Charalampous
 * A.M :    1115201500195
 * e-mail:  sdi1500195@di.uoa.gr
 */
/********************************/
#pragma once

/* Implementation of all algorithms that will be used in clustering */

/* Initializing algorithms */


/* Abstract type of initializing algorithm */
class cl_init_algorithm{
    private:

    public:
        /* All init algorithms must have this method */
        virtual void init_clusters() = 0;
};

/* Algorithm that selects randomly each centroid of cluster */
class cl_init_random: public cl_init_algorithm{
    private:

    public:
        void init_clusters(); 
        
};