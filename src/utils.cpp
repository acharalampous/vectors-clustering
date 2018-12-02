/*******************************/
/* utils.cpp */

/* Name:    Andreas Charalampous
 * A.M :    1115201500195
 * e-mail:  sdi1500195@di.uoa.gr
 */
/********************************/
#include <iostream>
#include <cmath>
#include <string>

#include "utils.h"
#include "metrics.h"

using namespace std;

template float vector_product(std::array<float,D>&, std::array<int ,D>&);
template float vector_product(std::array<float,D>&, std::array<double, D>&);
template double eucl_distance(vector_item<double>&, vector_item<double>&);
template double cs_distance(vector_item<double>&, vector_item<double>&);
template int get_second_best(vector_item<double>&, int, vector<cluster<double>*>&, dist_func&);
template double calculate_b(vector_item<double>&, cluster<double>*, dist_func&);
template float exchausting_s(dataset<int>&, vector_item<int>&, int);


/*  All functions implementions that are defined in utils.h */

int get_parameters(int argc, char** argv, string& input_file, string& queryset_file, string& output_file, int& k, int& L){
    for (int i = 1; i < argc; i += 2){ // get all parameters
        char par; // parameter given
            
        par = argv[i][1]; // get parameter
        switch(par){
            case 'd':{ // input file provided
                if(!input_file.empty()){ // input file(dataset) is provided twice
                    printf("Error in parameters! Input file [-d] is given more than once! Abort.\n");
                    return -2;
                }

                input_file = argv[i + 1];
                break;
            } // end case -d
            case 'q':{ // query file parameter
                if(!queryset_file.empty()){ // port given twice
                    printf("Error in parameters! Query file [-q] is given more than once! Abort.\n");
                    return -2;
                }
                
                queryset_file = argv[i + 1];
                break;
            } // end case -q
            case 'o':{ // output file parameter
                if(!output_file.empty()){ // port given twice
                    printf("Error in parameters! Output file [-o] is given more than once! Abort.\n");
                    return -2;
                }

                output_file = argv[i + 1];
                break;
            } // end case -o
            case 'k':{ // number of hash functions parameter
                if(k != -1){ // num of hash functions given twice
                    printf("Error in parameters! Number of hash functions [-k] is given more than once! Abort.\n");
                    return -2;
                }

                if(!isNumber(argv[i + 1])){ // num of hash functions is not a number
                    printf("Error in parameters! Number of hash functions [-k] given is not a number. Abort\n");
                    return -2;
                }

                k = atoi(argv[i + 1]);
                if(k <= 0){ // not positive non-zero number of threads given
                    printf("Error in parameters! Number of hash functions [-k] is not a positive non-zero number. Abort\n");
                    return -2;
                }

                break;
            } // end case -k
            case 'L':{ // number of hash tables to be created parameter
                if(L != -1){ // num of hash tables given twice
                    printf("Error in parameters! Number of hash tables [-L] is given more than once! Abort.\n");
                    return -2;
                }

                if(!isNumber(argv[i + 1])){ // num of hash tables is not a number
                    printf("Error in parameters! Number of hash tables [-L] given is not a number. Abort\n");
                    return -2;
                }

                L = atoi(argv[i + 1]);
                if(L <= 0){ // not positive non-zero number of threads given
                    printf("Error in parameters! Number of hash functions [-L] is not a positive non-zero number. Abort\n");
                    return -2;
                }

                break;
            } // end case -d
            default:{
                printf("Error in parameters. Unknown parameter given [%s]. Abort.\n", argv[i]);
                return -5;
            }
        } // end switch
    } // end for
    return 0;
}


int HC_get_parameters(int argc, char** argv, string& input_file, string& queryset_file, string& output_file, int& k, int& probes, int& M){
    for (int i = 1; i < argc; i += 2){ // get all parameters
        char par; // parameter given
            
        par = argv[i][1]; // get parameter
        switch(par){
            case 'd':{ // input file provided
                if(!input_file.empty()){ // input file(dataset) is provided twice
                    printf("Error in parameters! Input file [-d] is given more than once! Abort.\n");
                    return -2;
                }

                input_file = argv[i + 1];
                break;
            } // end case -d
            case 'q':{ // query file parameter
                if(!queryset_file.empty()){ // port given twice
                    printf("Error in parameters! Query file [-q] is given more than once! Abort.\n");
                    return -2;
                }
                
                queryset_file = argv[i + 1];
                break;
            } // end case -q
            case 'o':{ // output file parameter
                if(!output_file.empty()){ // port given twice
                    printf("Error in parameters! Output file [-o] is given more than once! Abort.\n");
                    return -2;
                }

                output_file = argv[i + 1];
                break;
            } // end case -o
            case 'k':{ // number of hash functions parameter
                if(k != -1){ // num of hash functions given twice
                    printf("Error in parameters! Number of hash functions [-k] is given more than once! Abort.\n");
                    return -2;
                }

                if(!isNumber(argv[i + 1])){ // num of hash functions is not a number
                    printf("Error in parameters! Number of hash functions [-k] given is not a number. Abort\n");
                    return -2;
                }

                k = atoi(argv[i + 1]);
                if(k <= 0){ // not positive non-zero number of threads given
                    printf("Error in parameters! Number of hash functions [-k] is not a positive non-zero number. Abort\n");
                    return -2;
                }

                break;
            } // end case -k
            case 'p':{ // number of probes parameter
                string par_str(argv[i]);
                if(par_str.compare("-probes")){
                    printf("Error in parameters. Unknown parameter given [%s]. Abort.\n", argv[i]);
                    return -5;
                }

                if(probes != -1){ // num of probes given twice
                    printf("Error in parameters! Number of probes [-probes] is given more than once! Abort.\n");
                    return -2;
                }

                if(!isNumber(argv[i + 1])){ // num of probes is not a number
                    printf("Error in parameters! Number of probes [-probes] given is not a number. Abort\n");
                    return -2;
                }

                probes = atoi(argv[i + 1]);
                if(probes <= 0){ // not positive non-zero number of probes given
                    printf("Error in parameters! Number of probes [-probes] is not a positive non-zero number. Abort\n");
                    return -2;
                }

                break;
            } // end case -d
            case 'M':{ // number of total points to be checked parameter
                if(M != -1){ // M given twice
                    printf("Error in parameters! Number of total points to be checked [-M] is given more than once! Abort.\n");
                    return -2;
                }

                if(!isNumber(argv[i + 1])){ // num of probes is not a number
                    printf("Error in parameters! Number of total points to be checked [-M] given is not a number. Abort\n");
                    return -2;
                }

                M = atoi(argv[i + 1]);
                if(M <= 0){ // not positive non-zero number of M given
                    printf("Error in parameters! Number of total points to be checked [-M] is not a positive non-zero number. Abort\n");
                    return -2;
                }

                break;
            } // end case -M
            default:{
                printf("Error in parameters. Unknown parameter given [%s]. Abort.\n", argv[i]);
                return -5;
            }
        } // end switch
    } // end for
    return 0;
}


void printValidParameters(){
    cout << "\n*Execute again providing (optionally) the following parameters:" << endl;
    cout << "\t-d inputFile" << endl;
    cout << "\t-q queryFile" << endl;
    cout << "\t-o outputFile" << endl;
    cout << "\t-k K" << endl;
    cout << "\t-L L" << endl;
    cout << "-inputFile: Path to the dataset file" << endl;
    cout << "-queryFule: Path to the query file" << endl;
    cout << "-outputFile: Path to the output file" << endl;
    cout << "-K: Number of hash functions for each hash table, [>=1]" << endl;
    cout << "-L: Number of hash tables, [>=1]" << endl;
}

void HC_printValidParameters(){
    cout << "\n*Execute again providing (optionally) the following parameters:" << endl;
    cout << "\t-d inputFile" << endl;
    cout << "\t-q queryFile" << endl;
    cout << "\t-o outputFile" << endl;
    cout << "\t-k K" << endl;
    cout << "\t-probes probes" << endl;
    cout << "\t-M M" << endl;
    cout << "-inputFile: Path to the dataset file" << endl;
    cout << "-queryFule: Path to the query file" << endl;
    cout << "-outputFile: Path to the output file" << endl;
    cout << "-K: Number of hash functions for each hash table, [>=1]" << endl;
    cout << "-probes: Number of neighbouring bucket to be checked, [>=1]" << endl;
    cout << "-M: Number of total points to be checked for a query, [>=1]" << endl;
}


int isNumber(char* str){
    char* temp = str;
    char ch = *temp;
    if(ch == '-') // negative number is provided
        ch = *(++temp);
    while(ch != '\0'){
        if(ch >= 48 && ch <= 57){
            ch = *(++temp);
        }
        else{
            return 0;
        }
    }
    return 1;
}

int get_int_len(int num){
    int num_digits = 1; // digits of number
    while(num >= 0){	
		num_digits++;
       	num = num >> 1;
	}

    return num_digits;
}

int get_combinations(int a, int b){
    /* Computes comb(a,b) = a! / (b! * (a - b)!) */
    int diff = a - b;

    return ( get_factorial(a) / ( get_factorial(b) * get_factorial(diff)) );


}

int get_factorial(int num){
    int total = 1;
    for(int i = 2; i <= num; i++){
        total *= i;
    }

    return total;
}

int hamming_dist(int x, int y){
    int x_or = x ^ y; // keep different bits, they will be represented as 1

    int dist = 0; // hamming_distance
    while(x_or > 0){ // check all bits
        if((x_or | 1) == x_or) // if true, means the leftmost bit is 1, hence different
            dist++;
        x_or = x_or >> 1;
    }
    return dist;
}

/* Compute inner product of two vectors */
template <class T>
float vector_product(std::array<float, D>& vec1, std::array<T, D>& vec2){
    float product = 0.0;

    for(unsigned int i = 0; i < D; i++)
        product += vec1[i] * vec2[i];

    return product;
}

long long int h_concantenate(vector<int>& hs){
    /* Creates a string with the concatenated ints, and then it is */
    /* returned converted back to int.                             */

    string res_str("");
    for(unsigned int i = 0; i < hs.size(); i++)
        res_str += to_string(abs(hs[i]));


    return stoll(res_str);
}

/* Computes and returns the modulo of the given values a mod b */
long long int my_mod(int a, long int b){
	long long int res = a - ((floor( (long double)a / (long double)b) ) * b); 
    return res;
}

int get_metrics(string& metrics){
    if(metrics.empty()) // no definition, use euclidean
        return 0;

    if(metrics.find("euclidean") != string::npos) // euclidean is defined
        return 1;
    
    if(metrics.find("cosine") != string::npos) // cosine is defined
        return 2;

    return 0; // nothing was defined, use euclidean
}

double get_radius(std::string& radius){
    if(radius.compare(0, 9, "Radius: <")) // check if radius is given
        return -1;

    /* Extract radius */
    int pos1 = radius.find('<');
    int pos2 = radius.find('>');

    /* Keep the part enclosed in <> */
    string temp = radius.substr(pos1 + 1, pos2 - pos1 - 1);
    
    double result = 0.0;
    try{
        result = stod(temp);
    }catch(invalid_argument){
        cout << "Invalid radius was given. Radius won't be used at all!" << endl;
        return 0.0;
    }

    return result;
}

template <class T>
double eucl_distance(vector_item<T>& vec1, vector_item<T>& vec2){
	/* Compute sqrt((Sum((ai - bi)^2))), for all points */

    if(vec1.get_size() != vec2.get_size()){
		cout << "Invalid dimensions of vectors" << endl;
		exit(0);
	}

	double dist = 0.0;

	array<T, D>& arr1 = vec1.get_points();
	array<T, D>& arr2 = vec2.get_points();


	for(unsigned int i = 0; i < arr1.size(); i++)
		dist += pow(arr1[i] - arr2[i], 2);

	dist = sqrt(dist);

	return dist;
}

template <class T>
double cs_distance(vector_item<T>& vec1, vector_item<T>& vec2){
	/* Compute (1 - ((vec1 * vec2) / (||vec1||*||vec2||))) */
 	if(vec1.get_size() != vec2.get_size()){
		cout << "Invalid dimensions" << endl;
		exit(0);
	}

	double dist;
	double euc_dist = 0.0;
	double arr1_norm = 0.0;
	double arr2_norm = 0.0;

	array<T, D>& arr1 = vec1.get_points();
	array<T, D>& arr2 = vec2.get_points();


	for(unsigned int i = 0; i < arr1.size(); i++){
		euc_dist += (arr1[i] * arr2[i]); // vec1 & vec2
		arr1_norm += arr1[i] * arr1[i]; // ||vec1||
		arr2_norm += arr2[i] * arr2[i]; // ||vec2||
	}

	arr1_norm = sqrt(arr1_norm);
	arr2_norm = sqrt(arr2_norm);

	dist = euc_dist / (arr1_norm * arr2_norm);
	dist = 1 - dist;

	return dist;
}

template <class T>
int get_second_best(vector_item<T>& query, int cl_num, vector<cluster<T>*>& clusters, dist_func& dist){
    int sec_best = cl_num; // index of the second best cluster for query to be placed
    double min_distance; // distance to second best
    int k = clusters.size();

    int flag = 0; // min_distance was initialized

    /* Check for all clusters */
    for(int i = 0; i < k; i++){
        if(i == cl_num) // if same cluster, skip
            continue;
        
        double curr_dist = dist(query, *(clusters[i]->get_centroid()));
        if(flag == 0){ // no minimum distance yet, so keep current distance
            sec_best = i;
            min_distance = curr_dist;
            flag = 1;
        }
        else{
            if(curr_dist <= min_distance){ // keep if less
                sec_best = i;
                min_distance = curr_dist;
            }
        }
    } // end for all clusters

    return sec_best;
}

template <class T>
double calculate_b(vector_item<T>& query, cluster<T>* cl, dist_func& dist){
    int num_of_vectors = cl->get_size();
    double b_value = 0.0;

    /* Calculate sum of distances of all vectors in cluster given */
    for(int i = 0; i < num_of_vectors; i++)
        b_value += dist(query, *(cl->get_vector(i)));

    if(cl->get_centroid_type() == 1){
        num_of_vectors++;

        b_value += dist(query, *(cl->get_centroid()));
    }

    if(num_of_vectors == 0)
        return 0.0;

    /* Get average distance and return */
    b_value = b_value / (double)num_of_vectors;
    
    return b_value;
}

double get_starting_r(vector<cluster<double>*>& clusters, dist_func& dist){
    int num_of_clusters = clusters.size(); // get number of clusters
    double min_distance = 0.0;
    vector<vector_item<double>*> centroids;

    /* Get all clusters from the begining, to save additional fetchs */ 
    for(int i = 0; i < num_of_clusters; i++)
        centroids.push_back(clusters[i]->get_centroid());

    /* Initialize min distance */
    min_distance = dist(*centroids[0], *centroids[1]);
    for(int i = 2; i < num_of_clusters; i++){
        double curr_dist = dist(*centroids[0], *centroids[i]);
        if(curr_dist <= min_distance && curr_dist != 0.0){
            min_distance = curr_dist;
        }
    }

    /* Check for all clusters the distance between centroids */ 
    for(int i = 1; i < num_of_clusters; i++){
        for(int j = i + 1; j < num_of_clusters; j++){
            double curr_dist = dist(*centroids[i], *centroids[j]);
            if(curr_dist <= min_distance && curr_dist != 0.0){
                min_distance = curr_dist;
            }
        }
    }


    return min_distance / 2.0;
}

int add_to_clusters(cluster<double>* cl, double& r, vector<vector_check*>& vectors_to_check, vector<cluster_info*>& vectors_info, int& vectors_left){
	int cluster_num = cl->get_cluster_num();
	int changes = 0;

	int num_of_vectors = vectors_to_check.size();

	/* Check all vectors that were collected before */
	for(int i = 0; i < num_of_vectors; i++){
		vector_item<double>& item = *(vectors_to_check[i]->item);

		/* Check if item was not checked already */
		double dist = vectors_to_check[i]->distance;
		if(dist <= r){
			int item_index = item.get_index();

			/* Vector is in radius, must check if its already assigned */
			if(vectors_info[item_index]->get_cluster_num() == -1){ // not assigned
				vectors_info[item_index]->set_cluster(cluster_num);
				vectors_info[item_index]->set_distance(dist);
				vectors_left--;
				changes++; // changes were made
			}
			else if(dist <= vectors_info[item_index]->get_distance()){ // vector is assigned, check if less distance
					vectors_info[item_index]->set_cluster(cluster_num);
					vectors_info[item_index]->set_distance(dist);
			}

			delete vectors_to_check[i];
			
			num_of_vectors--;
			i--;
			vectors_to_check.erase(vectors_to_check.begin() + i + 1); 	
		}
	}
	return changes;
}

void final_assign(cl_management<double>& cl_manage){
    dataset<double>* all_vectors = cl_manage.get_dataset();
    vector<cluster_info*>& vectors_info = cl_manage.get_vectors_info();
    vector<cluster<double>*>& clusters = cl_manage.get_clusters();
    
    dist_func dist_function = cl_manage.get_dist_func();

    int k = cl_manage.get_k(); // get number of total clusters
    int num_of_vectors = all_vectors->get_counter(); // get num of vectors

    /* Assign all non-centroids in dataset, to nearest centroid's cluster */ 
    for(int i = 0; i < num_of_vectors; i++){
        /* Check if not centroid */
        if(vectors_info[i]->get_is_centroid() == 1){
            continue;
        }

        vector_item<double>* curr_item = all_vectors->get_item(i);

        int cluster_num = vectors_info[i]->get_cluster_num();
        /* Check if cluster for current vector was found */
        if(cluster_num != -1){
            clusters[cluster_num]->add_vector(curr_item);
            continue;
        }

        /* Get current vector */
        double min_distance = 0.0;
        
        /* Check for all centroids and find the nearest */
        for(int j = 0; j < k; j++){
            vector_item<double>* curr_centroid = clusters[j]->get_centroid();

            /* Get distance and check if minimum */
            double dist = dist_function(*curr_item, *curr_centroid);
            if(j == 0 || dist <= min_distance){
                min_distance = dist;
                cluster_num = j;
            }
        }

        clusters[cluster_num]->add_vector(curr_item);
        
        vectors_info[i]->set_cluster(cluster_num);
        vectors_info[i]->set_distance(min_distance);
    }
}



int new_execution(ifstream& input, ifstream& query, ofstream& output){
    string choice;

    while(1){
        cout << "Would you like to continue(y or n): ";
        fflush(stdout);
        
        getline(cin, choice);
        fflush(stdin);

        if(!choice.compare("y") || !choice.compare("Y") || !choice.compare("yes") ){
            cout << "Please choose one of the following(select number):" << endl;
            cout << "  0 Same files as before" << endl;
            cout << "  1 New Dataset" << endl;
            cout << "  2 New Query File" << endl;
            cout << "  3 New Output File" << endl;
            cout << "  4 New Dataset and Query File" << endl;
            cout << "  5 New Dataset and Output File" << endl;
            cout << "  6 New Query and Output File" << endl;
            cout << "  7 Entirely New Files" << endl;
            cout << "  8 I Changed My Mind. I Want To Exit\n" << endl;
            cout << "Note: If you use the same output file name, the previous one will be replaced!\n" << endl;
            cout << "My choice: ";

            fflush(stdout);
            getline(cin, choice);
            fflush(stdin);

            if(!choice.compare("0")) // Use same files
                return 0;
            else if(!choice.compare("1")){ // diff input{
                input.close();
                return 1;
            }
            else if(!choice.compare("2")){ // diff query
                query.close();
                return 2;
            }
            else if(!choice.compare("3")){ // diff out
                output.close();
                return 3;
            }
            else if(!choice.compare("4")){ // diff in+q
                input.close();
                query.close();
                return 4;
            }
            else if(!choice.compare("5")){ // diff in+out
                input.close();
                output.close();
                return 5;
            }
            else if(!choice.compare("6")){ // diff q+out
                query.close();
                output.close();
                return 6;
            }
            else if(!choice.compare("7")){ // diff ALL
                input.close();
                query.close();
                output.close();
                return 7;
            }
            else if(!choice.compare("8")) // EXIT
                return 8;
            else
                cout << "Invalid input given. Try again." << endl;
        }
        if(!choice.compare("n") || !choice.compare("N") || !choice.compare("no") )
            return 8; // EXIT
        else
            cout << "Invalid input given. Try again." << endl;
    }
}

int in_set(unordered_set<string>& check_set, string& id){
    unordered_set<string>::iterator it;

    it = check_set.find(id);
    if(it == check_set.end()) // id does not exist
        return 0;
    else // id is in set
        return 1; 
}

template <class T>
float exchausting_s(dataset<T>& data_set, vector_item<T>& item, int metric){
    float min_dist = -1.0;
    float cur_dist = 0.0;

    int data_sz = data_set.get_counter();
    if(metric == 1){ // euclidean distance

        /* Check all items in dataset */
        for(int i = 0; i < data_sz; i++){
            vector_item<T>* current_item = data_set.get_item(i);
            cur_dist = eucl_distance(*current_item, item);

            /* keep distance, if smaller than the minimum */
            if(cur_dist <= min_dist || i == 0)
                min_dist = cur_dist;
        }
    }
    else if(metric == 2){ // cosine distance

        /* Check all items in dataset */
        for(int i = 0; i < data_sz; i++){
            vector_item<T>* current_item = data_set.get_item(i);
            cur_dist = cs_distance(*current_item, item);

            /* keep distance, if smaller than the minimum */
            if(cur_dist <= min_dist || i == 0)
                min_dist = cur_dist;
        }
    }
    return min_dist;
}