#include <iostream>
#include <fstream>
#include <string>

using namespace std;

int main(void){
    ifstream conf("cl.conf");
    string line;
    while(getline(conf, line)){
        cout << line << endl;
        if(line.compare(0, 19, "number_of_clusters:") == 0){
            string par = line.substr(19, line.length() - 19);
            cout << "**Got number of clusters: " << par << endl;
        }
        else if(line.compare(0, 25, "number_of_hash_functions:") == 0){
            string par = line.substr(25, line.length() - 25);
            cout << "**Got number of hash functions: " << par << endl;
        }
        else if(line.compare(0, 22, "number_of_hash_tables:") == 0){
            string par = line.substr(22, line.length() - 22);
            cout << "**Got number of hash tables: " << par << endl;
        }
        else if(line.compare(0, 12, "max_updates:") == 0){
            string par = line.substr(12, line.length() - 12);
            cout << "**Got number of updates: " << par << endl;
        }
        else if(line.compare(0, 10, "hc_probes:") == 0){
            string par = line.substr(10, line.length() - 10);
            cout << "**Got number of probes(hc): " << par << endl;
        }
        else if(line.compare(0, 5, "hc_M:") == 0){
            string par = line.substr(5, line.length() - 5);
            cout << "**Got number of clusters: " << par << endl;
        }
    }
}