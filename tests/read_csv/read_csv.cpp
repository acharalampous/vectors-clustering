#include <iostream>
#include <fstream>
#include <unistd.h>


using namespace std;

int main(void){
    ifstream new_file("twitter_dataset_small_v2.csv");
    string line;
    string id;
    while(getline(new_file, line)){
        cout << line << endl;
        size_t prev = 0, pos;
        int i = 0;

        /* reading id */
        pos = line.find_first_of("\t,", prev);
        
        id = line.substr(prev, pos - prev);
        prev = pos + 1;

        cout << "ID = " << id << endl;
        sleep(1);

        cout << "COORDINATES:" << endl;

        while ((pos = line.find_first_of("\t,", prev)) != string::npos)
        {
            // sleep(1);
            if (pos > prev)
                cout << i << ". " << line.substr(prev, pos - prev) << endl;
            prev = pos+1;
            i++;
        }
        
        if (prev < line.length())
            cout << i << ". " << line.substr(prev, string::npos) << endl;
        sleep(3);
    }

}