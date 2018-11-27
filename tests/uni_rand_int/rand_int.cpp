#include <iostream>
#include <random>

using namespace std;

int main(void){
    random_device rd;
    mt19937 gen(rd()); // choose random vector as centroids

    int min = 0;
    int max = 2;

    uniform_real_distribution<> rand_Z(0.0, 0.03123);
	cout.precision(15);
    for(int i = 0; i < 10; i++)
        cout << rand_Z(gen) << endl;
}
