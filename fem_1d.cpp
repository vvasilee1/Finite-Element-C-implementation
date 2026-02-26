#include <iostream>
#include <vector>

class SymmetricMatrix{
    private:
        int N;
        std::vector<double> data;

    public:
        SymmetricMatrix(int n): N(n){
            data.assign(n*(n+1)/2,0.0);
        }

        int getIndex(int i,int j) const{
            if(i>j) std::swap(i,j);
            return j + (i *(2*N -i-1) / 2);
        }

        void addValue(int i,int j,double val) {
            data[getIndex(i,j)] += val;
        }

        double getValue(int i,int j){
            return data[getIndex(i,j)];
        }
        void printFull() {
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                printf("%6.2f ", getValue(i, j));
            }
            std::cout << "\n";
        }
    }
};

void AssembleSymmetric(SymmetricMatrix& K,int num_elements){
    for(int e=0; e < num_elements; ++e){
        double local_K[2][2] ={{1.0,-1.0},{-1.0,1.0}};
        int n1 = e;
        int n2 = e+1;

        K.addValue(n1,n1,local_K[0][0]);
        K.addValue(n1,n2,local_K[0][1]);
        K.addValue(n2,n1,local_K[1][0]);
        K.addValue(n2,n2,local_K[1][1]);

    }
}

int main(){
    int nodes = 5;
    SymmetricMatrix K(nodes);
    AssembleSymmetric(K,4);
    K.printFull();
    return 0;
}