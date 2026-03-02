#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>

// Διαχείριση Συμμετρικού Πίνακα
class SymmetricMatrix {
    int N;
    std::vector<double> data;
public:
    SymmetricMatrix(int n) : N(n) { data.assign(n * (n + 1) / 2, 0.0); }
    
    int idx(int i, int j) const {
        if (i > j) std::swap(i, j);
        return i * N - (i * (i - 1) / 2) + (j - i); // [cite: 4]
    }
    
    void add(int i, int j, double v) { data[idx(i, j)] += v; }
    double get(int i, int j) const { return data[idx(i, j)]; }
    void set(int i, int j, double v) { data[idx(i, j)] = v; }
};

// Επίλυση Συστήματος (Cholesky)
class Solver {
public:
    static void solve(SymmetricMatrix& K, const std::vector<double>& F, std::vector<double>& u, int n) {
        // Παραγοντοποίηση Cholesky [cite: 17, 19]
        for (int i = 0; i <= n; i++) {
            for (int j = 0; j <= i; j++) {
                double s = 0;
                for (int k = 0; k < j; k++) s += K.get(i, k) * K.get(j, k);
                if (i == j) K.set(i, i, std::sqrt(K.get(i, i) - s));
                else K.set(i, j, (1.0 / K.get(j, j)) * (K.get(i, j) - s));
            }
        }
        
        // Forward Substitution [cite: 22, 24]
        std::vector<double> y(n + 1);
        for (int i = 0; i <= n; i++) {
            double s = 0;
            for (int k = 0; k < i; k++) s += K.get(i, k) * y[k];
            y[i] = (F[i] - s) / K.get(i, i);
        }
        
        // Backward Substitution [cite: 25, 26]
        for (int i = n; i >= 0; i--) {
            double s = 0;
            for (int k = i + 1; k <= n; k++) s += K.get(k, i) * u[k];
            u[i] = (y[i] - s) / K.get(i, i);
        }
    }
};

class FEMSolver {
    int n_el;
    double L, h, f_src;
public:
    FEMSolver(int elements, double length, double source) 
        : n_el(elements), L(length), f_src(source), h(length/elements) {}

    void execute() {
        SymmetricMatrix K(n_el + 1);
        std::vector<double> F(n_el + 1, 0.0), u(n_el + 1, 0.0);

        // Συναρμολόγηση (Assembly) [cite: 10, 15]
        double xi[2] = {-0.57735, 0.57735}, w[2] = {1.0, 1.0};
        double J = h / 2.0;

        for (int e = 0; e < n_el; ++e) {
            for (int k = 0; k < 2; ++k) {
                for (int i = 0; i < 2; ++i) {
                    double Ni = (i == 0) ? 0.5*(1-xi[k]) : 0.5*(1+xi[k]);
                    F[e + i] += f_src * Ni * J * w[k];
                    for (int j = 0; j < 2; ++j) {
                        double dNi = (i == 0) ? -0.5 : 0.5;
                        double dNj = (j == 0) ? -0.5 : 0.5;
                        K.add(e + i, e + j, (dNi/J)*(dNj/J)*J*w[k]);
                    }
                }
            }
        }

        // Συνοριακές Συνθήκες (Penalty Method) 
        K.set(0, 0, 1e15); F[0] = 0;
        K.set(n_el, n_el, 1e15); F[n_el] = 0;

        Solver::solve(K, F, u, n_el);
        saveResults(u);
    }

    void saveResults(const std::vector<double>& u) {
        std::ofstream file("results.csv");
        file << "x,u\n";
        for (int i = 0; i <= n_el; ++i) 
            file << i * h << "," << u[i] << "\n";
        file.close();
        std::cout << "Τα αποτελέσματα αποθηκεύτηκαν στο results.csv" << std::endl;
    }
};

int main() {
    FEMSolver fem(10, 2.0, 10.0);
    fem.execute();
    return 0;
}