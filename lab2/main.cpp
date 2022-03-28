#include <iostream>
#include <vector>
#include <cmath>
#include "FORSYTHE.H"

double getNorm(const std::vector<double> &x) {
    double sum = 0.0;
    for (double i : x) {
        sum += std::pow(i, 2);
    }
    return std::sqrt(sum);
}

int main() {
    const std::vector<double> parameters = {1, 0.1, 0.01, 0.0001, 0.000001};

    for (double p : parameters) {

        const std::vector<std::vector<double>> matrix = {
                {p - 3, -4,  -4,  7,   2,  3,  8,  7},
                {0,     -15, -1,  5,   -3, 6,  6,  -6},
                {-4,    2,   -16, 7,   0,  8,  -7, 6},
                {0,     8,   -5,  -11, 1,  0,  4,  5},
                {8,     6,   -8,  4,   27, -7, -1, 5},
                {-4,    -2,  1,   2,   -8, 10,  7, 0},
                {0,     -1,  5,   2,   -8, 2,  -2, 0},
                {0,     -8,  -7,  3,   -7, -4, -8, 5}
        };

        std::vector<double> B = {
                {(2 * p + 54),
                 -72,
                 -33,
                 -15,
                 180,
                 -5,
                 -14,
                 -131}
        };

        std::vector<std::vector<double>> transMatrix (matrix.size(), std::vector<double>(matrix.size()));
        for (int i = 0; i < 8; i++) {
            for (int j = 0; j < 8; j++) {
                transMatrix[i][j] = matrix[j][i];
            }
        }

        std::vector<std::vector<double>> newMatrix (matrix.size(), std::vector<double>(matrix.size()));
        for (int i = 0; i < 8; i++) {
            for (int j = 0; j < 8; j++) {
                newMatrix[i][j] = 0;
                for (int z = 0; z < 8; z++) {
                    newMatrix[i][j] += transMatrix[i][z] * matrix[z][j];
                }
            }
        }

        std::vector<double> newB (B.size());
        for (int i = 0; i < 8; i++) {
            newB[i] = 0;
            for (int z = 0; z < 8; z++) {
                newB[i] += transMatrix[i][z] * B[z];
            }
        }

        DECOMP decompFirst(matrix.size(), matrix);
        decompFirst.Solve(B);

        DECOMP decompSecond(newMatrix.size(), newMatrix);
        decompSecond.Solve(newB);

        const double norm = getNorm(B);

        std::cout << "\nCOND (first matrix) - " << decompFirst.Cond() << '\n';
        std::cout << "COND (second matrix) - " << decompSecond.Cond() << '\n';
        std::cout << "σ - " << abs((norm - getNorm(newB)) / norm) << '\n';
    }
}
