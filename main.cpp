#include "vecmath.h"
#include <stdio.h>
#include <direct.h>

template <size_t WIDTH, size_t HEIGHT>
using Mtx = BaseMtx<float, WIDTH, HEIGHT>;

template <size_t SIZE>
using Vec = BaseVec<float, SIZE>;

int main(int argc, char** argv)
{
    _mkdir("out");

    Mtx<3, 5> data =
    {
        90, 60, 90,
        90, 90, 30,
        60, 60, 60,
        60, 60, 90,
        30, 30, 30
    };

    // Calculate the mean of the data
    Vec<Columns(data)> mean;
    for (int i = 0; i < Columns(data); ++i)
    {
        mean[i] = 0.0f;
        for (int j = 0; j < Rows(data); ++j)
            mean[i] = Lerp(mean[i], data[j][i], 1.0f / float(j + 1));
    }

    // Make the covariance matrix
    Mtx<Columns(data), Columns(data)> covariance;
    for (int i = 0; i < Columns(covariance); ++i)
    {
        for (int j = 0; j < Rows(covariance); ++j)
        {
            covariance[i][j] = 0.0f;
            for (int k = 0; k < Rows(data); ++k)
                covariance[i][j] += (data[k][i] - mean[i]) * (data[k][j] - mean[j]) / float(Rows(data));
        }
    }

    // Get the eigenvalues
    Vec<3> eigenValues = QRAlgorithm(covariance, 10000, 0.0001f, "out/error.csv");

    // TODO: solve for eigenvectors

    return 0;
}

/*
TODO:
- link to this: https://towardsdatascience.com/the-mathematics-behind-principal-component-analysis-fff2d7f4b643


PCA Algorithm:
* part of it is:
* get eigenvectors & eigenvalues from covariance matrix using QR algorithm with grant schmidt process
 * Gram Schmidt: https://en.wikipedia.org/wiki/Gram%E2%80%93Schmidt_process
 * QR Decomp: https://en.wikipedia.org/wiki/QR_decomposition
 * QR algorithm: https://en.wikipedia.org/wiki/QR_algorithm
 ! There are other eigenvector/value algorithms! and also other choices to make in QR decomp.
*/