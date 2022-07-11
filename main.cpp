#include "vecmath.h"
#include <stdio.h>
#include <algorithm>
#include <direct.h>

template <size_t WIDTH, size_t HEIGHT>
using Mtx = BaseMtx<double, WIDTH, HEIGHT>;

template <size_t SIZE>
using Vec = BaseVec<double, SIZE>;

template <size_t DATA_WIDTH, size_t DATA_HEIGHT>
void Report(
    const Mtx<DATA_WIDTH, DATA_HEIGHT>& data,
    const Mtx<DATA_WIDTH, DATA_HEIGHT> recoveredData[],
    const Vec<DATA_WIDTH>& eigenValues,
    const std::array<Vec<DATA_WIDTH>, DATA_WIDTH>& eigenVectors,
    const Vec<DATA_WIDTH>& RMSE,
    const char* fileNameBase)
{
    printf("%s.txt\n", fileNameBase);

    char fileName[1024];
    sprintf_s(fileName, "out/%s.txt", fileNameBase);
    FILE* file = nullptr;
    fopen_s(&file, fileName, "w+t");

    fprintf(file, "Eigen values and eigen vectors\n\n");
    for (int i = 0; i < DATA_WIDTH; ++i)
    {
        fprintf(file, "[%f] [", eigenValues[i]);
        for (int j = 0; j < DATA_WIDTH; ++j)
        {
            if (j > 0)
                fprintf(file, ", %f", eigenVectors[i][j]);
            else
                fprintf(file, "%f", eigenVectors[i][j]);
        }
        fprintf(file, "]\n");
    }
    fprintf(file, "\nRMSE for component counts:\n");

    for (int componentCount = 0; componentCount < DATA_WIDTH; ++componentCount)
    {
        if (eigenValues[componentCount] == 0.0f)
            continue;
        fprintf(file, "[%i] %f\n", componentCount + 1, RMSE[componentCount]);
    }
    fprintf(file, "\n");

    fprintf(file, "Data:\n");
    for (int i = 0; i < DATA_HEIGHT; ++i)
    {
        fprintf(file, "[");
        for (int j = 0; j < DATA_WIDTH; ++j)
        {
            if (j > 0)
                fprintf(file, ", %f", data[i][j]);
            else
                fprintf(file, "%f", data[i][j]);
        }
        fprintf(file, "]\n");
    }

    for (int componentCount = 0; componentCount < DATA_WIDTH; ++componentCount)
    {
        if (eigenValues[componentCount] == 0.0f)
            continue;

        fprintf(file, "\n========================\n%i Components\n========================\n\n", componentCount + 1);

        for (int i = 0; i < DATA_HEIGHT; ++i)
        {
            fprintf(file, "[");
            for (int j = 0; j < DATA_WIDTH; ++j)
            {
                if (j > 0)
                    fprintf(file, ", %f", recoveredData[componentCount][i][j]);
                else
                    fprintf(file, "%f", recoveredData[componentCount][i][j]);
            }
            fprintf(file, "]\n");
        }
    }

    fclose(file);

    if (DATA_WIDTH != 2)
        return;

    // make and run python file to show eigenvectors
    {
        printf("%s.py\n", fileNameBase);

        sprintf_s(fileName, "out/%s.py", fileNameBase);
        FILE* file = nullptr;
        fopen_s(&file, fileName, "w+t");

        fprintf(file,
            "import numpy as np\n"
            "import matplotlib.pyplot as plt\n"
            "\n"
            "fig = plt.figure()\n"
            "ax = fig.add_subplot(111)\n"
            "ax.set_aspect('equal', adjustable='box')\n"
            "\n"
        );

        static const char colors[2] = {'r', 'g'};
        for (int eigenIndex = 0; eigenIndex < DATA_WIDTH; ++eigenIndex)
            fprintf(file, "plt.arrow(0, 0, %f, %f, color='%c', label=\"Eigen value %f\")\n", eigenVectors[eigenIndex][0], eigenVectors[eigenIndex][1], colors[eigenIndex], eigenValues[eigenIndex]);

        fprintf(file,
            "plt.title(\"Eigen Vectors\")\n"
            "plt.xlim([-1.1, 1.1])\n"
            "plt.ylim([-1.1, 1.1])\n"
            "plt.legend()\n"
            "plt.tight_layout()\n"
            "\n"
            "fig.savefig(\"out/%s.png\")\n", fileNameBase);
        fclose(file);

        // run the python file
        sprintf_s(fileName, "python out/%s.py", fileNameBase);
        system(fileName);
    }

    // get the min and max of the recovered data, to frame the graphs
    Vec<DATA_WIDTH> rdmin, rdmax;
    std::fill(rdmin.begin(), rdmin.end(), FLT_MAX);
    std::fill(rdmax.begin(), rdmax.end(), -FLT_MAX);
    for (int componentCount = 0; componentCount < DATA_WIDTH; ++componentCount)
    {
        if (eigenValues[componentCount] == 0.0f)
            continue;

        for (int rowIndex = 0; rowIndex < DATA_HEIGHT; ++rowIndex)
        {
            for (int columnIndex = 0; columnIndex < DATA_WIDTH; ++columnIndex)
            {
                rdmin[columnIndex] = std::min(rdmin[columnIndex], recoveredData[componentCount][rowIndex][columnIndex]);
                rdmin[columnIndex] = std::min(rdmin[columnIndex], data[rowIndex][columnIndex]);

                rdmax[columnIndex] = std::max(rdmax[columnIndex], recoveredData[componentCount][rowIndex][columnIndex]);
                rdmax[columnIndex] = std::max(rdmax[columnIndex], data[rowIndex][columnIndex]);
            }
        }
    }

    // expand the min/max box a bit
    for (int columnIndex = 0; columnIndex < DATA_WIDTH; ++columnIndex)
    {
        double mid = (rdmin[columnIndex] + rdmax[columnIndex]) / 2.0f;
        double halfWidth = (rdmax[columnIndex] - rdmin[columnIndex]) / 2.0f;
        rdmin[columnIndex] = mid - halfWidth * 1.1f;
        rdmax[columnIndex] = mid + halfWidth * 1.1f;
    }

    // spit out the python scripts to make graphs
    for (int componentCount = 0; componentCount < DATA_WIDTH; ++componentCount)
    {
        if (eigenValues[componentCount] == 0.0f)
            continue;

        printf("%s.%i.py\n", fileNameBase, componentCount + 1);

        sprintf_s(fileName, "out/%s.%i.py", fileNameBase, componentCount + 1);
        FILE* file = nullptr;
        fopen_s(&file, fileName, "w+t");
        fprintf(file,
            "import matplotlib.pyplot as plt\n"
            "import numpy as np\n"
            "\n"
        );

        // Print out eigen values and eigen vectors
        fprintf(file, "# Eigen Values & Eigen (Principle Component) Vectors:\n");
        for (int eigenIndex = 0; eigenIndex <= componentCount; ++eigenIndex)
        {
            fprintf(file, "# [%f] : [", eigenValues[eigenIndex]);
            for (int componentIndex = 0; componentIndex < eigenVectors[eigenIndex].size(); ++componentIndex)
            {
                if (componentIndex > 0)
                    fprintf(file, ", %f", eigenVectors[eigenIndex][componentIndex]);
                else
                    fprintf(file, "%f", eigenVectors[eigenIndex][componentIndex]);
            }
            fprintf(file, "]\n");
        }
        fprintf(file, "\n");

        // data
        for (int columnIndex = 0; columnIndex < DATA_WIDTH; ++columnIndex)
        {
            fprintf(file, "data%i = [", columnIndex);
            for (int rowIndex = 0; rowIndex < DATA_HEIGHT; ++rowIndex)
            {
                if (rowIndex > 0)
                    fprintf(file, ", %f", data[rowIndex][columnIndex]);
                else
                    fprintf(file, "%f", data[rowIndex][columnIndex]);
            }
            fprintf(file, "]\n");
        }

        // recoveredData
        fprintf(file, "\n");
        for (int columnIndex = 0; columnIndex < DATA_WIDTH; ++columnIndex)
        {
            fprintf(file, "recoveredData%i = [", columnIndex);
            for (int rowIndex = 0; rowIndex < DATA_HEIGHT; ++rowIndex)
            {
                if (rowIndex > 0)
                    fprintf(file, ", %f", recoveredData[componentCount][rowIndex][columnIndex]);
                else
                    fprintf(file, "%f", recoveredData[componentCount][rowIndex][columnIndex]);
            }
            fprintf(file, "]\n");
        }

        fprintf(file,
            "\n"
            "fig = plt.figure()\n"
            "ax = fig.add_subplot(111)\n"
            "ax.set_aspect('equal', adjustable='box')\n"
            "\n"
            "plt.scatter(data0, data1, label=\"Points\")\n"
            "plt.scatter(recoveredData0, recoveredData1, marker=\"+\", label=\"PCA Points\")\n"
            "\n"
        );
        fprintf(file,
            "plt.title(\"%i of %i PCA Components\\nRMSE = %f\")\n"
            "plt.xlim([%f, %f])\n"
            "plt.ylim([%f, %f])\n"
            "plt.legend()\n"
            "plt.tight_layout()\n"
            "\n",
            componentCount + 1, (int)DATA_WIDTH,
            RMSE[componentCount],
            rdmin[0], rdmax[0],
            rdmin[1], rdmax[1]
        );
        fprintf(file, "fig.savefig(\"out/%s.%i.png\")\n", fileNameBase, componentCount + 1);
        fclose(file);

        // run the python file
        sprintf_s(fileName, "python out/%s.%i.py", fileNameBase, componentCount + 1);
        system(fileName);
    }
}

template <size_t DATA_WIDTH, size_t DATA_HEIGHT>
void DoTest(Mtx<DATA_WIDTH, DATA_HEIGHT> data, const char* fileNameBase, bool centerData)
{
    // Calculate the mean of the data
    Vec<DATA_WIDTH> mean;
    for (int i = 0; i < Columns(data); ++i)
    {
        mean[i] = 0.0f;
        for (int j = 0; j < Rows(data); ++j)
            mean[i] = Lerp(mean[i], data[j][i], 1.0 / double(j + 1));
    }

    // center the data if we should
    Vec<DATA_WIDTH> dataOffset{};
    if (centerData)
    {
        for (int j = 0; j < Rows(data); ++j)
            data[j] = data[j] - mean;
        dataOffset = mean * -1.0;
        mean = Vec<DATA_WIDTH>{};
    }

    // Make the covariance matrix
    Mtx<DATA_WIDTH, DATA_WIDTH> covariance;
    for (int i = 0; i < Columns(covariance); ++i)
    {
        for (int j = 0; j < Rows(covariance); ++j)
        {
            covariance[i][j] = 0.0f;
            for (int k = 0; k < Rows(data); ++k)
                //covariance[i][j] += (data[k][i]) * (data[k][j]) / float(Rows(data));
                covariance[i][j] += (data[k][i] - mean[i]) * (data[k][j] - mean[j]) / float(Rows(data));
        }
    }

    // Get the eigenvalues
    Vec<Columns(covariance)> eigenValues = QRAlgorithm(covariance, 10000, 0.0001);

    // Sort from largest to smallest
    std::sort(eigenValues.begin(), eigenValues.end(), [](double a, double b) {return b < a; });

    // Solve for eigenvectors for each eigenvalue
    std::array<Vec<Columns(covariance)>, Columns(covariance)> eigenVectors;
    for (int i = 0; i < eigenValues.size(); ++i)
    {
        if (eigenValues[i] == 0.0f)
        {
            eigenVectors[i] = Vec<Columns(covariance)>{};
            continue;
        }

        eigenVectors[i] = InverseIteration(covariance, (float)eigenValues[i], 10000);

        // Get a better eigen value
        Vec<Columns(covariance)> test = Multiply(covariance, eigenVectors[i]);
        for (int j = 0; j < Columns(covariance); ++j)
        {
            if (eigenVectors[i][j] != 0.0f)
            {
                eigenValues[i] = test[j] / eigenVectors[i][j];
                break;
            }
        }
    }

    // Transform the data into the PCA space
    Mtx<Columns(covariance), Columns(covariance)> WT;
    for (int i = 0; i < Columns(WT); ++i)
        SetColumn(WT, i, eigenVectors[i]);
    Mtx<DATA_WIDTH, DATA_HEIGHT> pcaData = Multiply(data, WT);

    // un center the data if we centered it before
    if (centerData)
    {
        for (int j = 0; j < Rows(data); ++j)
            data[j] = data[j] - dataOffset;
    }

    // Recover the data back from the pca data
    Mtx<DATA_WIDTH, DATA_HEIGHT> recoveredData[Columns(WT)];
    Vec<Columns(WT)> RMSE{};
    for (int componentCount = 0; componentCount < Columns(WT); ++componentCount)
    {
        RMSE[componentCount] = 0.0;

        recoveredData[componentCount] = Mtx<DATA_WIDTH, DATA_HEIGHT>{};

        for (int rowIndex = 0; rowIndex < DATA_HEIGHT; ++rowIndex)
        {
            recoveredData[componentCount][rowIndex] = dataOffset * -1.0;
            for (int component = 0; component <= componentCount; ++component)
            {
                recoveredData[componentCount][rowIndex] = recoveredData[componentCount][rowIndex] + pcaData[rowIndex][component] * eigenVectors[component];
            }

            double squaredError = LengthSq(data[rowIndex] - recoveredData[componentCount][rowIndex]);
            RMSE[componentCount] = Lerp(RMSE[componentCount], squaredError, 1.0 / double(rowIndex + 1));
        }
        RMSE[componentCount] = sqrt(RMSE[componentCount]);
    }

    Report(
        data,
        recoveredData,
        eigenValues,
        eigenVectors,
        RMSE,
        fileNameBase);
}

int main(int argc, char** argv)
{
    _mkdir("out");

    // box blur - the covariance matrix is all zeros so is strange
    {
        Mtx<3, 3> data;

        for (auto& i : data)
            for (auto& j : i)
                j = 1.0f / float(Rows(data) * Columns(data));

        DoTest(data, "boxU", false);
        DoTest(data, "boxC", true);
    }

    // gaussian blur. Sigma 0.3
    // 1D kernel is 0.0478,	0.9044,	0.0478
    // http://demofox.org/gauss.html
    {
        Mtx<3, 3> data =
        {
            0.0023,	0.0432,	0.0023,
            0.0432,	0.8180,	0.0432,
            0.0023,	0.0432,	0.0023
        };

        DoTest(data, "gaussU", false);
        DoTest(data, "gaussC", true);
    }

    // test from the website
    {
        Mtx<3, 5> data =
        {
            90, 60, 90,
            90, 90, 30,
            60, 60, 60,
            60, 60, 90,
            30, 30, 30
        };

        DoTest(data, "test3DU", false);
        DoTest(data, "test3DC", true);
    }

    // Basic data test
    {
        Mtx<2, 4> data =
        {
            90, 60,
            90, 90,
            60, 60,
            30, 30,
        };

        DoTest(data, "testU", false);
        DoTest(data, "testC", true);
    }

    // Test centering vs not.
    // This also shows how it assumes the data is centered. it doesn't have a second axis to fit the data perfectly.
    {
        Mtx<2, 3> data =
        {
            0, 10,
            10, 11,
            20, 12
        };

        DoTest(data, "centerTest1U", false);
        DoTest(data, "centerTest1C", true);
    }

    // Same, but now they aren't colinear
    {
        Mtx<2, 3> data =
        {
            0, 10,
            10, 11,
            20, 20
        };

        DoTest(data, "centerTest2U", false);
        DoTest(data, "centerTest2C", true);
    }
}
