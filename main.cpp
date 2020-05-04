#include <mpi.h>
#include <iostream>
#include <zconf.h>
#include <cmath>
#include <chrono>



int main(int argc, char* argv[])
{
    clock_t start = clock();

    auto start_p = std::chrono::high_resolution_clock::now();
    start_p = std::chrono::high_resolution_clock::now();
    int pid, np, elements_per_process, n_elements_recieved;

    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    std::cout<<np<<"\n";

    const long int n = 500000;
    elements_per_process = n / np;


    if (pid == 0) {

        double *v = (double *) malloc(n * sizeof(double));
        double *v1 = (double *) malloc(n * sizeof(double));
        double *v3 = (double *) malloc(n * sizeof(double));
        double *v4 = (double *) malloc(n * sizeof(double));
        double *v5 = (double *) malloc(n * sizeof(double));

        double norm = 0.;
        double norm1 = 0.;
        double norm3 = 0.;
        double norm4 = 0.;
        double norm5 = 0.;

        for (int i = 0; i < n; i++) {
            v[i] = cos(double(i)) * 10.;
            v1[i] = cos(double(i+1)) * 10.;
            v3[i] = cos(double(i+2)) * 10.;
            v4[i] = cos(double(i+3)) * 10.;
            v5[i] = cos(double(i+4)) * 10.;
        }

        int index, i;
        if (np > 1) {
            for (i = 1; i < np - 1; i++) {
                index = i * elements_per_process;

                MPI_Send(&elements_per_process,
                         1, MPI_INT,
                         i, 0,
                         MPI_COMM_WORLD);
                MPI_Send(&v[index],
                         elements_per_process,
                         MPI_DOUBLE, i, 0,
                         MPI_COMM_WORLD);
            }


            index = i * elements_per_process;
            int elements_left = n - index;

            MPI_Send(&elements_left,
                  1, MPI_INT,
                     i, 0,
                     MPI_COMM_WORLD);
            MPI_Send(&v[index],
                     elements_left,
                     MPI_DOUBLE, i, 0,
                     MPI_COMM_WORLD);
        }

        for (i = 0; i < elements_per_process; i++){
            norm += v[i] * v[i];
            norm1 += v1[i] * v1[i];
            norm3 += v3[i] * v3[i];
            norm4 += v4[i] * v4[i];
            norm5 += v5[i] * v5[i];}
        norm = sqrt(norm);
        norm1 = sqrt(norm1);
        norm3 = sqrt(norm3);
        norm4 = sqrt(norm4);
        norm5 = sqrt(norm5);



        double tmp;
        for (i = 1; i < np; i++) {
            MPI_Recv(&tmp, 1, MPI_DOUBLE,
                     MPI_ANY_SOURCE, 0,
                     MPI_COMM_WORLD,
                     &status);
            int sender = status.MPI_SOURCE;
            v[i] /= norm;
        }

        //printf("Sum of array is : %f\n", norm);
    }

    else {
        MPI_Recv(&n_elements_recieved,
                 1, MPI_INT, 0, 0,
                 MPI_COMM_WORLD,
                 &status);

        double v2[n_elements_recieved];

        MPI_Recv(v2, n_elements_recieved,
                 MPI_DOUBLE, 0, 0,
                 MPI_COMM_WORLD,
                 &status);

        double partial_sum = 0;
        for (int i = 0; i < n_elements_recieved; i++) {
            partial_sum += v2[i] * v2[i];
        }

        for (int i = 0; i < n_elements_recieved; i++){
            v2[i] /= partial_sum;
        }

        MPI_Send(&partial_sum, 1, MPI_DOUBLE,
                 0, 0, MPI_COMM_WORLD);
    }

    // cleans up all MPI state before exit of process
    MPI_Finalize();
    if (pid == 0) {
        clock_t end = clock();
        double seconds = (double)(end - start) / CLOCKS_PER_SEC;
        auto finish_p = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_p = (finish_p - start_p);

        std::cout << "Parallel time: " << elapsed_p.count() << " s\n";
    }
    return 0;
}


//int main(int argc, char **argv) {
//    auto start_p = std::chrono::high_resolution_clock::now();
//    start_p = std::chrono::high_resolution_clock::now();
//
//    const long int N = 1000000;
//	auto *v = (double*)malloc(N * sizeof(double));
//
//	//MPI - Start
//	int rank, size;
//	MPI_Init(&argc, &argv);
//    MPI_Status status;
//    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//	MPI_Comm_size(MPI_COMM_WORLD, &size);
//
//	const int nT = N / 4;
//    double vT[N], v1[nT];
//
//	if (rank == 0) {
//		for (double & i : vT)
//		{
//			i = cos(double(i)) * 10.;
//		}
//	}
//
//    int temp;
//    for (int i = 0; i < size - 1; i++) {
//        for (int j = 0; j < size - i - 1; j++) {
//            if (v1[j] > v1[j + 1]) {
//                // меняем элементы местами
//                temp = v1[j];
//                v1[j] = v1[j + 1];
//                v1[j + 1] = temp;
//            }
//        }
//    }
//
//	MPI_Scatter(&vT[0], nT, MPI_INT, &v1[0], nT, MPI_INT, 0, MPI_COMM_WORLD);
//
//	double norm = 0.;
//
//	// compute the norm of v
//	for (int i = 0; i < nT; i++)
//		norm += v1[i] * v1[i] ;
//
//
//	int a, dest, source;
//	for (int i = 1; i < size - 1; i *= 2)
//	{
//		a = 0;
//		dest = rank + i;
//		source = rank - i;
//		if (dest > size - 1) dest = MPI_PROC_NULL;
//		if (source < 0)  source = MPI_PROC_NULL;
//
//		MPI_Send(&norm, 1, MPI_INT, dest, i, MPI_COMM_WORLD);
//        MPI_Recv(&a, 1, MPI_INT, source, i, MPI_COMM_WORLD, &status);
//		norm = norm + a;
//	}
//	if (rank == size - 1) {
//		norm = sqrt(norm);
//	}
//
//
//
//	MPI_Scatter(&norm, 1, MPI_INT, &norm, 1, MPI_INT, size - 1, MPI_COMM_WORLD);
//
//	// normalize v
//	for (int i = 0; i < nT; i++)
//		v1[i] /= norm;    sleep(elapsed_p*10);

//
//
//
//	MPI_Gather(&v1[0], nT, MPI_INT, &vT[0], nT, MPI_INT, 0, MPI_COMM_WORLD);
//
//	if (rank == 0) {
//		auto finish_p = std::chrono::high_resolution_clock::now();
//		std::chrono::duration<double> elapsed_p = finish_p - start_p;
//		std::cout << "Parallel time: " << elapsed_p.count() << " s\n";
//    }
//
//	MPI_Finalize(); //MPI - Finish
//	return 0;
//}