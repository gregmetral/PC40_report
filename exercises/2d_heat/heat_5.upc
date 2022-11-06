#include <upc_relaxed.h>
#include <upc_collective.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <stdlib.h>

shared double dTmax_final;

void initialize(shared[] double* grid, shared[] double* new_grid, int n) {
    for (int i = 0; i < n+2; i++) {
        for (int j = 0; j < n+2; j++) {
            grid[i * (n+2) + j] = 0.0;
            new_grid[i * (n+2) + j] = 0.0;
        }
    }

    for (int j = 1; j < n+1; j++) {
        grid[j] = 1.0;
        new_grid[j] = 1.0;
    }
}

int main(int argc, char* argv[]) {


    int n = atoi(argv[1]);
	
    struct timeval ts_st;
	struct timeval ts_end;

    shared[] double* dTmax = upc_all_alloc(1, THREADS);
    shared[] double* ptr = upc_all_alloc((n+2) * (n+2) / THREADS, (n+2) * (n+2));
    shared[] double* new_ptr = upc_all_alloc((n+2) * (n+2) / THREADS, (n+2) * (n+2));

	double T;
	double dT;

    if (MYTHREAD == 0) {
        initialize(ptr, new_ptr, n);
    }

    upc_barrier;
    int finished = 0;
	double time_;

	double epsilon  = 0.0001;
    int nr_iter = 0;
	int i, j, max_i;

    gettimeofday( &ts_st, NULL );
    do {
        dTmax[MYTHREAD] = 0.0;
        max_i = (n+2) * (MYTHREAD + 1) / THREADS;
        if (max_i > n+1) 
			max_i = n+1;

        for (i = ((n+2) * MYTHREAD / THREADS) ; i < max_i; i++) {
            for (j = 1; j <= n; j++) {
                T = 0.25 * (ptr[(i+1) * (n+2) + j] + ptr[(i-1) * (n+2) + j] + ptr[i * (n+2) + (j-1)] + ptr[i * (n+2) + (j+1)]);
                dT = fabs(T - ptr[i * (n+2) + j]);
                new_ptr[i * (n+2) + j] = T;
                if (dTmax[MYTHREAD] < dT) 
					dTmax[MYTHREAD] = dT;
            }
        }

        upc_all_reduceD(&dTmax_final, dTmax, UPC_MAX, THREADS, 1, NULL, UPC_IN_ALLSYNC | UPC_OUT_ALLSYNC);

        if (dTmax_final < epsilon) 
		{ 
            finished = 1;
        } else {
		
            shared[] double* tmp = ptr;
            ptr = new_ptr;
            new_ptr = tmp;
			
        }
        nr_iter++;

    } while (finished == 0);
    
	gettimeofday( &ts_end, NULL );

    time_ = ts_end.tv_sec + (ts_end.tv_usec / 1000000.0);
    time_ -= ts_st.tv_sec + (ts_st.tv_usec / 1000000.0);

    printf("THREAD numb %d : %d iterations in %.3lf sec\n",MYTHREAD, nr_iter, time_);



	return 0;
}