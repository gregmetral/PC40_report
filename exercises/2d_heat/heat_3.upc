#include <upc_relaxed.h>
#include <bupc_collectivev.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>

#define N 30
#define BLOCKSIZE ((N+2) * (N+2) / THREADS)

shared[BLOCKSIZE] double grid[N+2][N+2];
shared[BLOCKSIZE] double new_grid[N+2][N+2];
shared double dTmax[THREADS];

shared[BLOCKSIZE] double *ptr[N+2];
shared[BLOCKSIZE] double *new_ptr[N+2];
shared[BLOCKSIZE] double *tmp_ptr;


void initialize(void)
{
    int j;

    for( j=1; j<N+1; j++ )
    {
        grid[0][j] = 1.0;
        new_grid[0][j] = 1.0;
    }
}

int main(void)
{
    struct timeval ts_st;
	struct timeval ts_end;
	
	double time_;
    double dTmax;
	double dT;
	double dTmax_final;
	double T;
	
    int i, j, k, l;



    if( MYTHREAD == 0 )
        initialize();

    for( i=0; i<N+2; i++ )
    {
        ptr[i] = grid[i];
        new_ptr[i] = new_grid[i];
    }


    double epsilon  = 0.0001;
    int finished = 0;
    int nr_iter = 0;

    upc_barrier;

    gettimeofday( &ts_st, NULL );

    do
    {
        dTmax = 0.0;
        for( i=1; i< N+1; i++)
        {
            upc_forall( j=1; j<N+1; j++; &grid[i][j])
            {
                T = 0.25 * (ptr[i+1][j] + ptr[i-1][j] + ptr[i][j-1] + ptr[i][j+1]);
                dT = fabs(T - ptr[i][j]);
                new_ptr[i][j] = T;
                if( dTmax < dT )
                    dTmax = dT;
            }
        }

		upc_barrier;
		dTmax_final = bupc_allv_reduce_all(double, dTmax, UPC_MAX); 

        if( dTmax_final < epsilon )
            finished = 1;
        else
        {
            for( k=0; k<N+2; k++ )
            {
                tmp_ptr    = ptr[k];
                ptr[k]     = new_ptr[k];
                new_ptr[k] = tmp_ptr;
            }
        }
        nr_iter++;
    } while( finished == 0 );

    gettimeofday( &ts_end, NULL );

    time_ = ts_end.tv_sec + (ts_end.tv_usec / 1000000.0);
    time_ -= ts_st.tv_sec + (ts_st.tv_usec / 1000000.0);

    printf("THREAD numb %d : %d iterations in %.3lf sec\n",MYTHREAD, nr_iter, time_);

    return 0;
}

