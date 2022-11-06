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
	double dTmax_final;
	double T;
	double dT;
	double dTmax;

	int i, j, k, l;

	
	if(MYTHREAD == 0) {
	initialize();
	}
 
    upc_barrier;
	
    /* Set the precision wanted */
    double epsilon  = 0.0001;
    int finished = 0;
    int nr_iter = 0;

    /* and start the timed section */
    gettimeofday( &ts_st, NULL );

    do
    {
        dTmax = 0.0; 
        for( i=1; i< N+1; i++)
        {
            upc_forall( j=1; j<N+1; j++; &grid[i][j] )
            {
                T = 0.25 *(grid[i+1][j] + grid[i-1][j] + grid[i][j-1] + grid[i][j+1]); 
                dT = fabs(T - grid[i][j]); /* local variation */
                new_grid[i][j] = T;
                if( dTmax < dT ){
					dTmax = dT; /* max variation in this iteration */
				}        
            }
        }
		upc_barrier;
		dTmax_final = bupc_allv_reduce_all(double, dTmax, UPC_MAX); 
		printf("%.31f et %.31f \n", epsilon, dTmax_final); 
        
        if( dTmax_final < epsilon ){
			
            finished = 1;
		}
        else
        {
            for( k=0; k<N+2; k++){      /* not yet ... Need to prepare */
                upc_forall( l=0; l<N+2; l++; &grid[k][l] ){				/* ourselves for doing a new */
                    grid[k][l] = new_grid[k][l]; /* iteration */
                }   
			}	
        }
        nr_iter++;
        upc_barrier;
    } while( finished == 0 );

    gettimeofday( &ts_end, NULL ); /* end the timed section */
	

    time_ = ts_end.tv_sec + (ts_end.tv_usec / 1000000.0);
    time_ -= ts_st.tv_sec + (ts_st.tv_usec / 1000000.0);

    printf("THREAD numb %d : %d iterations in %.3lf sec\n",MYTHREAD, nr_iter, time_);

    return 0;
}
