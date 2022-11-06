#include <upc_relaxed.h>
#include <bupc_collectivev.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>

#define N 30
#define BLOCKSIZE ((N+2) * (N+2) / THREADS)
#define LOCALROWS ((N+2) / THREADS)
#define LOCALSIZE (sizeof(double) * LOCALROWS * (N+2))

shared[BLOCKSIZE] double grid[N+2][N+2];
shared[BLOCKSIZE] double new_grid[N+2][N+2];
shared double dTmax[THREADS];

shared[BLOCKSIZE] double (*ptr)[N+2];
shared[BLOCKSIZE] double (*new_ptr)[N+2];



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
	
	double (*ptr_priv)[N+2] = malloc(LOCALSIZE);
	double (*new_ptr_priv)[N+2] = malloc(LOCALSIZE);

    int i, j, k, l;



    if( MYTHREAD == 0 )
        initialize();

	upc_barrier;

	upc_memget(ptr_priv, &grid[LOCALROWS * MYTHREAD], LOCALSIZE);
    upc_memget(new_ptr_priv, &new_grid[LOCALROWS * MYTHREAD], LOCALSIZE);

    double epsilon  = 0.0001;
    int finished = 0;
    int nr_iter = 0;

    gettimeofday( &ts_st, NULL );

    do
    {
        dTmax = 0.0;
		i = 0;
		
		if (i+(LOCALROWS*MYTHREAD) > 0) //first block
		{
			for( j=1; j<N+1; j++)
			{
				T = 0.25 * (ptr_priv[i+1][j] + ptr[(LOCALROWS*MYTHREAD)+i-1][j] + ptr_priv[i][j-1] + ptr_priv[i][j+1]);
				dT = fabs(T - ptr_priv[i][j]);
				new_ptr_priv[i][j] = T;
				if( dTmax < dT )
					dTmax = dT;
			}
		}
		
        for( i+=1; i< LOCALROWS-1 ; i++)//middle block
        {
            for( j=1; j<N+1; j++)
            {
                T = 0.25 * (ptr_priv[i+1][j] + ptr_priv[i-1][j] + ptr_priv[i][j-1] + ptr_priv[i][j+1]);
                dT = fabs(T - ptr_priv[i][j]);
                new_ptr_priv[i][j] = T;
                if( dTmax < dT )
                    dTmax = dT;
            }
        }
		
		if (i+(LOCALROWS*MYTHREAD) < N+1) //last block
		{
			for( j=1; j<N+1; j++)
			{
				T = 0.25 * (ptr[(LOCALROWS*MYTHREAD)+i+1][j] + ptr_priv[i-1][j] + ptr_priv[i][j-1] + ptr_priv[i][j+1]);
				dT = fabs(T - ptr_priv[i][j]);
				new_ptr_priv[i][j] = T;
				if( dTmax < dT )
					dTmax = dT;
			}
		}
		
		upc_barrier;

		upc_memput(&new_ptr[LOCALROWS * MYTHREAD], new_ptr_priv, LOCALSIZE);
		dTmax_final = bupc_allv_reduce_all(double, dTmax, UPC_MAX); 
		printf("%.31f et %.31f \n", epsilon, dTmax_final); 

        if( dTmax_final < epsilon )
            finished = 1;
        else
        {

            shared[BLOCKSIZE] double (*tmp_ptr)[N+2] = ptr;
            ptr = new_ptr;
            new_ptr = tmp_ptr;
				
			
	        double (*ptr_tmp_priv)[N+2] = ptr_priv;
            ptr_priv = new_ptr_priv;
            new_ptr_priv = ptr_tmp_priv;
        }
        nr_iter++;
    } while( finished == 0 );

    gettimeofday( &ts_end, NULL );

    time_ = ts_end.tv_sec + (ts_end.tv_usec / 1000000.0);
    time_ -= ts_st.tv_sec + (ts_st.tv_usec / 1000000.0);

    printf("THREAD numb %d : %d iterations in %.3lf sec\n",MYTHREAD, nr_iter, time_);

    return 0;
}

