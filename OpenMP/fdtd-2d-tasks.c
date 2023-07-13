/* Include benchmark-specific header. */
#include "fdtd-2d.h"
#include <omp.h>
double bench_t_start, bench_t_end;

static
double rtclock()
{
    struct timeval Tp;
    int stat;
    stat = gettimeofday (&Tp, NULL);
    if (stat != 0)
      printf ("Error return from gettimeofday: %d", stat);
    return (Tp.tv_sec + Tp.tv_usec * 1.0e-6);
}

void bench_timer_start()
{
  bench_t_start = rtclock ();
}

void bench_timer_stop()
{
  bench_t_end = rtclock ();
}

void bench_timer_print()
{
  printf ("Time in seconds = %0.6lf\n", bench_t_end - bench_t_start);
}


static
void init_array (int tmax,
   int nx,
   int ny,
   double ex[ nx][ny],
   double ey[ nx][ny],
   double hz[ nx][ny],
   double _fict_[ tmax])
{
  int i, j;
  
#pragma omp taskgroup
{  
  #pragma omp taskloop num_tasks(NUM_TASKS) private(i) nogroup
  for (i = 0; i < tmax; i++)
    _fict_[i] = (double) i;

  #pragma omp taskloop num_tasks(NUM_TASKS) private(i, j) nogroup collapse(2)
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++)
      {
      ex[i][j] = ((double) i*(j+1)) / nx;
      ey[i][j] = ((double) i*(j+2)) / ny;
      hz[i][j] = ((double) i*(j+3)) / nx;
      }
}
}

static
void print_array(int nx,
   int ny,
   double ex[ nx][ny],
   double ey[ nx][ny],
   double hz[ nx][ny])
{
  int i, j;

  fprintf(stderr, "==BEGIN DUMP_ARRAYS==\n");
  fprintf(stderr, "begin dump: %s", "ex");
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++) {
      if ((i * nx + j) % 20 == 0) fprintf(stderr, "\n");
      fprintf(stderr, "%0.2lf ", ex[i][j]);
    }
  fprintf(stderr, "\nend   dump: %s\n", "ex");
  fprintf(stderr, "==END   DUMP_ARRAYS==\n");

  fprintf(stderr, "begin dump: %s", "ey");
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++) {
      if ((i * nx + j) % 20 == 0) fprintf(stderr, "\n");
      fprintf(stderr, "%0.2lf ", ey[i][j]);
    }
  fprintf(stderr, "\nend   dump: %s\n", "ey");

  fprintf(stderr, "begin dump: %s", "hz");
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++) {
      if ((i * nx + j) % 20 == 0) fprintf(stderr, "\n");
      fprintf(stderr, "%0.2lf ", hz[i][j]);
    }
  fprintf(stderr, "\nend   dump: %s\n", "hz");
}

static
void kernel_fdtd_2d(int tmax,
      int nx,
      int ny,
      double ex[ nx][ny],
      double ey[ nx][ny],
      double hz[ nx][ny],
      double _fict_[ tmax])
{
  int t, i, j;

  for(t = 0; t < tmax; t++)
  {
    #pragma omp taskgroup
    {
    #pragma omp taskloop num_tasks(NUM_TASKS) private(j) nogroup
    for (j = 0; j < ny; j++)
      ey[0][j] = _fict_[t];

    for (i = 1; i < nx; i++)
    #pragma omp task
      for (j = 1; j < ny; j++)
      {
        ey[i][j] = ey[i][j] - 0.5*(hz[i][j] - hz[i - 1][j]);
        ex[i][j] = ex[i][j] - 0.5*(hz[i][j] - hz[i][j - 1]);
      }

    #pragma omp taskloop num_tasks(NUM_TASKS) private(i) nogroup
    for (i = 1; i < nx; i++)
      ey[i][0] = ey[i][0] - 0.5*(hz[i][0] - hz[i - 1][0]);
    
    #pragma omp taskloop num_tasks(NUM_TASKS) private(j) nogroup
    for (j = 1; j < ny; j++)
      ex[0][j] = ex[0][j] - 0.5*(hz[0][j] - hz[0][j - 1]);
    }

    #pragma omp taskgroup
    {
    for (i = 0; i < nx - 1; i++)
    #pragma omp task
      for (j = 0; j < ny - 1; j++)
        hz[i][j] = hz[i][j] - 0.7* (ex[i][j + 1] - ex[i][j] + ey[i + 1][j] - ey[i][j]);
    }
  }
}

int main(int argc, char** argv)
{
  int tmax = TMAX;
  int nx = NX;
  int ny = NY;
  double (*ex)[nx][ny];
  ex = (double(*)[nx][ny])malloc ((nx) * (ny) * sizeof(double));

  double (*ey)[nx][ny];
  ey = (double(*)[nx][ny])malloc ((nx) * (ny) * sizeof(double));

  double (*hz)[nx][ny];
  hz = (double(*)[nx][ny])malloc ((nx) * (ny) * sizeof(double));

  double (*_fict_)[tmax];
  _fict_ = (double(*)[tmax])malloc ((tmax) * sizeof(double));

#pragma omp parallel num_threads(4)
#pragma omp single
{
  #pragma omp task
  init_array (tmax, nx, ny,
       *ex,
       *ey,
       *hz,
       *_fict_);
       
  #pragma omp taskwait

  bench_timer_start();
  
  #pragma omp task
  kernel_fdtd_2d (tmax, nx, ny,
    *ex,
    *ey,
    *hz,
    *_fict_);
}

  bench_timer_stop();

  bench_timer_print();

  if (argc > 42 && ! strcmp(argv[0], "")) print_array(nx, ny, *ex, *ey, *hz);
                             
  free((void*)ex);
  free((void*)ey);
  free((void*)hz);
  free((void*)_fict_);

  return 0;
}
