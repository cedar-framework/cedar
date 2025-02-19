/****************************************************/
/* MSG2.0 C functions and wallclock timer           */
/* written by Andrei Malevsky, version May 16, 1997 */
/*

   Reaord of changes for this version:

   two functions (MSG_SAVE_ARRAY and MSG_RESTORE_ARRAY) have
   been added April 7, 1997;
   the wrappers for MSG_TSETUP and MSG_TP_SETUP
   have been added
   and the function MSG_RESTORE_ARRAY has been modified
   April 22, 1997;
   include <stdlib.h> was added May 16, 1997
*/

#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>

#include "FC.h"

#define MAX_TIMERS 16

void msg_tsetup_int(int *, int *, int *, int *, int *, int *, int *, int *,
                    int *, int *, int *, int *, int *, int *, int *, int *,
                    int *, int *);
void msg_tp_setup_int(int *, int *, int *, int *, int *, int *, int *, int *,
                      int *, int *, int *, int *, int *, int *);


struct timeval start_tp[MAX_TIMERS] ;
struct timeval current_tp[MAX_TIMERS] ;
double myclock[MAX_TIMERS] ;
int *StorePtr;

void msg_tsetup(numproc,myproc,ptrn,
                grid_size,proc_size,overlap,ifp,
                nproc,proc,ipr,index,sfa,pfa,ier)
int *numproc, *myproc, *ptrn,
    *grid_size, *proc_size, *overlap, *ifp,
    *nproc, *proc, *ipr, *index, *sfa, *pfa, *ier;
{
   int nbytes, *gc_ld, *gc_eid, *la_size, *eid_s;
   nbytes = *numproc * 6 * sizeof(int);
   gc_ld = (int *) malloc(nbytes);
   gc_eid = (int *) malloc(nbytes);
   nbytes = *numproc * 3 * sizeof(int);
   la_size = (int *) malloc(nbytes);
   eid_s = (int *) malloc(nbytes);
   msg_tsetup_int(numproc,myproc,ptrn,
                grid_size,proc_size,overlap,ifp,
                nproc,proc,ipr,index,sfa,pfa,ier,
                gc_ld, gc_eid, la_size, eid_s);
   free(gc_ld);
   free(gc_eid);
   free(gc_eid);
   free(la_size);
   free(eid_s);
}


void msg_tp_setup(la_size, eid_s, gc_ld, gc_eid,
                  numproc, myproc, nproc, proc, ipr, index,
                  sfa, pfa, ier)
int *la_size, *eid_s, *gc_ld, *gc_eid,
    *numproc, *myproc, *nproc, *proc, *ipr, *index,
    *sfa, *pfa, *ier;
{
   int nbytes, *Periodic;
   nbytes = *numproc * 3 * sizeof(int);
   Periodic = (int *) malloc(nbytes);
   msg_tp_setup_int(la_size, eid_s, gc_ld, gc_eid,
                    numproc, myproc, nproc, proc, ipr, index,
                    sfa, pfa, ier, Periodic);
   free(Periodic);
}


void msg_save_array (Array, ArraySize)
int *Array, *ArraySize;
{
   int nbytes, i, *TmpPtr, *ArrayPtr;

   nbytes = *ArraySize * sizeof(int);
   StorePtr = (int *) malloc(nbytes);
   TmpPtr = StorePtr;
   ArrayPtr = Array;
   for (i=0; i<*ArraySize; i++)
   {
    *TmpPtr = *ArrayPtr;
    TmpPtr++; ArrayPtr++;
   }
}

void msg_restore_array (Array, ArraySize)
int *Array, *ArraySize;
{
   int i, *TmpPtr, *ArrayPtr;
   TmpPtr = StorePtr;
   ArrayPtr = Array;
   for (i=0; i<*ArraySize; i++)
   {
    *ArrayPtr = *TmpPtr;
    TmpPtr++; ArrayPtr++;
   }
   free (StorePtr);
}

void msg_timer_clear (timer)
int *timer;
{
     int res ;
     if ( *timer < 0 || *timer > MAX_TIMERS )
     {fprintf(stderr,"Wrong timer specified:%d\n",*timer);exit(1);}
     myclock[*timer] = 0.0;

}


void msg_timer_start (timer)
int *timer;
{
     int res ;
     if ( *timer < 0 || *timer > MAX_TIMERS )
     {fprintf(stderr,"Wrong timer specified:%d\n",*timer);exit(1);}
     res = gettimeofday( &start_tp[*timer] , NULL ) ;
     if ( res<0 )
     {fprintf(stderr,"Bad gettimeofday at start_wallclock_timer\n");exit(1);}
#ifdef _DEBUG
     printf("timer: %d\n",*timer);
     printf("start_tp.tv_sec  = %ld\n",start_tp[*timer].tv_sec );
     printf("start_tp.tv_usec = %ld\n",start_tp[*timer].tv_usec );
#endif

}

void msg_timer_stop (timer)
int *timer;
{
     int res ;
     double a , b ;
     if ( *timer < 0 || *timer > MAX_TIMERS )
     {fprintf(stderr,"Wrong timer specified:%d\n",*timer);exit(1);}
     res = gettimeofday( &current_tp[*timer] , NULL ) ;
     if ( res<0 )
     {fprintf(stderr,"Bad gettimeofday at stop_wallclock_timer\n");exit(1);}
     a = (double) (current_tp[*timer].tv_sec) * 1e6
         + (double) current_tp[*timer].tv_usec ;
     b = (double) (start_tp[*timer].tv_sec  ) * 1e6
         + (double) start_tp[*timer].tv_usec ;
     a = (a - b) * 1e-6 ;
     myclock[*timer] += a;
#ifdef _DEBUG
     printf("timer: %d\n",*timer);
     printf("current_tp.tv_sec  = %ld\n",current_tp[*timer].tv_sec  );
     printf("current_tp.tv_usec = %ld\n",current_tp[*timer].tv_usec );
#endif

}

double msg_timer (timer)
int *timer;
{
    if ( *timer < 0 || *timer > MAX_TIMERS )
    {fprintf(stderr,"Wrong timer specified:%d\n",*timer);exit(1);}
    return myclock[*timer];

}

/*  -------------------- test -------------------------------- */
#ifdef _DEBUG
main()
{
  double tribble ;
  int timer = MAX_TIMERS;

  msg_timer_clear (&timer) ;
  msg_timer_start (&timer) ;
  sleep ( 5 ) ;
  msg_timer_stop (&timer) ;
  sleep ( 7 ) ;
  msg_timer_start (&timer) ;
  sleep ( 3 ) ;
  msg_timer_stop (&timer) ;

  tribble = msg_timer (&timer);
  printf("Elapsed time should be about 8 seconds, is reported to be %f\n",
  tribble );

}

#endif
