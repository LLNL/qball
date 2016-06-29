#ifndef __PROFILE_H__
#define __PROFILE_H__

#ifdef TAU
#include <Profile/Profiler.h>
#define TAU_FSTART(ARG)                                 \
    TAU_PROFILE_TIMER(timer##ARG, #ARG, "", TAU_USER);  \
    TAU_PROFILE_START(timer##ARG)

#define TAU_FSTOP(ARG)                                  \
    TAU_PROFILE_STOP(timer##ARG)
#define QB_Pstart(n, ARG) TAU_FSTART(ARG)
#define QB_Pstop(ARG) TAU_FSTOP(ARG)
#else
#define TAU_PROFILE(NAME,ARG,USER)
#define TAU_PROFILE_TIMER(ARG1, ARG2, ARG3, ARG4)
#define TAU_PROFILE_STOP(ARG)
#define TAU_PROFILE_START(ARG)
#define TAU_FSTART(ARG)
#define TAU_FSTOP(ARG)
#ifdef USE_MPIP
#define QB_Pstart(n, ARG) MPI_Pcontrol(n)
#define QB_Pstop(ARG) MPI_Pcontrol(0)
#else
#define QB_Pstart(n, ARG) 
#define QB_Pstop(ARG) 
#endif
#endif




#endif// __PROFILE_H__
