#ifndef __TIMER_H
#define __TIMER_H

#if defined(WIN32)
#include <olectl.h>
typedef struct _Timer {
	DWORD start;
} Timer;
void timer_start(Timer * t){
	SYSTEMTIME st;
	FILETIME ft;
	GetSystemTime(&st);
	SystemTimeToFileTime(&st,&ft);
	t->start= ft.dwLowDateTime;
}
double timer_elapsed(Timer * t){
	SYSTEMTIME st;
	FILETIME ft;
	DWORD res;
	GetSystemTime(&st);
	SystemTimeToFileTime(&st,&ft);
	res = ft.dwLowDateTime - t->start;
	return ((double) res)/10000000.0;
}

#else

#include <stdlib.h>
#include <sys/time.h>

typedef struct __Timer{
	struct timeval start;
} Timer;

#define timer_start(t) gettimeofday(&(t)->start, NULL) 
double timer_elapsed(Timer * t){
	struct timeval now;
	double res; 
	gettimeofday(&now,NULL);
	timersub(&now, &(t->start), &now);
	res = now.tv_sec + now.tv_usec/1000000.0;
	return res;
}

#endif

#define timer_clear(t) bzero(t, sizeof(Timer))

#endif //__TIMER_H
