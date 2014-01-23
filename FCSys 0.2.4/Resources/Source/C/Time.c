// Prevent multiple inclusion.
#ifndef _FCSYS_TIME_INCLUDED_
#define _FCSYS_TIME_INCLUDED_

//#include <stdio.h>
#include <time.h>

int get_time(void)
{
    //printf ("%ld seconds since January 1, 1970", time(NULL));
    return (int) time(NULL);
}

#endif /* _FCSYS_TIME_INCLUDED_ */
