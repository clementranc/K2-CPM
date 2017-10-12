#ifndef __CPMFORF_H_
#define __CPMFORF_H_

#include "matrix.h"

void linear_least_squares(Table*, Table*, const Table*, const Table*, double*, Table*);
void fit_target(const Table&, Table&, const Table&, double*, Table&);

#ifdef __cplusplus
    extern"C" { 
#endif
void run_cpm_(double*, double*, double*, double*, double*, int&, int&, 
int&, int&, double&, double&, double&, int&, char*, double*, double*, double*);
#ifdef __cplusplus
    }
#endif

#endif
