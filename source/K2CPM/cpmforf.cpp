//==================================================================//
// Author: Clement Ranc
//==================================================================//
//
// Name of the code: CPM+LCM
//
// This code is a C++ adaptation of the K2-CPM code [1][2].
// This version is written to work with Bennett's code [3].
//
// References
// ----------
// [1] Wang, D., Hogg, D. W., Foreman-Mackey, D. & Schölkopf, B. A Causal,
//     Data-driven Approach to Modeling the Kepler Data. Publications of the
//     Astronomical Society of the Pacific 128, 94503 (2016).
// [2] https://github.com/jvc2688/K2-CPM
// [3] Bennett, D. 2010, An Efficient Method for Modeling High-magnification 
//     Planetary Microlensing Events, ApJ, 716, 1408-1422.
//
//==================================================================//

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <stdio.h>
#include <ctype.h>
#include <sstream>
#include <ctime>

#include "cpmforf.h"

using namespace std;

//==================================================================//
// Functions
//==================================================================//
void linear_least_squares(Table* a, Table* y, const Table* yvar,
    const Table* l2_tab, double* tprpr, Table* result){
/*
    Solver of linear systems using Cholesky's decomposition. Let's define
    a matrix A (dimension n1 x n2) and two vectors X (dimension n2) and Y
    (dimension n1). This function solves the equation AX = Y. The
    following steps are considered:

    1/ Add the observational uncertainties to the data.

    2/ Compute the square matrix A^T A and the vector A^T Y.

    3/ Find the square matrix L so that A^T A = L L^T, where L^T is the
    transposition of L (L is triangular and lower).

    4/ Find the solution X of the system A^T A X = L L^T X = A^T Y.

    The solver assumes that the matrix A^T A is a square matrix, symmetric
    and positive-definite.

    Inputs
    ------
    a -- Table *, dimension n_data x n_predictors.
        The basis matrix. Will be overwritten.
    y -- Table *, dimension n_data.
        The observations. Will be overwritten.
    yvar -- Table *, dimension n_data.
        The observational variance of the points y.
    l2 -- Table *, dimension n_predictors.
        The L2 regularization strength.
    tprpr -- double*, dimension (n_pre-1) x (n_pre-1)
        Auxiliary array that does not change for a given pixel and
        given predictors.
    result -- Table *, dimension n_predictors.
        The solution will be written in this Table.
*/

    // Declarations and initialization
    // -------------------------------
    int i, j, k, dim1a, dim2a, x;
    double s;

    dim1a = a->get_size1();
    dim2a = a->get_size2();
    x = y->get_size1();
    assert(dim1a == x);
    x = y->get_size2();
    assert(x == 1);
    x = y->get_size3();
    assert(x == 1);

    x = yvar->get_size1();
    assert(dim1a == x);
    x = yvar->get_size2();
    assert(x == 1);
    x = yvar->get_size3();
    assert(x == 1);

    x = l2_tab->get_size1();
    assert(dim2a == x);
    x = l2_tab->get_size2();
    assert(x == 1);
    x = l2_tab->get_size3();
    assert(x == 1);

    Matrix ata(dim2a);
    Table cia(dim1a, dim2a), at(dim2a, dim1a), ciy(dim1a), b(dim2a);

    // Incorporate observational uncertainties
    // ---------------------------------------
    for(i=0; i<dim1a; ++i){
        for(j=0; j<dim2a; ++j){
            cia.set(i, j) = (*a)(i, j) / (*yvar)(i);
            at.set(j, i) = (*a)(i, j);  // compute transpose of a
        }
    }
    ciy = (*y) / (*yvar);

    // Compute the pre-factor
    // ----------------------
    for(i = 0; i < dim2a; ++i){
        for(j = 0; j < dim2a; ++j){
            if((i<dim2a-1) && (j<dim2a-1)) ata.set(i,j) = tprpr[i*(dim2a-1)+j];
            else{
                s = 0;
                for(k = 0; k < dim1a; ++k) s += at(i, k) * cia(k, j);
                ata.set(i, j) = s;
            }
        }
    }

    for(i=0; i<dim2a; ++i){
        s = 0;
        for(j=0; j<dim1a; ++j) s += at(i, j) * ciy(j);
        b.set(i) = s;
    }

    // Incorporate any L2 regularization
    // ---------------------------------
    for(i = 0; i < dim2a; ++i) { ata.set(i, i) += (*l2_tab)(i); }

    // Solve the equations overwriting the matrix and tables
    // -----------------------------------------------------
    ata.cholesky_factor();
    ata.cholesky_solve(b);

    *result = b;
}
//==================================================================//
void fit_target(const Table& tpf_timeserie, Table& pre_matrix,
    const Table& l2_tab, double* tprpr, Table& result){
/*
    Find the best combination of the predictors to describe the pixel
    observations.

    Input
    -----
    tpf_timeserie -- Table &, dimension (n_dates x 3).
        First column is the date, second column is the target flux, third
        column is the error on the flux.
    pre_matrix -- Table &, dimension n_dates x n_pre.
        The flux of nearby stars used in the fitting process. Here, n_pre is
        the predictors number plus the polynomial order plus the model. 
        This Table is overwritten.
    l2_tab -- Table &, dimension n_pre.
        Array of L2 regularization strength.
    tprpr -- double*, dimension (n_pre-1) x (n_pre-1)
        Auxiliary array that does not change for a given pixel and
        given predictors.
    result -- Table &, dimension n_pre
        Result of the fit will be written in this Table.
*/

    // Declarations and initializations
    // --------------------------------
    int i, n_dates, n_pre;

    n_dates = tpf_timeserie.get_size1();
    n_pre = pre_matrix.get_size2();

    // Fit
    // ---
    Table y(n_dates), yvar(n_dates);
    for(i=0; i<n_dates; ++i) {
        y.set(i) = tpf_timeserie(i, 1);
        yvar.set(i) = pow(tpf_timeserie(i, 2), 2);
    }
    // yvar = 1.0; // EDIT

    linear_least_squares(&pre_matrix, &y, &yvar, &l2_tab, tprpr, &result);
}

//==================================================================//
void run_cpm_(double* pxdt, double* pxflx, double* pxsg, double* prflx,
double* tprpr, int& npxdt2, int& nprdt2, int& npr, int& npr2, double& l2_cpm,
double& sfcpm, double& chi2_cpm, int& fllc, char* fitname, double* cpmflx,
double* sgcpmflx, double* cpmrs){

    time_t start_cpm_part2_fit, end_cpm_part2_fit;
    cout << "CPM+LCM: C++ lib running..." << endl;

    // Declaration and initialisations
    // -------------------------------
    int i, j;
    int n_dates, n_pre, n_pre2, n_pre_dates;
    double x, l2, chi2;
    string fit_id;

    // Convert types from fortran
    // --------------------------
    n_dates = npxdt2;
    n_pre_dates = nprdt2;
    n_pre = npr;
    n_pre2 = npr2;
    l2 = l2_cpm;

    // Caution: here same number of date for predictors than for pixel.
    Table tpf_timeserie(n_dates, 3), ml_model(n_pre_dates);
    Table pre_matrix(n_pre_dates,n_pre2);

    for (i=0; i<n_dates; ++i){
        tpf_timeserie.set(i,0) = pxdt[i];
        tpf_timeserie.set(i,1) = pxflx[i];
        tpf_timeserie.set(i,2) = pxsg[i];
        for (j=0; j<n_pre2; ++j){
            pre_matrix.set(i,j)=prflx[i*n_pre2+j];
//            if(j==n_pre2-1) cout << pre_matrix(i,j) << " ";
        }
    }

    // Prepare regularization
    // ----------------------
    Table l2_tab(n_pre2);
    l2_tab = l2;
    if (n_pre < n_pre2) for(i=n_pre; i<n_pre2; ++i) l2_tab.set(i) = 0.0;

    // Fit target
    // ----------
    clock_t sfit = clock();
    Table result(n_pre2);
    fit_target(tpf_timeserie, pre_matrix, l2_tab, tprpr, result);
    Table flux_fit(n_dates);
    for(i=0; i<n_dates; ++i){
        x = 0;
        for(j=0; j<n_pre2; ++j) x += pre_matrix(i, j) * result(j);
        flux_fit.set(i) = x;
    }

    // Resulting light-curve and chi2
    // ------------------------------
    Table dif(n_dates);
    for(i=0; i<n_dates; ++i) dif.set(i) = tpf_timeserie(i, 1) - flux_fit(i);
        
    chi2 = 0;
    for (i=0; i<n_dates; ++i) {
        if(tpf_timeserie(i, 2)<1e-10) {
            cout << "Error equal to zero. Result infinite." << endl;
        }
        chi2 += pow(dif(i) / tpf_timeserie(i, 2), 2);
    }

    chi2_cpm = chi2;
    cout << setprecision(12);
    cout << "chi2 of K2 observations: " << chi2 << endl;

    sfcpm = result(n_pre2-1);

    // Save the light curve
    // --------------------
    // fllc=1;
    if(fllc){
        Table flux_fit2(n_dates);
        for(i=0; i<n_dates; ++i){
            x = 0;
            for(j=0; j<n_pre2-1; ++j) x += pre_matrix(i, j) * result(j);
            flux_fit2.set(i) = x;
        }

        Table dif2(n_dates);
        for(i=0; i<n_dates; ++i) {
            //dif2.set(i) = (tpf_timeserie(i, 1) - flux_fit2(i)) / result(n_pre2-1);
            dif2.set(i) = tpf_timeserie(i, 1) - flux_fit2(i);
            cpmrs[i] = dif(i);
            cpmflx[i] = dif2(i);
            sgcpmflx[i] = tpf_timeserie(i,2);
        }

//        fit_id = "resid.k2.";
//        for (i=0; i<80;++i){
//            if(fitname[i]==' ') continue;
//            else fit_id += fitname[i];  // Caution: the name is wrong, to much incremented.
//        } 
//        ofstream file1(fit_id);
//        if (file1.is_open()){
//            file1 << setprecision(8);
//            file1 << "# Source flux: " << result(n_pre2-1) << endl;
//            file1 << "# Date TPF_Flux Error_TPF_Flux CPM+LCM_Flux Error_CPM+LCM_Flux Residuals_Flux Magnification Error_Magnification Model_Magnification" << endl;
//            for (i=0; i<n_dates; ++i) {
//                file1 << fixed << setprecision(10);
//                file1 << tpf_timeserie(i,0)-2450000 << " ";
//                file1 << fixed << setprecision(5);
//                file1 << tpf_timeserie(i,1) << " ";
//                file1 << fixed << setprecision(3);
//                file1 << tpf_timeserie(i,2) << " ";
//                file1 << fixed << setprecision(5);
//                file1 << dif2(i) << " ";
//                file1 << fixed << setprecision(3);
//                file1 << tpf_timeserie(i, 2) << " ";
//                file1 << fixed << setprecision(8);
//                file1 << dif(i) << " ";
//                file1 << dif2(i)/result(n_pre2-1) << " ";
//                file1 << tpf_timeserie(i, 2)/result(n_pre2-1) << " ";
//                file1 << pre_matrix(i,n_pre2-1) << endl;
//            }
//            file1.close();
//        }
    }

    clock_t efit = clock();
    cout << "CPM+LCM in " << 1e-6*(efit-sfit) << " seconds. DONE." << endl;
}

