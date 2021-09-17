#include <iostream>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include "math/gsl_helpers.h"
#include "Eigenvalue.h"

using namespace std;

Eigenvalue::Eigenvalue(gsl_matrix* matrix1, gsl_matrix* matrix2):
        Jacobianmatrix(matrix1),
        Hessiancartesian(matrix2),
        m(matrix1->size1),
        n(matrix1->size2),
        Hessiantorsionangle(gsl_matrix_alloc(n,n)),
        eigenvalue(gsl_vector_alloc (n)),
        eigenvector(gsl_matrix_alloc (n, n))
{
    gsl_matrix* Hessiantorsionangle1 = gsl_matrix_calloc(n,m);

    gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, Jacobianmatrix, Hessiancartesian, 0.0, Hessiantorsionangle1);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Hessiantorsionangle1, Jacobianmatrix, 0.0, Hessiantorsionangle);

    gsl_matrix_free(Hessiantorsionangle1);
    gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (n);

    gsl_eigen_symmv (Hessiantorsionangle, eigenvalue, eigenvector, w);

    gsl_eigen_symmv_free (w);

    gsl_eigen_symmv_sort (eigenvalue, eigenvector, GSL_EIGEN_SORT_ABS_ASC);
}


Eigenvalue::~Eigenvalue(){
    gsl_matrix_free(Jacobianmatrix);
    gsl_matrix_free(Hessiancartesian);
    gsl_matrix_free(Hessiantorsionangle);
    gsl_vector_free(eigenvalue);
    gsl_matrix_free(eigenvector);
}

gsl_matrix* Eigenvalue::getHessiantorsionangle() const{
    return Hessiantorsionangle;
}

gsl_vector* Eigenvalue::getEigenvalue() const{
    return eigenvalue;
}

