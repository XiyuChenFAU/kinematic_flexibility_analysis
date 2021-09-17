#include <iostream>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sort_vector.h>
#include "math/gsl_helpers.h"
#include "Eigenvalue.h"
#include <math/math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

using namespace std;

Eigenvalue::Eigenvalue(gsl_matrix* matrix1, gsl_matrix* matrix2, gsl_matrix* matrix3):
        Jacobianmatrix(matrix1),
        Hessiancartesian(matrix2),
        Massmatrix(matrix3),
        m(matrix1->size1),
        n(matrix1->size2),
        Hessiantorsionangle(gsl_matrix_calloc(n,n)),
        eigenvaluealpha(gsl_vector_complex_calloc (n)),
        singularvalue(gsl_vector_calloc (n)),
        eigenvaluebeta(gsl_vector_calloc (n)),
        eigenvector(gsl_matrix_complex_calloc (n,n))
{
    setSingularvalue();
}


Eigenvalue::~Eigenvalue(){
    //gsl_matrix_free(Jacobianmatrix);
    //gsl_matrix_free(Hessiancartesian);
    //gsl_matrix_free(Massmatrix);
    gsl_matrix_free(Hessiantorsionangle);
    gsl_vector_complex_free(eigenvaluealpha);
    gsl_vector_free(singularvalue);
    gsl_vector_free(eigenvaluebeta);
    gsl_matrix_complex_free(eigenvector);
}

gsl_matrix* Eigenvalue::times(gsl_matrix* matrix1, gsl_matrix* matrix2) const{
    gsl_matrix* Hessiantorsionangle1 = gsl_matrix_calloc(n,m);
    gsl_matrix* Hessiantorsionangle2 = gsl_matrix_calloc(n,n);
    gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, matrix1, matrix2, 0.0, Hessiantorsionangle1);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Hessiantorsionangle1, matrix1, 0.0, Hessiantorsionangle2);

    gsl_matrix_free(Hessiantorsionangle1);

    return Hessiantorsionangle2;
}

gsl_matrix* Eigenvalue::getHessiantorsionangle() const{
    return Hessiantorsionangle;
}


gsl_vector* Eigenvalue::getSingularvalue() const{
    return singularvalue;
}

void Eigenvalue::setSingularvalue() const{

    gsl_eigen_genv_workspace *w = gsl_eigen_genv_alloc(n);

    gsl_matrix_memcpy(Hessiantorsionangle, times(Jacobianmatrix,Hessiancartesian));

    gsl_eigen_genv(Hessiantorsionangle, times(Jacobianmatrix,Massmatrix), eigenvaluealpha, eigenvaluebeta, eigenvector, w);

    gsl_eigen_genv_free(w);
    /*
    gsl_eigen_symmv_sort (eigenvalue, eigenvector, GSL_EIGEN_SORT_ABS_ASC);

     * gsl_eigen_gensymmv_workspace *w = gsl_eigen_gensymmv_alloc(n);

    gsl_eigen_gensymmv(Hessiantorsionangle, eig->getHessiantorsionangle(), eigenvalue, eigenvector, w);

    gsl_eigen_gensymmv_free(w);

    gsl_eigen_symmv_sort (eigenvalue, eigenvector, GSL_EIGEN_SORT_ABS_ASC);
     */

    for(int i=0; i<n; i++){
        double value=gsl_vector_get(eigenvaluebeta,i);
        if(value!=0){
            gsl_complex betavalue = gsl_complex_rect(value, 0.0);
            gsl_complex dividevalue = gsl_complex_div(gsl_vector_complex_get(eigenvaluealpha,i), betavalue);
            gsl_complex singular = gsl_complex_sqrt(dividevalue);
            if(GSL_REAL(singular)>1e-12){
                gsl_vector_set(singularvalue,i,GSL_REAL(singular));
            }
        }
    }
    gsl_sort_vector(singularvalue);

}

