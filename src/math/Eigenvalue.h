
#ifndef KGS_EIGENVALUE_H
#define KGS_EIGENVALUE_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>

class Eigenvalue {
    protected:
        const int m, n; ///< Dimensions of matrix


    public:
        Eigenvalue(gsl_matrix* matrix1, gsl_matrix* matrix2);
        gsl_matrix * const Jacobianmatrix;  //TODO: Make private
        gsl_matrix * const Hessiancartesian;       //TODO: Make private
        gsl_matrix * const Hessiantorsionangle;       //TODO: Make private
        gsl_vector * const eigenvalue;
        gsl_matrix * const eigenvector;

        virtual ~Eigenvalue();

        gsl_matrix* getHessiantorsionangle() const;

        gsl_vector* getEigenvalue() const;

    };


#endif //KGS_EIGENVALUE_H
