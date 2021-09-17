
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
        Eigenvalue(gsl_matrix* matrix1, gsl_matrix* matrix2, gsl_matrix* matrix3);
        gsl_matrix * const Jacobianmatrix;  //TODO: Make private
        gsl_matrix * const Hessiancartesian;       //TODO: Make private
        gsl_matrix * const Hessiantorsionangle;       //TODO: Make private
        gsl_matrix * const Massmatrix;
        gsl_vector_complex * const eigenvaluealpha;
        gsl_vector * const eigenvaluebeta;

        gsl_matrix_complex * const eigenvector;
        gsl_matrix* times(gsl_matrix* matrix1, gsl_matrix* matrix2) const;

        virtual ~Eigenvalue();

        gsl_matrix* getHessiantorsionangle() const;

        void setSingularvalue() const;

        gsl_vector* getSingularvalue() const;

        gsl_vector * const singularvalue;
    };


#endif //KGS_EIGENVALUE_H
