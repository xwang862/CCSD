#include <libmints/mints.h>
#include <mymatrix.h>

using namespace boost;

namespace psi{ namespace main{

    void init_matrix2d(double**& matrix, int size1, int size2){
        if (!size1 || !size2){
            printf("\n\n\tNULL Matrix\n");
            matrix = NULL;
        }
        else{
            matrix = new double* [size1];
            for (int i = 0; i < size1; i++){
                matrix[i] = new double [size2];
                for (int j = 0; j < size2; j++){
                    matrix[i][j] = 0.0;
                }
            }
        } 
    }

    void free_matrix2d(double**& matrix, int size1, int size2){
        if (!size1 || !size2){
            printf("\n\n\tNULL Matrix\n"); 
        }
        else{
            for (int i = 0; i < size1; i++){
                delete[] matrix[i];
            }
            delete[] matrix;
        }
    }

    void init_matrix4d(double****& matrix, int size1, int size2, int size3, int size4)
    {
        if ( !size1 || !size2 || !size3 || !size4){
            printf("\n\n\tNULL Matrix\n");
            matrix = NULL;
        }
        else{
            matrix = new double*** [size1];
            for (int i = 0; i < size1; i++){
                matrix[i] = new double** [size2];
                for (int j = 0; j < size2; j++){
                    matrix[i][j] = new double* [size3];
                    for (int k = 0; k < size3; k++){
                        matrix[i][j][k] = new double [size4];
                        for (int l = 0; l < size4; l++){
                            matrix[i][j][k][l] = 0.0;
                        }
                    }
                }
            }
        }
    }

    void free_matrix4d(double****& matrix, int size1, int size2, int size3, int size4)
    {
        if (!size1 || !size2 || !size3 || !size4){
            printf("\n\n\tNULL Matrix\n");
        }
        else{
            for (int i = 0; i < size1; i++){
                for (int j = 0; j < size2; j++){
                    for (int k = 0; k < size3; k++){
                        delete[] matrix[i][j][k];
                    }
                    delete[] matrix[i][j];
                }
                delete[] matrix[i];
            }
            delete[] matrix;
        }
    }

}}//End namespace
