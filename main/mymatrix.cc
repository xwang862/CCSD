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

    void print_matrix2d(double**& matrix, int size1, int size2){
        if (!size1 || !size2){
            printf("\n\n\tNULL Matrix\n"); 
        }
        else{
            fprintf(outfile, "\n\t## Matrix size: %3d x %3d\n", size1, size2);
            for (int block = 0; block <= int(size2 / 5); block++){
                fprintf(outfile, "\n\t    ");
                for (int c = 0; c < 5; c++){
                    if ((block * 5 + c) < size2){
                        fprintf(outfile, "  %20d", block * 5 + c + 1);
                    }
                }
                fprintf(outfile, "\n\n");

                for (int r = 0; r < size1; r++){
                    fprintf(outfile, "\t%4d", r+1);
                    for (int c = 0; c < 5; c++){
                        if ((block * 5 + c) < size2){
                            fprintf(outfile, "  %20.10f", matrix[r][block * 5 + c]);
                        }
                    }
                    fprintf(outfile, "\n");
                }
                fprintf(outfile, "\n");
            }
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
