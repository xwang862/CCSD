#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <ump2.h>
#include <ctime>

namespace psi{ namespace main{

UMP2::UMP2(Options &options, int nso_uhf, SharedVector& Ft_Evals, SharedMatrix& C_MO, double****& matrix):
  options_(options)  
{
    ndocc_ = options.get_int_array("DOCC")[0];

    nso_ = nso_uhf;

    init_matrix(tei_mo_, nso_, nso_, nso_, nso_);

    transform_AO2MO(C_MO, matrix);

    Epsilon_ = SharedVector(new Vector("Orbital Energies", nso_));

    for (int i = 0; i < nso_; ++i){
        double val = Ft_Evals->get(0, i);
        Epsilon_->set(0, i, val);
    }
}

UMP2::~UMP2()
{
  free_matrix(tei_mo_, nso_, nso_, nso_, nso_);
}

void UMP2::init_matrix(double****& matrix, int size1, int size2, int size3, int size4){
    if ((size1 == 0) || (size2 == 0) || (size3 == 0) || (size4 == 0)){
        printf("\n\n\tNULL Matrix\n");
        matrix = NULL;
    }
    else{
        matrix = new double***[size1];
        for ( int i = 0; i < size1; i++){
            matrix[i] = new double**[size2];
            for (int j = 0; j < size2; j++){
                matrix[i][j] = new double*[size3];
                for ( int k = 0; k < size3; k++){
                    matrix[i][j][k] = new double[size4];
                    for ( int l = 0; l < size4; l++){
                        matrix[i][j][k][l] = 0.0;
                    }
                }
            }
        }
    }
}

void UMP2::free_matrix( double****& matrix, int size1, int size2, int size3, int size4){
    if ((size1 == 0) || (size2 == 0) || (size3 == 0) || (size4 == 0)){
        printf("\n\n\tNULL Matrix\n");
    }
    if (matrix != NULL){
        for ( int i = 0; i < size1; i++){
            for ( int j = 0; j < size2; j++){
                for ( int k = 0; k < size3; k++){
                    delete[] matrix [i][j][k];
                }
                delete[] matrix[i][j];
            }
            delete[] matrix[i];
        }
        delete[] matrix;
    }
}

void UMP2::transform_AO2MO(SharedMatrix& C_MO, double****& matrix)
{
    // Start the timer
    const std::clock_t ump2_start_transform = std::clock();

    //Transform 2-e integrals from AO to MO basis:
    double**** tempr;
    init_matrix(tempr, nso_, nso_, nso_, nso_);

    // 1st quarter tranform, (pq|rs) -> (pq|rl)
    for ( int p = 0; p < nso_; ++p){
        for ( int q = 0; q < nso_; ++q){
            for ( int r = 0; r < nso_; ++r){
                for (int l = 0; l < nso_; ++l){
                    double pqrl = 0.0;
                    for ( int s = 0; s < nso_; ++s){
                        pqrl += C_MO->get(0, s, l) * matrix[p][q][r][s];
                    }
                    tempr[p][q][r][l] = pqrl;
                }
            }
        }
    }

    // 2nd quarter transform, (pq|rl) -> (pq|kl)  
    for (int p = 0; p < nso_; ++p){
        for (int q = 0; q < nso_; ++q){
            for (int k = 0; k < nso_; ++k){
                for (int l = 0; l < nso_; ++l){
                    double pqkl = 0.0;
                    for (int r = 0; r < nso_; ++r){
                        pqkl += C_MO->get(0, r, k) * tempr[p][q][r][l];
                    }
                    tei_mo_[p][q][k][l] = pqkl;
                }
            }
        }
    }

    // 3rd quarter transform, (pq|kl) -> (pj|kl)
    for (int p = 0; p < nso_; ++p){
        for ( int j = 0; j < nso_; ++j){
            for ( int k = 0; k < nso_; ++k){
                for ( int l = 0; l < nso_; ++l){
                    double pjkl = 0.0;
                    for ( int q = 0; q < nso_; ++q){
                        pjkl += C_MO->get(0, q, j) * tei_mo_[p][q][k][l];
                    }
                    tempr[p][j][k][l] = pjkl;
                }
            }
        }
    }

    // 4th quarter tranform, (pj|kl) -> (ij|kl)
    for ( int i = 0; i < nso_; ++i){
        for ( int j = 0; j < nso_; ++j){
            for ( int k = 0; k < nso_; ++k){
                for ( int l = 0; l < nso_; ++l){
                    double ijkl = 0.0;
                    for ( int p = 0; p < nso_; ++p){
                        ijkl += C_MO->get(0, p, i) * tempr[p][j][k][l];
                    }
                    tei_mo_[i][j][k][l] = ijkl;
                }
            }
        }
    }

    free_matrix(tempr, nso_, nso_, nso_, nso_);

    // End the timer
    const std::clock_t ump2_end_transform = std::clock();

    timer_AO2MO_ = double(ump2_end_transform - ump2_start_transform);
}

double UMP2::compute_MP2_energy()
{

    // Start the timer
    const std::clock_t ump2_start_comput = std::clock();

    double e_MP2 = 0.0;

    for (int i = 0; i < 2*ndocc_; ++i){
        for (int j = 0; j < 2*ndocc_; ++j){
            for (int a = 2*ndocc_; a < nso_; ++a){
                for (int b = 2*ndocc_; b < nso_; ++b){
                    e_MP2 += pow(tei_mo_[i][a][j][b] - tei_mo_[i][b][j][a], 2.0)/(Epsilon_->get(0, i) + Epsilon_->get(0, j) - Epsilon_->get(0, a) - Epsilon_->get(0, b));
                }
            }
        }
    }

    e_MP2 *= 0.25;

    // End the timer
    const std::clock_t ump2_end_comput = std::clock();
    
    timer_comput_ = double(ump2_end_comput - ump2_end_comput);

    // Print the time consumed by this module
    fprintf(outfile, "\n\t================= Timings for UMP2 Module =================\n\n");
    fprintf(outfile, "\tAO->MO 2-eletron integral transformation:  ...  %15.3f sec \n", timer_AO2MO_/CLOCKS_PER_SEC);
    fprintf(outfile, "\tUMP2 energy calculation:                   ...  %15.3f sec \n\n", timer_comput_/CLOCKS_PER_SEC);
    fprintf(outfile, "\tTotal Time:                                ...  %15.3f sec \n\n", (timer_AO2MO_ + timer_comput_)/CLOCKS_PER_SEC);

    return e_MP2;

}

double**** UMP2::return_tei_mo()
{
    return tei_mo_;
}

}} // End namespace
