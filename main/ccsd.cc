#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <ccsd.h>
#include <mymatrix.h>
#include <ctime>

namespace psi{ namespace main{

CCSD::CCSD(Options &options, int nso_uhf, SharedVector Epsilon, double**** tei_MO, SharedMatrix Hcore_AO, SharedMatrix C_MO)
{
    fprintf(outfile, "\n\t================= CCSD Module =================\n");
    // Transfer maximum CC iterations
    maxiter_ = options.get_int("CC_MAXITER");

    // Energy convergence for TCSCF iterations
    e_convergence_ = options.get_double("E_CONVERGENCE");

    // The number of occupied spin orbitals
    nsocc_ = 2 * options.get_int_array("DOCC")[0];

    // The number of spin orbitals
    nso_ = nso_uhf;

    // Orbital energies
    Epsilon_ = SharedVector(Epsilon);

    // 2-e integrals in MSO (chemists' notation)
    tei_MO_ = tei_MO;  

    // Initialize core Hamiltonian matrix (MSO)
    Hcore_MO_ = SharedMatrix(new Matrix("Core Hamiltonian (MSO)", nso_, nso_));
    // Transforms core Hamiltonian from ASO to MSO
    transform_Hcore_AO2MO(Hcore_AO, C_MO); 

    // Initialize Fock matrix (MSO)
    init_matrix2d(Fock_, nso_, nso_);
//    Fock_MO_ = SharedMatrix(new Matrix("Fock matrix (MSO)", nso_, nso_));
    // Form Fock matrix (MSO)
    form_Fock_MO();
    
    // Initialize T1 and T2 amplitudes
    init_matrix2d(T1_, nso_, nso_);
    init_matrix4d(T2_, nso_, nso_, nso_, nso_);

}

CCSD::~CCSD()
{

}

void CCSD::transform_Hcore_AO2MO(SharedMatrix Hcore_AO, SharedMatrix C_MO)
{
    for (int p = 0; p < nso_; p++){
        for (int q = 0; q < nso_; q++){
            double buff = 0.0;
            for (int u = 0; u < nso_; u++){
                for (int v = 0; v < nso_; v++){
                    buff += C_MO->get(0, u, p) * C_MO->get(0, v, q) * Hcore_AO->get(0, u, v);
                }
            }
            Hcore_MO_->set(0, p, q, buff); 
        }
    }   
}

void CCSD::form_Fock_MO()
{
    for (int p = 0; p < nso_; p++){
        for (int q = 0; q < nso_; q++){
            double buff = Hcore_MO_->get(0, p, q);
            for (int m = 0; m < nsocc_; m++){
                buff += tei_MO_[p][q][m][m] - tei_MO_[p][m][m][q];
            }
            Fock_[p][q] = buff;
            //Fock_MO_->set(0, p, q, buff);
        }
    }
}

void CCSD::build_intermediates()
{
    // Initialize the effective doubles and 2- and 4- index intermediates
    init_matrix4d(tau_, nso_, nso_, nso_, nso_);
    init_matrix4d(tau_tilde_, nso_, nso_, nso_, nso_);
    init_matrix2d(F_, nso_, nso_);
    init_matrix4d(W_, nso_, nso_, nso_, nso_);

    // Form effective doubles
    for (int i = 0; i < nsocc_; i++){
        for (int j = 0; j < nsocc_; j++){
            for (int a = nsocc_; a < nso_; a++){
                for (int b = nsocc_; b < nso_; b++){
                    tau_tilde_[i][j][a][b] = T2_[i][j][a][b] + 0.5 * (T1_[i][a] * T1_[j][b] - T1_[i][b] * T1_[j][a]);
                    tau_[i][j][a][b] = T2_[i][j][a][b] + T1_[i][a] * T1_[j][b] - T1_[i][b] * T1_[j][a];
                }
            }
        }
    } 

    // Form the 2-index intermediates: F_[a][e]
    for (int a = nsocc_; a < nso_; a++){
        for (int e = nsocc_; e < nso_; e++){
            // Some summations
            double sum1 = 0.0;
            for (int m = 0; m < nsocc_; m++){
                sum1 += Fock_[m][e] * T1_[m][a];
            }

            double sum2 = 0.0;
            for (int m = 0; m < nsocc_; m++){
                for (int f = nsocc_; f < nso_; f++){
                    sum2 += T1_[m][f] * (tei_MO_[m][f][a][e] - tei_MO_[m][e][a][f]);
                }
            }

            double sum3 = 0.0;
            for (int m = 0; m < nsocc_; m++){
                for (int n = 0; n < nsocc_; n++){
                    for (int f = nsocc_; f < nso_; f++){
                        sum3 += tau_tilde_[m][n][a][f] * (tei_MO_[m][e][n][f] - tei_MO_[m][f][n][e]);
                    }
                }
            }

            // According to Eq (3)
            F_[a][e] = (1.0 - K_delta(a, e)) * Fock_[a][e] - 0.5 * sum1 + sum2 - 0.5 * sum3;
        }
    } 

    // Form the 2-index intermediates: F_[m][i]
    for (int m = 0; m < nsocc_; m++){
        for (int i = 0; i < nsocc_; i++){
            double sum1 = 0.0;
            for (int e = nsocc_; e < nso_; e++){
                sum1 += T1_[i][e] * Fock_[m][e];
            }
            
            double sum2 = 0.0;
            for (int e = nsocc_; e < nso_; e++){
                for (int n = 0; n < nsocc_; n++){
                    sum2 += T1_[n][e] * (tei_MO_[m][i][n][e] - tei_MO_[m][e][n][i]); 
                }
            }

            double sum3 = 0.0;
            for (int n = 0; n < nsocc_; n++){
                for (int e = nsocc_; e < nso_; e++){
                    for (int f = nsocc_; f < nso_; f++){
                        sum3 += tau_tilde_[i][n][e][f] * (tei_MO_[m][e][n][f] - tei_MO_[m][f][n][e]);
                    }
                }
            }

            //According to Eq (4)
            F_[m][i] = (1.0 - K_delta(m, i)) * Fock_[m][i] + 0.5 * sum1 + sum2 + 0.5 * sum3;
        }
    }

    // Form the 2-index intermediates: F_[m][e]
    for (int m = 0; m < nsocc_; m++){
        for (int e = nsocc_; e < nso_; e++){
            double sum1 = 0.0;
            for (int n = 0; n < nsocc_; n++){
                for (int f = nsocc_; f < nso_; f++){
                    sum1 += T1_[n][f] * (tei_MO_[m][e][n][f] - tei_MO_[m][f][n][e]);
                }
            }
            //According to Eq (5)
            F_[m][e] = Fock_[m][e] + sum1;
        }
    }

    // Form the 4-index intermediates: W_[m][n][i][j]
}

double CCSD::K_delta(int p, int q)
{
    if (p == q) return 1.0;
    else return 0.0;
}

double CCSD::compute_energy()
{
    // Initial guess for T1( t[i][a] = 0.0 is already done ) and T2 amplitudes
    for (int i = 0; i < nsocc_; i++){
        for (int j = 0; j < nsocc_; j++){
            for (int a = nsocc_; a < nso_; a++){
                for (int b = nsocc_; b < nso_; b++){
                    T2_[i][j][a][b] = (tei_MO_[i][a][j][b] - tei_MO_[i][b][j][a]) / (Epsilon_->get(0, i) + Epsilon_->get(0, j) - Epsilon_->get(0, a) - Epsilon_->get(0, b));
                }
            }
        }
    }
    
    // Check for the initial-guess cluster amplitudes
    double E_mp2 = 0.0;
    for (int i = 0; i < nsocc_; i++){
        for (int j = 0; j < nsocc_; j++){
            for (int a = nsocc_; a < nso_; a++){
                for (int b = nsocc_; b < nso_; b++){
                    E_mp2 += (tei_MO_[i][a][j][b] - tei_MO_[i][b][j][a]) * T2_[i][j][a][b];
                }
            }
        }
    }
    E_mp2 *= 0.25;    
    fprintf(outfile, "\n\tCheck for initial-guess cluster amplitudes by mp2 energy: %20.15f\n", E_mp2);
    // End Check

    // Build the 2-index (F) and 4-index (W) intermediates as well as effective doubles
    build_intermediates();



    double E_ccsd = 0.0;

    return E_ccsd; 
}

}}//End of namespace
