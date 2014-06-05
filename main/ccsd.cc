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

    // Form Fock matrix (MSO)
    form_Fock_MO();
    
    // Initialize T1 and T2 amplitudes
    init_matrix2d(T1_, nso_, nso_);
    init_matrix4d(T2_, nso_, nso_, nso_, nso_);

    // Initialize the effective doubles and 2- and 4- index intermediates
    init_matrix4d(tau_, nso_, nso_, nso_, nso_);
    init_matrix4d(tau_tilde_, nso_, nso_, nso_, nso_);
    init_matrix2d(F_, nso_, nso_);
    init_matrix4d(W_, nso_, nso_, nso_, nso_);

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

void CCSD::calculate_intermediates()
{
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

            //Check for summations
            //if (a == 10 && e == 10){
            //    fprintf(outfile, "\n\tpart1 = %20.15f\n\tsum1 = %20.15f\n\tsum2 = %20.15f\n\tsum3 = %20.15f\n", (1.0 - K_delta(a, e)) * Fock_[a][e], sum1, sum2, sum3);
            //}
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
    for (int m = 0; m < nsocc_; m++){
        for (int n = 0; n < nsocc_; n++){
            for (int i = 0; i < nsocc_; i++){
                for (int j = 0; j < nsocc_; j++){
                    double sum11 = 0.0, sum12 = 0.0;
                    for (int e = nsocc_; e < nso_; e++){
                        sum11 += T1_[j][e] * (tei_MO_[m][i][n][e] - tei_MO_[m][e][n][i]);
                        sum12 += T1_[i][e] * (tei_MO_[m][j][n][e] - tei_MO_[m][e][n][j]);
                    }

                    double sum2 = 0.0;
                    for (int e = nsocc_; e < nso_; e++){
                        for (int f = nsocc_; f < nso_; f++){
                            sum2 += tau_[i][j][e][f] * (tei_MO_[m][e][n][f] - tei_MO_[m][f][n][e]);
                        }
                    }

                    //According to Eq (6)
                    W_[m][n][i][j] = (tei_MO_[m][i][n][j] - tei_MO_[m][j][n][i]) + (1.0 - (sum11 - sum12)) + 0.25 * sum2;
                }
            }
        }
    }

    // Form the 4-index intermediates: W_[a][b][e][f]
    for (int a = nsocc_; a < nso_; a++){
        for (int b = nsocc_; b < nso_; b++){
            for (int e = nsocc_; e < nso_; e++){
                for (int f = nsocc_; f < nso_; f++){
                    double sum11 = 0.0, sum12 = 0.0;
                    for (int m = 0; m < nsocc_; m++){
                        sum11 += T1_[m][b] * (tei_MO_[a][e][m][f] - tei_MO_[a][f][m][e]);
                        sum12 += T1_[m][a] * (tei_MO_[b][e][m][f] - tei_MO_[b][f][m][e]);
                    }

                    double sum2 = 0.0;
                    for (int m = 0; m < nsocc_; m++){
                        for (int n = 0; n < nsocc_; n++){
                            sum2 += tau_[m][n][a][b] * (tei_MO_[m][e][n][f] - tei_MO_[m][f][n][e]);
                        }
                    }

                    //According to Eq (7)
                    W_[a][b][e][f] = (tei_MO_[a][e][b][f] - tei_MO_[a][f][b][e]) - (1.0 - (sum11 - sum12)) + 0.25 * sum2;
                }
            }
        }
    }

    // Form the 4-index intermediates: W_[m][b][e][j]
    for (int m = 0; m < nsocc_; m++){
        for (int b = nsocc_; b < nso_; b++){
            for (int e = nsocc_; e < nso_; e++){
                for (int j = 0; j < nsocc_; j++){
                    double sum1 = 0.0;
                    for (int f = nsocc_; f < nso_; f++){
                        sum1 += T1_[j][f] * (tei_MO_[m][e][b][f] - tei_MO_[m][f][b][e]);
                    }

                    double sum2 = 0.0;
                    for (int n = 0; n < nsocc_; n++){
                        sum2 += T1_[n][b] * (tei_MO_[m][e][n][j] - tei_MO_[m][j][n][e]);
                    }

                    double sum3 = 0.0;
                    for (int n = 0; n < nsocc_; n++){
                        for (int f = nsocc_; f < nso_; f++){
                            sum3 += (0.5 * T2_[j][n][f][b] + T1_[j][f] * T1_[n][b]) * (tei_MO_[m][e][n][f] - tei_MO_[m][f][n][e]);
                        }
                    }

                    //According to Eq (8)
                    W_[m][b][e][j] = (tei_MO_[m][e][b][j] - tei_MO_[m][j][b][e]) + sum1 - sum2 - sum3;
                }
            }
        }
    }
}

void CCSD::update_amplitudes()
{
    // Update T1 amplitudes
    for (int i = 0; i < nsocc_; i++){
        for (int a = nsocc_; a < nso_; a++){
            double sum1 = 0.0;
            for (int e = nsocc_; e < nso_; e++){
                sum1 += T1_[i][e] * F_[a][e];
            }

            double sum2 = 0.0;
            for (int m = 0; m < nsocc_; m++){
                sum2 += T1_[m][a] * F_[m][i];
            }

            double sum3 = 0.0;
            for (int m = 0; m < nsocc_; m++){
                for (int e = nsocc_; e < nso_; e++){
                    sum3 += T2_[i][m][a][e] * F_[m][e];
                }
            }

            double sum4 = 0.0;
            for (int n = 0; n < nsocc_; n++){
                for (int f = nsocc_; f < nso_; f++){
                    sum4 += T1_[n][f] * (tei_MO_[n][i][a][f] - tei_MO_[n][f][a][i]);
                }
            }

            double sum5 = 0.0;
            for (int m = 0; m < nsocc_; m++){
                for (int e = nsocc_; e < nso_; e++){
                    for (int f = nsocc_; f < nso_; f++){
                        sum5 += T2_[i][m][e][f] * (tei_MO_[m][e][a][f] - tei_MO_[m][f][a][e]);
                    }
                }
            }

            double sum6 = 0.0;
            for (int m = 0; m < nsocc_; m++){
                for (int e = nsocc_; e < nso_; e++){
                    for (int n  = 0; n < nsocc_; n++){
                        sum6 += T2_[m][n][a][e] * (tei_MO_[n][e][m][i] - tei_MO_[n][i][m][e]);
                    }
                }
            }

            //According to Eq (1)
            T1_[i][a] = (Fock_[i][a] + sum1 - sum2 + sum3 - sum4 - 0.5 * sum5 - 0.5 * sum6) / (Epsilon_->get(0, i) - Epsilon_->get(0, a));
        }
    }
    
    // Update T2 amplitudes
    for (int i = 0; i < nsocc_; i++){
        for (int j = 0; j < nsocc_; j++){
            for (int a = nsocc_; a < nso_; a++){
                for (int b = nsocc_; b < nso_; b++){
                    // Since T2 equation is long, let's break it into 5 lines as in the paper*
                    // *** Line 1 ***
                    double sum11 = 0.0, sum12 = 0.0;
                    for (int e = nsocc_; e < nso_; e++){
                        double sum11_inner = 0.0, sum12_inner = 0.0;
                        for (int m = 0; m < nsocc_; m++){
                            sum11_inner += T1_[m][b] * F_[m][e];
                            sum12_inner += T1_[m][a] * F_[m][e];
                        }
                        sum11 += T2_[i][j][a][e] * (F_[b][e] - 0.5 * sum11_inner);
                        sum12 += T2_[i][j][b][e] * (F_[a][e] - 0.5 * sum12_inner);
                    }
                    double line1 = (tei_MO_[i][a][j][b] - tei_MO_[i][b][j][a]) + 1.0 - (sum11 - sum12);

                    //*** Line 2 ***
                    double sum21 = 0.0, sum22 = 0.0;
                    for (int m = 0; m < nsocc_; m++){
                        double sum21_inner = 0.0, sum22_inner = 0.0;
                        for (int e = nsocc_; e < nso_; e++){
                            sum21_inner += T1_[j][e] * F_[m][e];
                            sum22_inner += T1_[i][e] * F_[m][e];
                        }
                        sum21 += T2_[i][m][a][b] * (F_[m][j] + 0.5 * sum21_inner);
                        sum22 += T2_[j][m][a][b] * (F_[m][i] + 0.5 * sum22_inner);
                    }
                    double line2 = 0.0 - (1.0 - (sum21 - sum22));

                    // *** Line 3 ***
                    double sum3 = 0.0;
                    for (int m = 0; m < nsocc_; m++){
                        for (int n = 0; n < nsocc_; n++){
                            sum3 += tau_[m][n][a][b] * W_[m][n][i][j];
                        }
                    }
                    double sum4 = 0.0;
                    for (int e = nsocc_; e < nso_; e++){
                        for (int f = nsocc_; f < nso_; f++){
                            sum4 += tau_[i][j][e][f] * W_[a][b][e][f];
                        }
                    }
                    double line3 = 0.5 * sum3 + 0.5 * sum4;
                
                    // *** Line 4 ***
                    double sum51 = 0.0, sum52 = 0.0, sum53 = 0.0, sum54 = 0.0;
                    for (int m = 0; m < nsocc_; m++){
                        for (int e = nsocc_; e < nso_; e++){
                            sum51 += T2_[i][m][a][e] * W_[m][b][e][j] - T1_[i][e] * T1_[m][a] * (tei_MO_[m][e][b][j] - tei_MO_[m][j][b][e]);
                            sum52 += T2_[i][m][b][e] * W_[m][a][e][j] - T1_[i][e] * T1_[m][b] * (tei_MO_[m][e][a][j] - tei_MO_[m][j][a][e]);
                            sum53 += T2_[j][m][a][e] * W_[m][b][e][i] - T1_[j][e] * T1_[m][a] * (tei_MO_[m][e][b][i] - tei_MO_[m][i][b][e]);
                            sum54 += T2_[j][m][b][e] * W_[m][a][e][i] - T1_[j][e] * T1_[m][b] * (tei_MO_[m][e][a][i] - tei_MO_[m][i][a][e]); 
                        }
                    }
                    double line4 = 1 + sum51 - sum52 - sum53 + sum54;

                    // *** Line 5 ***
                    double sum61 = 0.0, sum62 = 0.0;
                    for (int e = nsocc_; e < nso_; e++){
                        sum61 += T1_[i][e] * (tei_MO_[a][e][b][j] - tei_MO_[a][j][b][e]);
                        sum62 += T1_[j][e] * (tei_MO_[a][e][b][i] - tei_MO_[a][i][b][e]);
                    }
                    double sum71 = 0.0, sum72 = 0.0;
                    for (int m = 0; m < nsocc_; m++){
                        sum71 += T1_[m][a] * (tei_MO_[m][i][b][j] - tei_MO_[m][j][b][i]);
                        sum72 += T1_[m][b] * (tei_MO_[m][i][a][j] - tei_MO_[m][j][a][i]);
                    }
                    double line5 = 1 - (sum61 - sum62) - (1 - (sum71 - sum72));

                    // According to Eq (2)
                    double denom = Epsilon_->get(0, i) + Epsilon_->get(0, j) - Epsilon_->get(0, a) - Epsilon_->get(0, b);
                    T2_[i][j][a][b] = (line1 + line2 + line3 + line4 + line5) / denom;
                }
            }
        }
    }
}

void CCSD::calculate_cc_energy()
{
    double sum1 = 0.0;
    for (int i = 0; i < nsocc_; i++){
        for (int a = nsocc_; a < nso_; a++){
            sum1 += Fock_[i][a] * T1_[i][a];
        }
    }

    double sum2 = 0.0;
    for (int i = 0; i < nsocc_; i++){
        for (int j = 0; j < nsocc_; j++){
            for (int a = nsocc_; a < nso_; a++){
                for (int b = nsocc_; b < nso_; b++){
                    sum2 += (tei_MO_[i][a][j][b] - tei_MO_[i][b][j][a]) * T2_[i][j][a][b];
                }
            }
        }
    }

    double sum3 = 0.0;
    for (int i = 0; i < nsocc_; i++){
        for (int j = 0; j < nsocc_; j++){
            for (int a = nsocc_; a < nso_; a++){
                for (int b = nsocc_; b < nso_; b++){
                    sum3 += (tei_MO_[i][a][j][b] - tei_MO_[i][b][j][a]) * T1_[i][a] * T1_[j][b];
                }
            }
        }
    }

    E_CCSD_ = sum1 + 0.25 * sum2 + 0.5 * sum3;
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

    // Calculate the 2-index (F) and 4-index (W) intermediates as well as effective doubles
    calculate_intermediates();

    // Update cluster amplitudes
    update_amplitudes();

    // Calculate the CC energy
    calculate_cc_energy();

    // Setup for the CC iterations
    int iter = 1;
    double E_old = 0.0, E_new = E_CCSD_;    
    bool converged = false;

    fprintf(outfile, "\n\tCC Iterations Start: \n\n");
    fprintf(outfile, "\t*===========================================*\n");
    fprintf(outfile, "\t* Iter       CC Energy            delta E   *\n");
    fprintf(outfile, "\t*===========================================*\n");
    fprintf(outfile, "\t*   0 %20.14f                  *\n", E_new);

    while( !converged && iter < maxiter_) {
        E_old = E_new;
        calculate_intermediates();
        update_amplitudes();
        calculate_cc_energy();
        E_new = E_CCSD_;

        // Check for amplitudes and intermediates
        //if (iter < 4){
        //    print_matrix2d(T1_, nso_, nso_);
        //    print_matrix2d(F_, nso_, nso_);
        //}

        // Energy convergence
        double dE = E_new - E_old;
        converged = (fabs(dE) < e_convergence_);
        fprintf(outfile, "\t* %3d %20.14f   %9.2e      *\n", iter, E_new, dE);

        iter++;
    }

    fprintf(outfile, "\t*===========================================*\n");

    return E_new; 
}

}}//End of namespace
