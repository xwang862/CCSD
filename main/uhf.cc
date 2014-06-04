#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <uhf.h>
#include <ctime>

namespace psi{ namespace main{

UHF::UHF(Options &options):
  options_(options)
{
    // Start the timer
    const std::clock_t scf_start_init = std::clock();

    print_ = options_.get_int("PRINT");
    maxiter_ = options_.get_int("SCF_MAXITER");
    e_convergence_ = options_.get_double("E_CONVERGENCE");
    d_convergence_ = options_.get_double("D_CONVERGENCE");

    init_integrals(); // What is it??

    // End the timer
    const std::clock_t scf_end_init = std::clock();

    timer_init_ = double(scf_end_init - scf_start_init);
}


UHF::~UHF()
{
  free_matrix(tei_, nso_, nso_, nso_, nso_);
}

void UHF::init_matrix(double****& matrix, int size1, int size2, int size3, int size4){
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

// void UHF::free_matrix( double****& matrix, int size1, int size2, int size3, int size4){
//     if ((size1 == 0) || (size2 == 0) || (size3 == 0) || (size4 == 0)){
//         printf("\n\n\tNULL Matrix\n");
//     }
//     if (matrix != NULL){
//         for ( int i = 0; i < size1; i++){
//             for ( int j = 0; j < size2; j++){
//                 for ( int k = 0; k < size3; k++){
//                     delete[] matrix [i][j][k];
//                 }
//             }
//         }
//         for ( int i = 0; i < size1; i++){
//             for ( int j = 0; j < size2; j++){
//                 delete[] matrix[i][j];
//             }
//         }
//         for ( int i = 0; i < size1; i++){
//             delete[] matrix[i];
//         }
//         delete[] matrix;
//     }
// }

void UHF::free_matrix( double****& matrix, int size1, int size2, int size3, int size4){
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

void UHF::init_integrals()
{
    // This grabs the current molecule
    boost::shared_ptr<Molecule> molecule = Process::environment.molecule();

    // Form basis object:
    // Create a basis set parser object
    boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
    // Construct a new basis set
    boost::shared_ptr<BasisSet> aoBasis = BasisSet::construct(parser, molecule, "BASIS");

    // The integral factory oversees the creation of integrals objects
    boost::shared_ptr<IntegralFactory> integral(new IntegralFactory
        (aoBasis, aoBasis, aoBasis, aoBasis));

    // Determine the number of electrons in the system
    // N.B. This should be done after the basis has been built, because the geometry hasnot been 
    // fully initialized until this moment
    int charge = molecule->molecular_charge();
    int nelec = 0;
    for (int i = 0; i < molecule->natom(); ++i)
        nelec += (int)molecule->Z(i);
    nelec -= charge;
    if (nelec % 2)
        throw PSIEXCEPTION("This is only an RHF code, but you gave it an odd number of electrons. Try again!"); //What is it?
    ndocc_ = nelec / 2;

    nsocc_ = nelec;

    fprintf(outfile, "\tThere are %d doubly occupied orbitals\n", ndocc_);
    molecule->print();
    if (print_ > 1){
        aoBasis->print_detail();
        options_.print();
    }

    // Number of AO basis functions
    nao_ = aoBasis->nbf();
    // Number of Spin Orbitals (SO)
    nso_ = 2 * aoBasis->nbf();

    e_nuc_ = molecule->nuclear_repulsion_energy();
    fprintf(outfile, "\n   Nuclear repulsion energy: %20.10f\n\n", e_nuc_);

    // These don't need to be declared since they belong to the class
    S_ = SharedMatrix(new Matrix("Overlap matrix (SO)", nso_, nso_));
    H_ = SharedMatrix(new Matrix("Core Hamiltonian (SO)", nso_, nso_)); 
//    // These don't belong to the class, so we have to define them as having type SharedMatrix
//    SharedMatrix T = SharedMatrix(new Matrix("Kinetic integrals matrix (SO)", nso_, nso_));
//    SharedMatrix V = SharedMatrix(new Matrix("Potential integrals matrix (SO)", nso_, nso_));
    
    // Define corresponding integrals in AO basis    
    SharedMatrix S_AO = SharedMatrix(new Matrix("Overlap matrix (AO)", nao_, nao_));
    SharedMatrix H_AO = SharedMatrix(new Matrix("Core Hamiltonian (AO)", nao_, nao_));
    SharedMatrix T_AO = SharedMatrix(new Matrix("Kinetic integrals matrix (AO)", nao_, nao_));
    SharedMatrix V_AO = SharedMatrix(new Matrix("Potential integral matrix (AO)", nao_, nao_));

    // Form the one-electron integral objects from the integral factory
    boost::shared_ptr<OneBodyAOInt> sOBI(integral->ao_overlap());
    boost::shared_ptr<OneBodyAOInt> tOBI(integral->ao_kinetic());
    boost::shared_ptr<OneBodyAOInt> vOBI(integral->ao_potential());
    // Compute the one electron integrals, telling each object where to store the results
    sOBI->compute(S_AO);
    tOBI->compute(T_AO);
    vOBI->compute(V_AO);

    // Form h = T + V by first cloning T and the adding V
    H_AO->copy(T_AO);
    H_AO->add(V_AO);

    // Form H_ and S_ matrices in SO basis
    for (int i = 0; i < nao_; ++i){
        for (int j = 0; j < nao_; ++j){
            int I = i + nao_;
            int J = j + nao_;
            H_->set(0, i, j, H_AO->get(0, i, j));
            H_->set(0, I, J, H_AO->get(0, i, j));
            S_->set(0, i, j, S_AO->get(0, i, j));
            S_->set(0, I, J, S_AO->get(0, i, j));
        }
    }
     
    if (print_ > 3){
        S_->print();
        T_AO->print();
        V_AO->print();
        H_->print();
    }

    fprintf(outfile, "\tForming Two-electron Integrals\n\n");
    
    // Allocate some storage for the integrals
    init_matrix(tei_, nso_, nso_, nso_, nso_);

    // Now, the two-electron integrals
    boost::shared_ptr<TwoBodyAOInt> eri(integral->eri());
    // The buffer will hold the integrals for each shell, as they're computed
    const double *buffer = eri->buffer();
    // The iterator conveniently lets us iterate over functions within shells
    AOShellCombinationsIterator shellIter = integral->shells_iterator(); //What is it?

    int count = 0;
    for ( shellIter.first(); shellIter.is_done() == false; shellIter.next()){
        // Compute quartet
        eri->compute_shell(shellIter);
        // From the quartet get all the integrals
        AOIntegralsIterator intIter = shellIter.integrals_iterator();

        for (intIter.first(); intIter.is_done() == false; intIter.next()){
            int p = intIter.i();
            int q = intIter.j();
            int r = intIter.k();
            int s = intIter.l();
            int P = p + nao_;
            int Q = q + nao_;
            int R = r + nao_;
            int S = s + nao_;

            double val = buffer[intIter.index()];

            // Store ( p q | r s ): Chemists' Notation
            tei_[p][q][r][s] = tei_[p][q][s][r] = tei_[q][p][r][s] = tei_[q][p][s][r] =
              tei_[r][s][p][q] = tei_[s][r][p][q] = tei_[r][s][q][p] = tei_[s][r][q][p] =
              tei_[P][Q][R][S] = tei_[P][Q][S][R] = tei_[Q][P][R][S] = tei_[Q][P][S][R] =
              tei_[R][S][P][Q] = tei_[S][R][P][Q] = tei_[R][S][Q][P] = tei_[S][R][Q][P] =
              tei_[p][q][R][S] = tei_[p][q][S][R] = tei_[q][p][R][S] = tei_[q][p][S][R] =
              tei_[r][s][P][Q] = tei_[s][r][P][Q] = tei_[r][s][Q][P] = tei_[s][r][Q][P] =
              tei_[P][Q][r][s] = tei_[P][Q][s][r] = tei_[Q][P][r][s] = tei_[Q][P][s][r] =
              tei_[R][S][p][q] = tei_[S][R][p][q] = tei_[R][S][q][p] = tei_[S][R][q][p] = val;

            if ( print_ > 4)
                fprintf(outfile, "\t(%2d %2d | %2d %2d) = %20.15f\n", p, q, r, s, val);

            ++count;
        }
    }
    fprintf(outfile, "\n\tThere are %d unique integrals.\n\n", count);
}

double UHF::compute_electronic_energy()
{
    Matrix HplusF;
    HplusF.copy(H_);
    HplusF.add(F_);
    return 0.5 * D_->vector_dot(HplusF);
}

void UHF::form_density()
{
    for ( int p = 0; p < nso_; ++p){
      for ( int q = 0; q < nso_; ++q){
        double val = 0.0;
        for( int i = 0; i < nsocc_; ++i){
          val += C_->get(p, i) * C_->get(q, i);
        }
        D_->set(p, q, val);
      }
    }
}

double UHF::compute_energy()
{

  // Start the timer
  const std::clock_t scf_start_comput = std::clock();

  // Allocate some matrices
  X_ = SharedMatrix(new Matrix("S^-1/2", nso_, nso_));
  F_ = SharedMatrix(new Matrix("Fock matrix", nso_, nso_));
  Ft_ = SharedMatrix(new Matrix("Transformed Fock Matrix", nso_, nso_));
  C_ = SharedMatrix(new Matrix("MO Coefficients", nso_, nso_));
  D_ = SharedMatrix(new Matrix("The Density Matrix", nso_, nso_));
  SharedMatrix Temp1(new Matrix("Temporary Array 1", nso_, nso_));
  SharedMatrix Temp2(new Matrix("Temporary Array 2", nso_, nso_));
  SharedMatrix FDS(new Matrix("FDS", nso_, nso_));
  SharedMatrix SDF(new Matrix("SDF", nso_, nso_));
  SharedMatrix Evecs(new Matrix("Eigenvectors", nso_, nso_));

  // Eigenvalues don't need to be declared because it's in .h file
  Evals_ = SharedVector(new Vector("Eigenvalues", nso_));

  // Form the X_ Matrix (S^-1/2)
  S_->diagonalize(Evecs, Evals_);
  for (int p = 0; p < nso_; ++p){
    double val = 1.0 / sqrt(Evals_->get(p));
    Evals_->set(p, val);
  }
  Temp1->set_diagonal(Evals_);
  Temp2->gemm(false, true, 1.0, Temp1, Evecs, 0.0);
  X_->gemm(false, false, 1.0, Evecs, Temp2, 0.0);

  F_->copy(H_);
  Ft_->transform(F_, X_);//What is it?
  Ft_->diagonalize(Evecs, Evals_);

  C_->gemm(false, false, 1.0, X_, Evecs, 0.0);
  form_density();
  if(print_ > 1){
    fprintf(outfile, "MO Coefficients and density from Core Hamiltonian guess:\n");
    C_->print();
    D_->print();
  }

  int iter = 1;
  bool converged = false;
  double e_old;
  double e_new = e_nuc_ + compute_electronic_energy();

  fprintf(outfile, "\tEnergy from core Hamiltonian guess: %20.16f\n\n", e_new);

  fprintf(outfile, "\t*=======================================================*\n");
  fprintf(outfile, "\t* Iter       Energy            delta E    ||gradient||  *\n");
  fprintf(outfile, "\t*=======================================================*\n");

  while(!converged && iter < maxiter_){
    e_old = e_new;

    // Add the core Hamiltonian term to the Fock operator
    F_->copy(H_);
    // Add the two electron terms to the Fock operator
    for ( int p = 0; p < nso_; ++p){
      for ( int q = 0; q < nso_; ++q){
        double J = 0.0;
        double K = 0.0;
        for ( int r = 0; r < nso_; ++r){
          for ( int s = 0; s < nso_; ++s){
            J += tei_[p][q][r][s] * D_->get(r, s);
            K += tei_[p][r][q][s] * D_->get(r, s);
          }
        }
        F_->add(p, q, J - K);
      }
    }

    // Transform the Fock operator and diagonalize it
    Ft_->transform(F_, X_);
    Ft_->diagonalize(Evecs, Evals_);

    // Form the orbitals from the eigenvectors of the transformed Fock matrix
    C_->gemm(false, false, 1.0, X_, Evecs, 0.0);

    // Rebuild the density using the new orbitals
    form_density();

    // Compute the energy
    e_new = e_nuc_ + compute_electronic_energy();
    double dE = e_new - e_old;

    // Compute the orbital gradient, FDS - SDF
    Temp1->gemm(false, false, 1.0, D_, S_, 0.0);
    FDS->gemm(false, false, 1.0, F_, Temp1, 0.0);
    Temp1->gemm(false, false, 1.0, D_, F_, 0.0);
    SDF->gemm(false, false, 1.0, S_, Temp1, 0.0);
    Temp1->copy(FDS);
    Temp1->subtract(SDF);
    double dRMS = Temp1->rms();

    if(print_ > 1){
      Ft_->print();
      Evecs->print();
      C_->print();
      D_->print();
      FDS->print();
      SDF->print();
      Temp1->set_name("Orbital gradient");
      Temp1->print();
    }

    converged = (fabs(dE) < e_convergence_) && (dRMS < d_convergence_);

    fprintf(outfile, "\t* %3d %20.14f    %9.2e    %9.2e    *\n", iter, e_new, dE, dRMS);
    iter++;
  }

  fprintf(outfile, "\t*=======================================================*\n");

  if (!converged)
    throw PSIEXCEPTION("The SCF iterations did not converge.");

  Evals_->set_name("Orbital Energies");
  Evals_->print();


  // End the timer
  const std::clock_t scf_end_comput = std::clock();

  timer_comput_ = double(scf_end_comput - scf_start_comput);

  // Print the time consumed by this module
  fprintf(outfile, "\n\t================= Timings for UHF Module =================\n\n");
  fprintf(outfile, "\tMatrix Initialization:  ...  %15.3f sec \n", timer_init_/CLOCKS_PER_SEC);
  fprintf(outfile, "\tSCF Iteration:          ...  %15.3f sec \n\n", timer_comput_/CLOCKS_PER_SEC);
  fprintf(outfile, "\tTotal Time:             ...  %15.3f sec \n\n", (timer_init_ + timer_comput_)/CLOCKS_PER_SEC);

  return e_new;
}

SharedVector UHF::return_Ft_Evals()
{
    return Evals_;
}

SharedMatrix UHF::return_MO_C()
{
    return C_;
}

int UHF::return_nso()
{
    return nso_;
}

double**** UHF::return_tei()
{
    return tei_;
}

SharedMatrix UHF::return_Hcore()
{
    return H_;
}


}} // End namespaces

