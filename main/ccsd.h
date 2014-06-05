#include <psi4-dec.h>
#include <libmints/typedefs.h>

namespace psi{
    class Options;
    namespace main{

class CCSD{
  public:
    /// The constructor
    CCSD(Options &options, int nso_uhf, SharedVector Epsilon, double**** tei_MO, SharedMatrix Hcore_AO, SharedMatrix C_MO);
    /// The destructor
    ~CCSD();
    /// Compute the CCSD energy, and returns it
    double compute_energy();
  protected:
    /// Maximum iterations of CC optimization
    int maxiter_;
    /// Energy convergence for CC iterations
    double e_convergence_; 
    /// Number of spin orbitals
    int nso_;
    /// Number of occupied spin orbitals
    int nsocc_;
    /// Orbital energies
    SharedVector Epsilon_;
    /// Core Hamiltonian matrix (ASO)
    SharedMatrix Hcore_MO_;
    /// 2-e integrals in MSO (chemists' notation)
    double**** tei_MO_;
    /// Fock matrix (MSO)
//    SharedMatrix Fock_MO_;
    double** Fock_;
    /// A 2d array holds T1 amplitudes
    double** T1_;
    /// A 4d array holds T2 amplitudes
    double**** T2_;
    /// Effective double excitation operators
    double**** tau_tilde_;
    double**** tau_; 
    /// 2-index (F) intermediates
    double** F_;
    /// 4-index (W) intermediates
    double**** W_;
    /// CCSD energy
    double E_CCSD_;

    /// Transforms core Hamiltonian from ASO to MSO
    void transform_Hcore_AO2MO(SharedMatrix Hcore_AO, SharedMatrix C_MO);
    /// Form Fock matrix in MSO
    void form_Fock_MO();
    /// Build the 2-index (F) and 4-index (W) intermediates as well as effective doubles
    void calculate_intermediates();
    /// Update cluster amplitudes
    void update_amplitudes();
    /// Calculate CC energy using current T1 and T2 amplitudes
    void calculate_cc_energy();
    /// Kroneckder delta function
    double K_delta(int p, int q);
};

}} //End namespace
