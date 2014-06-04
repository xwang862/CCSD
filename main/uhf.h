#include <psi4-dec.h>
#include <libmints/typedefs.h>

namespace psi{
    class Options;
    namespace main{

class UHF{
  public:
    /// The constructor
    UHF(Options &options);
    /// The destructor
    ~UHF();
    /// Computes the SCF energy, and returns it
    double compute_energy();
    /// Returns the number of spin orbitals
    int return_nso();
    /// Returns the transformed (canonical) Fock matrix
    SharedVector return_Ft_Evals();
    /// Returns the MO coefficient matrix
    SharedMatrix return_MO_C();
    /// Returns values of the 2-e integrals in AO
    double**** return_tei();
    /// Returns core Hamiltonian matrix in AO
    SharedMatrix return_Hcore();
  protected:
    /// The Options object, to interact with the input file
    Options &options_;
    /// The amount of information to print to the output file
    int print_;
    /// The number of doubly occupied orbitals
    int ndocc_;
    /// The number of singly occupied spin orbitals
    int nsocc_;
    /// The number of AO basis functions (delete if not used by member functions other than init_integral)
    int nao_;
    /// The number of symmetrized spin orbitals
    int nso_;
    /// The maximum number of SCF iterations
    int maxiter_;
    /// The nuclear repulsion energy
    double e_nuc_;
    /// The convergence criterion for the density
    double d_convergence_;
    /// The convergence criterion for the energy
    double e_convergence_;
    /// The one electron integrals
    SharedMatrix H_;
    /// The overlap matrix
    SharedMatrix S_;
    /// The inverse square root of the overlap matrix
    SharedMatrix X_;
    /// The Fock Matrix
    SharedMatrix F_;
    /// The transformed Fock Matrix
    SharedMatrix Ft_;
    /// The MO coefficients
    SharedMatrix C_;
    /// The density matrix
    SharedMatrix D_;
    /// A 4d array containing all two electron integrals
    double ****tei_;
    /// Eigenvalues of Fock matrix
    SharedVector Evals_;
    /// Computes the electronic part of the SCF energy, and returns it
    double compute_electronic_energy();
    /// Initializes the integrals object
    void init_integrals();
    /// Form the density matrix form the MO coefficients
    void form_density();
    /// Initializes a 4 dimensional array, setting the elements to zero
    void init_matrix(double****& matrix, int dim1, int dim2, int dim3, int dim4);
    /// Frees the memory used for a 4D array
    void free_matrix(double****& matrix, int dim1, int dim2, int dim3, int dim4);
    /// Timers for checking how effienct each part is
    double timer_init_;
    double timer_comput_;
};

}} //End namespaces
