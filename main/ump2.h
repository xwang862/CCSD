#include <psi4-dec.h>
#include <libmints/typedefs.h>

namespace psi{
    class Options;
    namespace main{

class UMP2{
  public:
    /// The constructor
    UMP2(Options &options, int nso_uhf, SharedVector& Ft_Evals, SharedMatrix& C_MO, double****& matrix);
    /// The destructor
    ~UMP2();
    /// Computes the MP2 energy, and returns it 
    double compute_MP2_energy();
    /// Return the 2-e integrals in MO
    double**** return_tei_mo();
  protected:
    /// The Options object, to interact with the input file
    Options &options_;
    /// The number of doubly occupied orbitals
    int ndocc_;
    /// The number of spin orbitals
    int nso_;
    /// A 4d array containing all two electron integrals in MO basis using spin orbitals (SO)
    double ****tei_mo_;
    /// Eigenvalues of Fock matrix (Orbital Energies)
    SharedVector Epsilon_;
    /// 2-electron integral transformation from AO to MO basis
    void transform_AO2MO(SharedMatrix& C_MO, double****& matrix);
    /// Initialize a 4d array, setting elements to zero
    void init_matrix(double****& matrix, int dim1, int dim2, int dim3, int dim4);
    /// Free the memory used for a 4d array
    void free_matrix(double****& matrix, int dim1, int dim2, int dim3, int dim4);
    /// Timers for checking how effienct each part is
    double timer_AO2MO_;
    double timer_comput_;

};

}} //End namespaces
