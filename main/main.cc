#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>
#include <ctime>
#include <uhf.h>
#include <ump2.h>
#include <ccsd.h>

INIT_PLUGIN

using namespace boost;

namespace psi{ namespace main {

extern "C" 
int read_options(std::string name, Options& options)
{
    if (name == "MAIN"|| options.read_globals()) {
        /*- The amount of information printed to the output file -*/
        options.add_int("PRINT", 1);
        // Energy convergence
        options.add_double("E_CONVERGENCE", 1.0E-6);
        // Density convergence
        options.add_double("D_CONVERGENCE", 1.0E-8);
        // Maximum SCF iterations
        options.add_int("SCF_MAXITER", 300);
        // This array contains the number of doubly occupied MOs per irrep
        options.add_array("DOCC");
        // Maximum CC iterations
        options.add_int("CC_MAXITER", 300);
    }

    return true;
}

extern "C" 
PsiReturnType main(Options& options)
{
    int print = options.get_int("PRINT");

    /* Your code goes here */

    // Start the timer
    const std::clock_t wall_time_start = std::clock();

    // Construct a UHF object
    UHF uhf(options);

    // Run UHF code, and return UHF energy
    double E_uhf = uhf.compute_energy();

    // In order to build a UMP2 object, we need transfer some values as arguments
    // Transfer number of spin orbitals
    int nso_uhf = uhf.return_nso();

    // Form a vector to transfer HF orbital energies
    SharedVector Epsilon(uhf.return_Ft_Evals());

    // Form a matrix to transfer MO coefficients
    SharedMatrix C_MO(uhf.return_MO_C());

    // Form a 4d array to transfer 2-e integrals (ASO)
    double**** tei_AO = uhf.return_tei();

    // Construct a UMP2 object
    UMP2 ump2(options, nso_uhf, Epsilon, C_MO, tei_AO);

    // Run UMP2 code and return UMP2 correlation energy
    double E_mp2 = ump2.compute_MP2_energy();

    // In order to build a CCSD object, we need some values as arguments
    // Form a 4d array to transfer 2-e integrals (MSO)
    double**** tei_MO = ump2.return_tei_mo();

    // Form a matrix to transfer core Hamiltonian (ASO)
    SharedMatrix Hcore_AO(uhf.return_Hcore());

    // Construct a CCSD object
    CCSD ccsd(options, nso_uhf, Epsilon, tei_MO, Hcore_AO, C_MO);

    // Compute CCSD energy and return it
    double E_ccsd = ccsd.compute_energy(); 

    // Energy Summary 
    fprintf(outfile, "\n\n\t================= Summary =================\n");
    fprintf(outfile, "\n\tThe UHF Energy:           %20.16f a.u.\n", E_uhf);
    fprintf(outfile, "\n\tThe UMP2 Energy:          %20.16f a.u.\n", E_mp2);
    fprintf(outfile, "\n\tThe Total UMP2 Energy:    %20.16f a.u.\n", E_uhf + E_mp2);
    fprintf(outfile, "\n\tThe CCSD Energy:          %20.16f a.u.\n", E_ccsd);
    fprintf(outfile, "\n\tThe Total CCSD Energy:    %20.16f a.u.\n", E_uhf + E_ccsd);

    // End the timer
    const std::clock_t wall_time_end = std::clock();

    // Print the wall time
    fprintf(outfile, "\n\n\tWall Time:                %20.3f sec \n\n", double(wall_time_end - wall_time_start)/CLOCKS_PER_SEC);

    return Success;
}

}} // End namespaces

