#include "gmx_CNC-HOH.h"

/* Proposed workflow : 
    1. VMD Stamp Structural alignment of reference structure to 1LFD Chain B
    2. Python script to use Kabsch algorithm to align GTPase of reference
        structure to gro structure
    3. trjconv to align trajectory to the aligned gro structure
    4. Make a .tpr file for 1LFD reference structure (NOT 1LFD.pdb from pdb.org
        due to stamp structural alignment moving it!)
    5. Make an ndx file with the Ral CA to align for each part
    6. Run this!
*/

int main(int argc, char * argv[])
{
    gmx_CNCHOH(argc, argv);
    return 0;
}
