#include "CSC_Hamiltonian/CSC_Hamiltonian.hpp"

using namespace std;

int main()
{
    CSC_Hamiltonian H("jset.txt", 4);
    H.csc_hamiltonian();
    H.csc_lanczos(100, 'V');
    cout << H << endl;
}
