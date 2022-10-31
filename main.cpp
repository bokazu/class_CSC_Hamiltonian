#include "CSC_Hamiltonian/CSC_Hamiltonian.hpp"

using namespace std;

int main()
{
    CSC_Hamiltonian H("jset_16site_pb.txt", 16);
    H.csc_hamiltonian();
    H.csc_lanczos(200);
    cout << H << endl;
}
