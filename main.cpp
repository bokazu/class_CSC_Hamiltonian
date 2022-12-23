#include "CSC_Hamiltonian/CSC_Hamiltonian.hpp"
#include "CSC_Sub_Hamiltonian/CSC_Sub_Hamiltonian.hpp"

using namespace std;

int main()
{

    /*全状態空間での対角化*/
    // CSC_Hamiltonian H("model/jset_4site.txt", 4);
    // H.csc_hamiltonian();
    // H.csc_lanczos(10);
    // cout << H << endl;

    /*部分状態空間での対角化*/
    CSC_Sub_Hamiltonian H("model/jset_4site.txt", 4, 2);
    H.subspace_check();
    H.csc_sub_hamiltonian();
    // for (int i = 0; i < 6; i++)
    // {
    //     cout << "bm[" << i << "] = " << H.bm[i] << endl;
    //     cout << "gbm[bm[" << i << " ]] = " << H.gbm[H.bm[i]] << endl;
    // }

    // for (int i = 0; i < H.nnz; i++)
    // {
    //     cout << "row[" << i << "] = " << H.row[i] << endl;
    //     // cout << "val[" << i << "] = " << H.val[i] << endl;
    // }
    // for (int i = 0; i < H.sub_mat_dim + 1; i++)
    // {
    //     cout << "col_ptr[" << i << "] = " << H.col_ptr[i] << endl;
    // }
    H.csc_sub_lanczos(20, 'N', 'y');
    cout << H << endl;
}
