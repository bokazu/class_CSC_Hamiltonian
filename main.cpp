#include "CSC_Hamiltonian/CSC_Hamiltonian.hpp"
// #include "CSC_Sub_Hamiltonian/CSC_Sub_Hamiltonian.hpp"　//各磁化ごとの部分空間における対角化を行う場合はこちらを使用する

using namespace std;

int main()
{

    /*-------------------全状態空間での対角化------------------------*/
    CSC_Hamiltonian H("model/jset_4site.txt", 4); //オブジェクトの構築
    H.csc_hamiltonian(); //ハミルトニアン行列要素の計算
    H.csc_lanczos(10); //ランチョス法を用いた数値対角化の実行
    cout << H << endl; //結果の出力
    /*-----------------------------------------------------------*/

    /*------------各磁化ごとの部分空間における対角化を行う場合にはこちらを使用する------------*/
    // CSC_Sub_Hamiltonian H("model/jset12site_radder.txt", 12, 11);
    // H.subspace_check();
    // H.csc_sub_hamiltonian();
    // for (int i = 0; i < 6; i++)
    // {
    //     cout << "bm[" << i << "] = " << H.bm[i] << endl;
    //     cout << "gbm[bm[" << i << " ]] = " << H.gbm[H.bm[i]] << endl;
    // }
    /*-------------------------------------------------------------------------*/

    // for (int i = 0; i < H.nnz; i++)
    // {
    //     cout << "row[" << i << "] = " << H.row[i] << endl;
    //     // cout << "val[" << i << "] = " << H.val[i] << endl;
    // }
    // for (int i = 0; i < H.sub_mat_dim + 1; i++)
    // {
    //     cout << "col_ptr[" << i << "] = " << H.col_ptr[i] << endl;
    // }
    // H.csc_sub_lanczos(100, 'V');
    // cout << H << endl;
}
