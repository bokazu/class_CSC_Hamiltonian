#include "CSC_Sub_Hamiltonian.hpp"

#include <mkl.h>

#include <boost/dynamic_bitset.hpp>
#include <fstream>
#include <iomanip>
#include <random>

using namespace std;

// row, col_ptr,valの初期化を行う
void CSC_Sub_Hamiltonian::init()
{
    for (int i = 0; i < nnz; i++)
    {
        row[i] = 0;
        val[i] = 0.0;
    }

    for (int i = 0; i < sub_mat_dim + 1; i++)
    {
        col_ptr[i] = 0;
    }
}

// CSC形式の配列のrow,valの要素数を変更し、すべて0で初期化する
void CSC_Sub_Hamiltonian::set_sub_CscElem(int i)
{
    delete[] row;
    delete[] val;

    row = new int[i];
    val = new double[i];

    this->init();
}

// 二項係数の計算
int CSC_Sub_Hamiltonian::comb(int n, int r)
{
    vector<vector<int>> v(n + 1, vector<int>(n + 1, 0));
    for (int i = 0; i < v.size(); i++)
    {
        v[i][0] = 1;
        v[i][i] = 1;
    }

    for (int k = 2; k < v.size(); k++)
    {
        for (int j = 1; j < k; j++)
        {
            v[k][j] = (v[k - 1][j - 1] + v[k - 1][j]);
        }
    }
    return v[n][r];
}

// 部分空間を張るスピン状態を確認する
void CSC_Sub_Hamiltonian::subspace_check()
{
    delete[] gbm;

    int itr = 0;
    int gbm_size;

    // gbmのサイズ確認と配列のメモリ確保
    for (int m = mat_dim - 1; m >= 0; m--)
    {
        boost::dynamic_bitset<> ket_m(tot_site_num, m);
        bool flag = false; // このflagはfor loopの外に出さないと多分意味がない
        if (ket_m.count() == bit)
        {
            if (itr == 0)
            {
                gbm_size = m + 1;
                gbm = new int[gbm_size];
                vec_sub_init(gbm_size, gbm);
                flag = true;
            }
        }
        if (flag == true)
            break;
    }

    // bm、gbmへの要素の格納
    for (int m = 0; m < gbm_size; m++)
    {
        boost::dynamic_bitset<> ket_m(tot_site_num, m);
        if (ket_m.count() == bit)
        {
            bm[itr] = m;
            gbm[m] = itr;
            itr++;
        }
    }

    if (bm[sub_mat_dim - 1] != gbm_size - 1)
    {
        cout << "error::bm & gbm bad array. " << endl;
    }
}

// 部分空間におけるHamiltonian行列の非ゼロ要素をカウントする
void CSC_Sub_Hamiltonian::count_sub_nnz()
{
    nnz = 0;
    int m;
    int Jset_line = J.get_line();
    double szz;

    for (int bm_itr = 0; bm_itr < sub_mat_dim; bm_itr++)
    {
        m = bm[bm_itr];
        szz = 0.;
        for (int site_i = 0; site_i < tot_site_num; site_i++)
        {
            boost::dynamic_bitset<> ket_m(tot_site_num, m);
            bool bit_check0, bit_check1;

            int j_line = 0;
            for (int l = j_line; l < Jset_line; l++)
            {
                if (J.index(0, l) == site_i)
                {
                    int site_j = J.index(1, l);
                    bit_check0 = ket_m.test(site_i);
                    bit_check1 = ket_m.test(site_j);

                    // 注目する2サイトのスピン状態についての場合分け
                    if (bit_check0 == bit_check1)
                    {
                        szz += 0.25 * J.val(l);
                    }
                    else
                    {
                        nnz++;
                        szz -= 0.25 * J.val(l);
                    }
                }
            }
        }
        if (szz != 0)
        {
            nnz++;
        }
    }
}

// 部分空間におけるspin演算子と状態ベクトルの計算を行う
void CSC_Sub_Hamiltonian::csc_sub_spin(int bm_itr, int site_i, int &row_index, int &col_ptr_val, double &szz)
{
    int m = bm[bm_itr];
    boost::dynamic_bitset<> ket_m(tot_site_num, m);
    bool bit_check0, bit_check1;

    int j_line = 0;
    int Jset_line = J.get_line();

    for (int l = j_line; l < Jset_line; l++)
    {
        if (J.index(0, l) == site_i)
        {
            int site_j = J.index(1, l);
            bit_check0 = ket_m.test(site_i);
            bit_check1 = ket_m.test(site_j);

            // 注目する2サイトのスピン状態を場合分けする
            if (bit_check0 == bit_check1)
            {
                szz += 0.25 * J.val(l);
            }
            else
            {
                boost::dynamic_bitset<> ket_m1(tot_site_num, m);

                ket_m1.flip(site_i);
                ket_m1.flip(site_j);

                int n = (int)(ket_m1.to_ulong());
                int bm_ctr = gbm[n]; // spin状態nが部分空間における何番目の状態かを取得

                row[row_index] = bm_ctr;
                val[row_index] = 0.5 * J.val(l);
                row_index++;
                col_ptr_val++;

                szz -= 0.25 * J.val(l);
            }
        }
    }
}

// 部分空間におけるHamiltonian行列をCSC形式で格納する
void CSC_Sub_Hamiltonian::csc_sub_hamiltonian()
{
    // 特定の磁化における部分空間での非ゼロ要素数の算出
    count_sub_nnz();
    // CSC形式の配列row, valのメモリの再確保
    set_sub_CscElem(nnz);

    int row_index = 0;
    int col_ptr_index = 0;
    int col_ptr_val = 0;

    double szz;
    for (int bm_itr = 0; bm_itr < sub_mat_dim; bm_itr++)
    {
        szz = 0.;
        col_ptr[col_ptr_index] = col_ptr_val;
        col_ptr_index++;

        for (int site_i = 0; site_i < tot_site_num; site_i++)
        {
            csc_sub_spin(bm_itr, site_i, row_index, col_ptr_val, szz);
        }

        if (szz != 0.0)
        {
            row[row_index] = bm_itr;
            val[row_index] = szz;
            row_index++;
            col_ptr_val++;
        }
    }

    col_ptr[sub_mat_dim] = nnz;
}

/*-------------Lanczos algorithm関係の関数--------------------*/
// ベクトルの規格化用関数
void sub_sdz(int dim, double *vec)
{
    double a = 1. / cblas_dnrm2(dim, vec, 1);
    cblas_dscal(dim, a, vec, 1);
}

/*部分空間におけるCSC形式でのMV積の計算を行う*/
void CSC_Sub_Hamiltonian::csc_sub_mvprod(double *u_i, double *u_j)
{
    for (int i = 0; i < sub_mat_dim; i++)
    {
        for (int j = col_ptr[i]; j < col_ptr[i + 1]; j++)
        {
            u_j[i] += val[j] * u_i[row[j]];
        }
    }
}

void CSC_Sub_Hamiltonian::csc_sub_lanczos(int tri_mat_dim, char c, char info_ls)
{
    ls_count = 0;
    double eps = 1.0;        // step間の誤差を代入するための変数
    double err = 1.0e-15;    // 要求精度
    bool err_checker = true; // 誤差が要求精度の範囲内に収まっているかを確認するためのflag

    // 固有ベクトルを格納するための配列の要素数を変更する
    if (c == 'V')
    {
        Eig.evec_elem(sub_mat_dim);
    }

    /*-----------------初期ベクトルの用意------------------*/
    double **u = new double *[2];
    for (int k = 0; k < 2; k++)
    {
        u[k] = new double[sub_mat_dim];
    }

    random_device rand;
    mt19937 mt(rand());
    uniform_real_distribution<> rand1(0, 1);
    for (int k = 0; k < sub_mat_dim; k++)
    {
        u[0][k] = rand1(mt);
        u[1][k] = 0.0;
    }
    sub_sdz(sub_mat_dim, u[0]);

    if (c == 'V')
    {
        cblas_dcopy(sub_mat_dim, u[0], 1, Eig.data(), 1);
    }
    /*-----------------------------------------------------*/

    // 三重対角行列の主対角成分
    double *alpha = new double[tri_mat_dim];
    vec_sub_init(tri_mat_dim, alpha);

    // 三重対角行列の次対角成分
    double *beta = new double[tri_mat_dim - 1];
    vec_sub_init(tri_mat_dim - 1, beta);

    // ls = 偶数stepでの近似固有値
    double *eval_even = new double[tri_mat_dim];
    vec_sub_init(tri_mat_dim, eval_even);

    // ls = 奇数stepでの近似固有値
    double *eval_odd = new double[tri_mat_dim];
    vec_sub_init(tri_mat_dim, eval_odd);

    // LAPACKに三重対角行列の主対角成分を渡す用の配列
    double *diag = new double[tri_mat_dim];
    vec_sub_init(tri_mat_dim, diag);

    // LAPACKに三重対角行列の主対角成分を渡す用の配列
    double *sub_diag = new double[tri_mat_dim - 1];
    vec_sub_init(tri_mat_dim - 1, sub_diag);

    // LAPACKに渡し、c = 'N'なら参照されず、'V'なら固有ベクトルが格納される
    double *tri_diag_evec;

    if (c == 'V')
    {
        tri_diag_evec = new double[tri_mat_dim * tri_mat_dim];
        vec_sub_init(tri_mat_dim * tri_mat_dim, tri_diag_evec);
    }

    /*=======Lanczos Algorithm=======*/
    for (int ls = 0; ls < tri_mat_dim; ls++)
    {
        if (err_checker)
        {
            ls_count = ls;
            /*省メモリのためのlanczosベクトルの更新*/
            if (ls > 0)
            {
                if (ls % 2 == 0)
                {
                    cblas_dscal(sub_mat_dim, -beta[ls - 1], u[1], 1);
                }
                else
                {
                    cblas_dscal(sub_mat_dim, -beta[ls - 1], u[0], 1);
                }
            }

            if (ls % 2 == 0)
            {
                // ls = 偶数step
                csc_sub_mvprod(u[0], u[1]);
                alpha[ls] = cblas_ddot(sub_mat_dim, u[1], 1, u[0], 1);
                cblas_daxpy(sub_mat_dim, -alpha[ls], u[0], 1, u[1], 1);
                if (ls != tri_mat_dim - 1)
                {
                    beta[ls] = cblas_dnrm2(sub_mat_dim, u[1], 1);
                    cblas_dscal(sub_mat_dim, 1. / beta[ls], u[1], 1);
                }
            }
            else
            {
                // ls = 奇数step
                csc_sub_mvprod(u[1], u[0]);
                alpha[ls] = cblas_ddot(sub_mat_dim, u[1], 1, u[0], 1);
                cblas_daxpy(sub_mat_dim, -alpha[ls], u[1], 1, u[0], 1);
                if (ls != tri_mat_dim - 1)
                {
                    beta[ls] = cblas_dnrm2(sub_mat_dim, u[0], 1);
                    cblas_dscal(sub_mat_dim, 1. / beta[ls], u[0], 1);
                }
            }
            /*===========================三重対角行列の数値対角(LAPACK)=============================*/
            vec_sub_init(tri_mat_dim, diag);
            vec_sub_init(tri_mat_dim - 1, sub_diag);
            int info = 0;

            if (ls % 2 == 0)
            {
                // 偶数step
                cblas_dcopy(tri_mat_dim, alpha, 1, diag, 1);
                cblas_dcopy(tri_mat_dim - 1, beta, 1, sub_diag, 1);

                if (ls < tri_mat_dim - 1)
                {
                    sub_diag[ls] = 0.;
                    if (c == 'N')
                    {
                        // 固有値のみを計算する場合
                        info =
                            LAPACKE_dstev(LAPACK_COL_MAJOR, 'N', ls + 1, diag,
                                          sub_diag, tri_diag_evec, ls + 1);
                    }
                    else
                    { // 固有ベクトルも計算する場合
                        info =
                            LAPACKE_dstev(LAPACK_COL_MAJOR, 'V', ls + 1, diag,
                                          sub_diag, tri_diag_evec, ls + 1);
                    }

                    if (info != 0)
                    {
                        std::cout << "@ls = " << ls
                                  << " , LAPACKE_detev errored." << std::endl;
                    }
                }
                else
                {
                    if (c == 'N')
                    {
                        // 固有値のみを計算する場合
                        info =
                            LAPACKE_dstev(LAPACK_COL_MAJOR, 'N', ls + 1, diag,
                                          sub_diag, tri_diag_evec, ls + 1);
                    }
                    else
                    {
                        // 固有ベクトルを計算する場合
                        info =
                            LAPACKE_dstev(LAPACK_COL_MAJOR, 'V', ls + 1, diag,
                                          sub_diag, tri_diag_evec, ls + 1);
                    }

                    if (info != 0)
                    {
                        std::cout << "@ls = " << ls
                                  << " , LAPACKE_detev errored." << std::endl;
                    }
                }
                cblas_dcopy(tri_mat_dim, diag, 1, eval_even, 1);
                if (info_ls == 'y')
                {
                    cout << "@ls = " << ls << " : eigen value = " << eval_even[0] << endl;
                }
                else if (info_ls == 's')
                {
                    std::cout << "@ls = " << ls << std::endl;
                }
            }
            else
            {
                cblas_dcopy(tri_mat_dim, alpha, 1, diag, 1);
                cblas_dcopy(tri_mat_dim - 1, beta, 1, sub_diag, 1);

                if (ls < tri_mat_dim - 1)
                {
                    sub_diag[ls] = 0.;
                    if (c == 'N')
                    {
                        // 固有値のみを計算する場合
                        info =
                            LAPACKE_dstev(LAPACK_COL_MAJOR, 'N', ls + 1, diag,
                                          sub_diag, tri_diag_evec, ls + 1);
                    }
                    else
                    {
                        // 固有ベクトルのみを計算する場合
                        info =
                            LAPACKE_dstev(LAPACK_COL_MAJOR, 'V', ls + 1, diag,
                                          sub_diag, tri_diag_evec, ls + 1);
                    }

                    if (info != 0)
                    {
                        std::cout << "@ls = " << ls
                                  << " , LAPACKE_detev's error." << std::endl;
                    }
                }
                else
                {
                    int info =
                        LAPACKE_dstev(LAPACK_COL_MAJOR, 'N', ls + 1, diag,
                                      sub_diag, tri_diag_evec, ls + 1);
                    if (info != 0)
                    {
                        std::cout << "@ls = " << ls
                                  << " , LAPACKE_detev's error." << std::endl;
                    }
                }
                cblas_dcopy(tri_mat_dim, diag, 1, eval_odd, 1);
                if (info_ls == 'y')
                {
                    std::cout << "@ls = " << ls
                              << " : eigen value = " << eval_odd[0]
                              << std::endl;
                }
                else if (info_ls == 's')
                {
                    cout << "@ls = " << ls << endl;
                }
            }

            /*============================収束状況の確認==============================*/
            if (ls > 0)
            {
                eps = abs(eval_even[0] - eval_odd[0]);
                if (info_ls == 'y')
                {
                    cout << "eps = " << std::setprecision(16) << eps << endl;
                }

                if (eps > err)
                    err_checker = true;
                else
                {
                    err_checker = false;
                    lanczos_check = true;
                }
            }
        }
        else
        {
            --ls_count;
            break;
        }
    }

    // /*========================基底状態の固有値===========================*/
    if (ls_count % 2 == 0)
        Eig.set_eval(eval_even[0]);
    else
        Eig.set_eval(eval_odd[0]);
    /*========================配列リソースのリリース part1===================*/
    delete[] eval_even;
    delete[] eval_odd;

    // 固有ベクトルを計算する
    if (c == 'V')
    {
        vec_sub_init(sub_mat_dim, u[0]);
        vec_sub_init(sub_mat_dim, u[1]);
        cblas_dcopy(sub_mat_dim, Eig.data(), 1, u[0], 1);

        for (int ls = 0; ls < ls_count + 2; ls++)
        {
            if (ls % 2 == 0)
            {
                if (ls == 0)
                    cblas_dscal(sub_mat_dim, tri_diag_evec[ls], Eig.data(), 1);
                else
                    cblas_daxpy(sub_mat_dim, tri_diag_evec[ls], u[0], 1, Eig.data(),
                                1);
            }
            else
            {
                cblas_daxpy(sub_mat_dim, tri_diag_evec[ls], u[1], 1, Eig.data(), 1);
            }

            if (ls > 0)
            {
                if (ls % 2 == 0)
                    cblas_dscal(sub_mat_dim, -beta[ls - 1], u[1], 1);
                else
                    cblas_dscal(sub_mat_dim, -beta[ls - 1], u[0], 1);
            }

            if (ls % 2 == 0)
            {
                csc_sub_mvprod(u[0], u[1]);
                cblas_daxpy(sub_mat_dim, -alpha[ls], u[0], 1, u[1], 1);
                if (ls != tri_mat_dim - 1)
                    cblas_dscal(sub_mat_dim, 1. / beta[ls], u[1], 1);
            }
            else
            {
                csc_sub_mvprod(u[1], u[0]);
                cblas_daxpy(sub_mat_dim, -alpha[ls], u[1], 1, u[0], 1);
                if (ls != tri_mat_dim - 1)
                    cblas_dscal(sub_mat_dim, 1. / beta[ls], u[0], 1);
            }
        }
        sub_sdz(sub_mat_dim, Eig.data());
    }
    /*==========================配列リソースのリリース part2
     * ====================*/
    for (int i = 0; i < 2; i++)
    {
        delete[] u[i];
    }
    delete[] u;
    delete[] alpha;
    delete[] beta;

    delete[] diag;
    delete[] sub_diag;
    if (c == 'V')
        delete[] tri_diag_evec;
}

std::string CSC_Sub_Hamiltonian::to_string() const
{
    std::ostringstream s;
    s << "\n\n";
    s << "Information of Hamiltonian\n";
    s << "---------------------------------------\n";
    s << "@Number of site   : " << tot_site_num << std::endl;
    s << "up spin site      : " << bit << endl;
    s << "Magnetization val : " << bit - tot_site_num * 0.5 << endl;
    s << "@Matrix dimension : " << sub_mat_dim << std::endl;
    s << "@non zero elements: " << nnz << endl;
    s << J.to_string() << std::endl;
    s << "@lanczos check    : " << lanczos_check << endl;
    s << "@lanczos step     : " << ls_count << endl;
    s << Eig.to_string() << endl;
    s << "------------------------------------------\n";

    return s.str();
}

ostream &operator<<(ostream &s, const CSC_Sub_Hamiltonian &h)
{
    return s << h.to_string();
}
