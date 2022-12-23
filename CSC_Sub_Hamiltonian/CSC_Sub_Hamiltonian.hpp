#ifndef ___Class_CSC_Sub_Hamiltonian
#define ___Class_CSC_Sub_Hamiltonian

#include <mkl.h>
#include <iomanip>

#include "../EIGEN/EIGEN.hpp"
#include "../Jset/Jset.hpp"

class CSC_Sub_Hamiltonian
{
public:
    std::string jset_filename;
    int tot_site_num;
    int bit;
    int mat_dim;
    int sub_mat_dim;
    int nnz;
    int *row;
    int *col_ptr;
    double *val;
    int *bm;
    int *gbm;
    int ls_count;
    bool lanczos_check;
    Jset J;
    EIGEN Eig;

    // コンストラクタ
    CSC_Sub_Hamiltonian(std::string filename, int site, int b)
        : jset_filename(filename),
          tot_site_num(site),
          bit(b),
          mat_dim(1 << tot_site_num),
          sub_mat_dim(comb(site, b)),
          nnz(0),
          row(new int[1]),
          col_ptr(new int[sub_mat_dim + 1]),
          val(new double[1]),
          bm(new int[sub_mat_dim]),
          gbm(new int[1]),
          ls_count(0),
          lanczos_check(false),
          J(filename),
          Eig(1)
    {
        std::cout << comb(site, b) << std::endl;
        set_J();
        std::cout << "CSC_Sub_Hamiltonian::constructed.\n";
    }

    // デストラクタ
    ~CSC_Sub_Hamiltonian()
    {
        delete[] row;
        delete[] col_ptr;
        delete[] val;
        delete[] bm;
        delete[] gbm;
        std::cout << "CSC_Hamiltonian::destructed.\n";
    }

    /*----------------------ゲッタ-----------------------*/
    // Jsetの情報を書き込んだファイルの名前を返す
    std::string jsetfile() const { return jset_filename; }

    // 系のサイト数
    int site() const { return tot_site_num; }

    // 行列の次元を返す
    int sub_dim() const { return sub_mat_dim; }

    // HamiltonianのNon zero要素の個数を返す
    int num_nnz() const { return nnz; }

    // Hamiltonianのrow[i]の値を返す
    int at_row(int i) const { return row[i]; }

    // Hamiltonianのcol[i]の値を返す
    int at_col_ptr(int i) const { return col_ptr[i]; }

    // Hamiltonianのval[i]の値を返す
    double at_val(int i) const { return val[i]; }

    // lanczos法での反復回数を返す
    int num_ls() const { return ls_count; }
    /*----------------------------------------------------------------*/

    /*------------------------------------その他メンバ関数------------------------------*/
    // row valの初期化を行う
    void init();

    // CSC形式のrowとvalの要素数をiに変更し、0で初期化する
    void set_sub_CscElem(int i);

    // 部分空間を張るスピン状態を確認する
    void subspace_check();

    // スピン演算子
    void csc_sub_spin(int m, int site, int &row_index, int &col_ptr_val, double &szz);

    // 部分空間におけるHamiltonian行列の非ゼロ要素の個数を数え上げる
    void count_sub_nnz();

    // 部分空間におけるHamiltonian行列要素の計算と配列への格納を行う
    void csc_sub_hamiltonian();

    void csc_sub_lanczos(int tri_mat_dim, char c = 'N', char info_ls = 'n');

    // CSC形式での行列ベクトル積の計算を行う
    void csc_sub_mvprod(double *u_i, double *u_j);

    // Jsetを設定する
    void set_J() { J.set(); }

    int comb(int n, int r);
    /*-------------------------------------------------------------------------------------*/

    /*-------------------------------------標準出力関係---------------------------------*/
    // 標準出力において、倍精度を何桁まで出力するかを指定する
    int precision = 5;

    // precisionのsetter
    void set_precision(int i) { precision = i; }

    // Hamiltonian行列の全要素を標準出力する
    void print() const;

    // COO_Hamiltonianオブジェクトの文字列表現を返却する
    std::string to_string() const;
    /*--------------------------------------------------------------------*/
};

// 出力ストリームにhを挿入する
std::ostream &operator<<(std::ostream &s, const CSC_Sub_Hamiltonian &h);

template <typename T>
void vec_sub_init(int dim, T *vec)
{
    for (int i = 0; i < dim; i++)
    {
        vec[i] = 0;
    }
}

#endif