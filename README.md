# CSC_Hamiltonian.hpp

## 概要
`CSC_Hamiltonian.hpp`はスピン1/2 Heisenberg模型

$$\mathcal{H} = \sum_{i,j} \boldsymbol{S}_i \cdot \boldsymbol{S}_j$$

に対し、

- Hamiltonian行列を作成し、疎行列格納形式の1つであるCSC形式で保存する
- Hamiltonian行列に対してlanczos法を用いた厳密対角化による基底状態の固有値、固有ベクトルの計算

を行うためのライブラリです。動作についての注意事項があります。詳しくは動作環境の欄を御覧ください。
また、疎行列の格納形式については[疎行列フォーマット.pdf](疎行列フォーマット.pdf)

## アップデート予定
- 相互作用の情報を記入するファイルのフォーマットについて説明する欄を追加します。
- 使用できる関数について、その使い方と機能を説明する欄を追加します。
- OpenMPを用いたコードに修正し、並列処理を用いた計算の高速化を目指します。

## 動作環境
以下の環境で正常に動作することが確認できています。

- OS : Ubuntu 18.04 LTS
- コンパイラ : Intel C++ Compiler (GNU GCCでも動作します)

なお、使用するにあたっては以下の注意が必要です。
- 系のサイズに対して消費メモリ量がベキ的に増加するため、計算にあたってはご自身の計算環境でのメモリ量への十分な配慮が必要です。消費メモリ量については後述します。
- boost (boostのdynammic bitsetを使用します)のインストールが必要です。
- CBLAS、LAPACKのインストールと、適切なコンパイルオプションの指定が必要です。

## 使用方法
1. まずはライブラリを適当なディレクトリにダウンロードします。<br>
   今回は、`dir`というディレクトリを作成し、そこに
   ```
   git clone https://github.com/bokazu/class_CSC_Hamiltonian.git
   ```
   でダウンロードします。ディレクトリ構造は以下のようになります。なお、`EigenTest.cpp`、`HamiltonianTest.cpp`といったファイルは削除していただいてかまいません、下に記したファイルが存在すれば大丈夫です。
    ```
    dir
    ├── CSC_Hamiltonian
    │   ├── CSC_Hamiltonian.cpp
    │   └── CSC_Hamiltonian.hpp
    ├── EIGEN
    │   ├── EIGEN.cpp
    │   ├── EIGEN.hpp
    ├── Jset
    │   ├── Jset.cpp
    │   └── Jset.hpp
    ├── README.md
    ├── jset1.txt
    ├── jset2.txt
    └── makefile
    ``` 
2. 適当なcppファイルを作成します。個々では例としてmain.cppというフォルダを作成したとして説明します。
3. 系の相互作用の情報を記述したテキストファイルを用意します。ここでは、そのファイルに`jset.txt`という名前をつけて説明します。
4. main.cppに以下のように書き込みます。 
    1. Hamiltonian行列の作成だけを行う場合 

        ```cpp
        /*main.cpp*/

        #include "CSC_Hamiltonian/CSC_Hamiltonian.hpp"

        using namespace std;

        int main(){
            int site_num = 2; //系のサイト数
            CSC_Hamiltonian H("jset.txt" , site_num);

            H.hamiltonian(); //Hamiltonian行列を作成
        }
        ```

        これでHamiltonian行列が求まりました。もし、行列を表示したい場合は上記の`H.hamiltonian();`に続けて、
        ```
        H.print();
        ```
        も加えます。
    2. Hamiltonian行列を作成し、固有値も求める。
        
        ```cpp
        /*main.cpp*/

        #include "CSC_Hamiltonian/CSC_Hamiltonian.hpp"

        using namespace std;

        int main(){
                
            int site_num = 2; //系のサイト数
            CSC_Hamiltonian H("jset.txt" , site_num);

            H.hamiltonian(); // Hamiltonian行列を作成

            int tri_mat_dim = 4; //Lanczos法における最大反復回数を設定する
            H.csc_lanczos(tri_mat_dim);
        }
        ```

        CSC_Hamiltonianクラスのオブジェクト`H`の情報は以下のようにして標準出力することができます。

        ```cpp
         std::cout << H << std::endl;
        ```

        実行結果は以下のようになります。
        ```
        Information of Hamiltonian
        ----------------------------------------------------
        @Number of site    : 2
        @Matrix dimension : 4
        @nnz              : 6

        file name : jset1.txt
        ======================================
        i   j   J[i][j]
        --------------------------------------
        0   1   1
        ======================================

        @eigen_value = -0.75
        @eigen_vector
        eigen_vec[0] = 0
        ```

    3. Hamiltonian行列を作成し、固有値、固有ベクトルを求める
        ```cpp
        /*main.cpp*/

        #include "CSC_Hamiltonian/CSC_Hamiltonian.hpp"

        using namespace std;

        int main(){
                
            int site_num = 2; //系のサイト数
            CSC_Hamiltonian H("jset.txt" , site_num);

            H.hamiltonian(); // Hamiltonian行列を作成

            int tri_mat_dim = 4; //Lanczos法における最大反復回数を設定する
            H.csc_lanczos(tri_mat_dim, 'V');//引数として'V'を加える
        }
        ```
 5. 実行する<br>
    makefileを作成し、以下のように書き込みます。コンパイラ等は個人の環境に合わせてよしなに設定してください。
    ```makefile
    l_b = -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl

    program : main.o CSC_Hamiltonian.o Jset.o EIGEN.o
	icpc -o $@ $^ $(l_b)

    main.o : main.cpp
	icpc -c $< $(l_b)

    CSC_Hamiltonian.o : CSC_Hamiltonian/DNS_Hamiltonian.cpp
	icpc -c $< $(l_b)

    Jset.o : Jset/Jset.cpp
	icpc -c $< $(l_b)
	
    EIGEN.o : EIGEN/EIGEN.cpp
	icpc -c $< $(l_b)

    run : program
	    ./program

    clean:
	    rm -f ./program

    .PHONY : run clean
    ```
    記入し終えたらターミナルで
    ```bash
    make run
    ```
    と打つと実行されます。
## 消費メモリ
固有値計算を行う場合、計算する場合の各サイト数での消費メモリ量はおおよそ以下の表のようになります。

※`反復回数`とは<br>
...本ライブラリではlanczos法を用いた数値対角化を行いますが、この手法は反復法であり、固有値の精度が必要精度に達するまで繰り返しループを回します。そのループを抜けるまでに要した回数が反復回数となります。計算する前から何回で計算が収束するかはもちろんわかりませんが、本コードを動作させる場合は経験則から決め打ちします。


1. 固有値のみを計算する場合の消費メモリ量
    $x$をHamiltonian行列の非ゼロ要素数、 $n$を系の格子点数、 $l$をlanczos法における反復回数とすると

    $$消費メモリ = (12x + 48l + 20*2^n - 12) * 10^{-9} \text{GB}$$


2. 固有ベクトルも計算する場合のメモリ消費量
    $x$をHamiltonian行列の非ゼロ要素数、$n$を系の格子点数、$l$をlanczos法における反復回数とすると

    $$消費メモリ = (12x + 48l + 28*2^n - 12) * 10^{-9} \text{GB}$$
