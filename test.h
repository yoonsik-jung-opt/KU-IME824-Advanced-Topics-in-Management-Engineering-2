//
// Created by yoonsikjung on 2022/10/06.
//

#ifndef IME824_TEST_H
#define IME824_TEST_H

#include "LinAlg.h"
#include "utils.h"

using namespace LinearAlgebra;

void LUEx1(){
    Matrix A;
    A.addRow(vect{2, 3, 1, 5});
    A.addRow(vect{6, 13, 5, 19});
    A.addRow(vect{2, 19, 10, 23});
    A.addRow(vect{4, 10, 11, 31});
    A.LUDecomposition();

}

void LUEx2(){
    // https://courses.physics.illinois.edu/cs357/sp2020/notes/ref-9-linsys.html
    Matrix A;
    A.addRow({1,2,2});
    A.addRow({4,4,2});
    A.addRow({4,6,4});
    A.addRow({2,2,2});
    A.LUDecomposition();
}

void linearSystemEx1(){
    Matrix A;
    A.addRow(vect{1, 2, -4, 5});
    A.addRow(vect{2, 1, -6, 8});
    A.addRow(vect{4, -1, -12, 13});
    A.solveLinearEquation();
}

void linearSystemEx2(){
    Matrix A;
    A.addRow(vect{8, 7, 9, 5});
    A.addRow(vect{6, 7, 9, 8});
    A.addRow(vect{4, 3, 3, 1});
    A.addRow(vect{2, 1, 1, 0});
    A.solveLinearEquation();
}

void outerProductEx(){
    matx X = outerProduct({2, 2, -1}, {2, 2, -1});
    Matrix Z(X);
    Z.print();
}

void matrixMulEx(){
    Matrix A;
    A.addRow({2./3,2./3,-1./3});
    A.addRow({2./3,-1./3,2./3});
    A.addRow({-1./3,2./3,2./3});

    Matrix B(A.T());

    Matrix C(matMul(A.getMatrix(), B.getMatrix()));
    C.print();

}

void eigenEx1(){
    Matrix A;
    A.addRow({2,1});
    A.addRow({1,2});
//    Matrix B(GramSchmidt(A.getMatrix()));
//    B.print();
//    A.QRDecomposition();
    A.eigenValue();
}


void eigenEx2(){
    Matrix A;
//    A.addRow(vect{8,7,9,5});
//    A.addRow(vect{6,7,9,8});
//    A.addRow(vect{4,3,3,1});
//    A.addRow(vect{2,1,1,4});

//    A.addRow({3,4,-2});
//    A.addRow({1,4,-1});
//    A.addRow({2,6,-1});



//    A.addRow({2./3,2./3,-1./3});
//    A.addRow({2./3,-1./3,2./3});
//    A.addRow({-1./3,2./3,2./3});
    Matrix B(matMul(A.T(), A.getMatrix()));
    B.print();
    B.eigenValue();
    B.eigenVector();

}

void operTest(){
    vect a = {0, 1, 2, 3};
    vect b = {1, 2, 3, 4};
    vect sum = vectorSum(a, b);
    vect min = vectorMinus(a, b);

    vect c = {2, 1};
    vect d = {1, 2};
    vect proj = projection(c, d);
}

void SVDTest(){
    Matrix A;
//    A.addRow({3,4,-2});
//    A.addRow({1,4,-1});
//    A.addRow({2,6,-1});

    A.addRow({2,1,6});
    A.addRow({3,3,4});
    A.addRow({1,5,-1});

//    A.addRow(vect{8,7,9,5});
//    A.addRow(vect{6,7,9,8});
//    A.addRow(vect{4,3,3,1});
//    A.addRow(vect{2,1,1,4});
//    A.addRow(vect{1,1,1,1});

    A.SVD(true);
}

void armadiloTest(){
    mat temp;
    temp.load("wisc.dat");

    mat X = temp.submat(0, 0, temp.n_rows-1, temp.n_cols-2);
    vec y = temp.col(9);

    mat coeff;
    mat score; // projected data
    vec latent; // eigenValues
    vec tsquared;

    princomp(coeff, score, latent, tsquared, X);
//    X.print();
//    y.print();
//    cout << endl;
//    coeff.print();
//    cout << endl;
//    score.print();
//    cout << endl;
    latent.print();
//    cout << endl;
//    tsquared.print();
//    mat sub = score.submat(0,0,score.n_rows-1, 2);
//    sub.save("PCA_3feat.csv", csv_ascii);
//    y.save("PCA_y.csv", csv_ascii);
//    mat X(5, 4, fill::randu);
//
//    mat coeff;
//    mat score;
//    vec latent;
//    vec tsquared;
//
//    princomp(coeff, score, latent, tsquared, X);

}

void efficiencyTest(){
    clock_t start, end, start1, end1;

    // LU Decomposition
    mat A;
    A.load("wisc.dat");
    A = A * A.t();
//    A = A.submat(0,0,A.n_rows-1, A.n_cols-2);
//    A={{1,2,3,4,5}, {2,1,4,5,6}, {8,1,2,6,4}, {2,8,9,0,1}, {5,7,2,3,4}};
    mat L, U, P;
    start = clock();
    lu(L, U, P, A);
    end = clock();

    matx tempB;

    readData<double>(tempB, "wisc.dat", '\t');
    Matrix B(matMul(tempB, Matrix(tempB).T()));

//    B.addRow({1,2,3,4,5});
//    B.addRow({2,1,4,5,6});
//    B.addRow({8,1,2,6,4});
//    B.addRow({2,8,9,0,1});
//    B.addRow({5,7,2,3,4});

    start1 = clock();
    B.LUDecomposition();
    end1 = clock();
    cout << "LU" << endl;
    cout << "Armadillo : " <<(double) (end - start) / CLOCKS_PER_SEC << endl;
    cout << "Implementation : " <<(double) (end1 - start1) / CLOCKS_PER_SEC << endl;

    // Qr

    mat Q, R;
    start = clock();
    qr(Q, R, A);
    end = clock();

    start1 = clock();
    B.QRDecomposition();
    end1 = clock();
    cout << "QR" << endl;
    cout << "Armadillo : " <<(double) (end - start) / CLOCKS_PER_SEC << endl;
    cout << "Implementation : " << (double) (end1 - start1) / CLOCKS_PER_SEC << endl;
    // EigenValue

    cx_vec eigval;
    cx_mat eigvec;

    start = clock();
    eig_gen(eigval, eigvec, A);
    end = clock();

    start1 = clock();
    B.eigenValue();
    B.eigenVector();
    end1 = clock();

    cout << "Eigen" << endl;
    cout << "Armadillo : " <<(double) (end - start) / CLOCKS_PER_SEC << endl;
    cout << "Implementation : " << (double) (end1 - start1) / CLOCKS_PER_SEC << endl;

    // SVD

    mat UU;
    vec s;
    mat V;

    start = clock();
    svd(UU,s,V,A);
    end = clock();

    start1 = clock();
    B.SVD();
    end1 = clock();

    cout << "SVD" << endl;
    cout << "Armadillo : " <<(double) (end - start) / CLOCKS_PER_SEC << endl;
    cout << "Implementation: " << (double) (end1 - start1) / CLOCKS_PER_SEC << endl;

}

#endif //IME824_TEST_H
