//
// Created by yoonsikjung on 2022/10/06.
//

#ifndef IME824_TEST_H
#define IME824_TEST_H

#include "LinAlg.h"

using namespace LinearAlgebra;

void LUEx1(){
    Matrix A;
    A.addRow(vec{2,3,1,5});
    A.addRow(vec{6,13,5,19});
    A.addRow(vec{2,19,10,23});
    A.addRow(vec{4,10,11,31});
    A.LUDecomposition();

}

void LUEx2(){
    // https://courses.physics.illinois.edu/cs357/sp2020/notes/ref-9-linsys.html
    Matrix A;
    A.addRow({1,2,2});
    A.addRow({4,4,2});
    A.addRow({4,6,4});
    A.LUDecomposition();
}

void linearSystemEx1(){
    Matrix A;
    A.addRow(vec{1,2,-4,5});
    A.addRow(vec{2,1,-6,8});
    A.addRow(vec{4,-1, -12, 13});
    A.solveLinearEquation();
}

void linearSystemEx2(){
    Matrix A;
    A.addRow(vec{8,7,9,5});
    A.addRow(vec{6,7,9,8});
    A.addRow(vec{4,3,3,1});
    A.addRow(vec{2,1,1,0});
    A.solveLinearEquation();
}

void outerProductEx(){
    mat X = outerProduct({2, 2, -1}, {2, 2, -1});
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
    A.addRow(vec{8,7,9,5});
    A.addRow(vec{6,7,9,8});
    A.addRow(vec{4,3,3,1});
    A.addRow(vec{2,1,1,4});
//    A.addRow({2./3,2./3,-1./3});
//    A.addRow({2./3,-1./3,2./3});
//    A.addRow({-1./3,2./3,2./3});
    A.eigenValue();

}

void operTest(){
    vec a = {0,1,2,3};
    vec b = {1,2,3,4};
    vec sum = vectorSum(a, b);
    vec min = vectorMinus(a, b);

    vec c = {2, 1};
    vec d = {1, 2};
    vec proj = projection(c, d);


}

#endif //IME824_TEST_H
