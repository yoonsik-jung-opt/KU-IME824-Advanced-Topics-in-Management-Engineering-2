//
// Created by yoonsikjung on 2022/09/21.
//

#ifndef IME824_LINALG_H
#define IME824_LINALG_H
#include <vector>
#include <limits>
#include <iostream>
#include <algorithm>
#include <cmath>

using namespace std;
typedef vector<vector<double>> mat;
typedef vector<double> vec;

namespace LinearAlgebra{
    mat gaussElimination(mat);
    tuple<mat, mat> LUDecomposition(mat A);
    double innerProduct(vec &vec1, vec &vec2);
    void solveLinearEquation(mat Uc);
    mat outerProduct(vec vec1, vec vec2);
    mat matMul(mat, mat);
    mat GramSchmidt(mat A);


    class Matrix {
    private:
        mat matrix;
        vector<int> rowBasisInd;
        vector<int> colBasisInd;
        int nRows = 0;
        int nCols = 0;
        int rank = 0;

        mat L;
        mat U;
        mat Q;
        mat R;

        vector<double> eigneValues;
        vector<vector<double>> eigneVectors;

    public:
        Matrix() = default;

        Matrix(int rows, int cols, int initValue = 0) {
            this->nRows = rows;
            this->nCols = cols;

            mat matrix(rows, vec(cols, initValue));
            this->matrix = matrix;
        }

        Matrix(mat A){
            this->matrix = A;
            this->nRows = A.size();
            this->nCols = A[0].size();
        }


        vector<double> getRow(int rowInd) {
            return this->matrix[rowInd];
        }

        double getElement(int rowInd, int colInd) {
            return this->matrix[rowInd][colInd];
        }

        void addRow(vec row) {
            if(nCols == 0)
                nCols = row.size();
            assert(row.size() == nCols);
            this->nRows += 1;

            this->matrix.push_back(row);
        }

        void QRDecomposition(){
            Matrix Q(GramSchmidt(this->matrix));
            Matrix R(matMul(Q.T(), this->matrix));
            this->Q = Q.getMatrix();
            this->R = R.getMatrix();

//            Q.print();
//            R.print();
//            Matrix X(matMul(Q.getMatrix(), R.getMatrix()));
//            X.print();

        }

        void eigenValue(){
            // assure the square matrix
            assert(this->matrix.size() == this->matrix[0].size());

            if(Q.empty()){
                QRDecomposition();
            }
            // QR Method
            double eps = 1e-10;
            mat tempQ = Q;
            mat temp = this->matrix;
            double threshold = numeric_limits<double>::max();
            int iter = 0;
            int iterMax = 100000;

            // check 0 in diagonal
            for(int i = 0; i < nRows; i++){
                if(temp[i][i] == 0){
                    assert(false);
                }
            }

            while(threshold >= eps){
                tempQ = GramSchmidt(temp);
                temp = matMul(Matrix(tempQ).T(), temp);
                temp = matMul(temp, tempQ);

                double upperTri = 0;
                double lowerTri = 0;

                for(int i = 0; i < nRows; i++){
                    for(int j = 0; j < nCols; j++){
                        if(i > j){
                            // lower
                            lowerTri += temp[i][j] * temp[i][j];
                        }
                        else if(i < j){
                            // upper
                            upperTri += temp[i][j] * temp[i][j];
                        }
                    }
                }
                threshold = upperTri > lowerTri ? lowerTri : upperTri;

                iter++;
                if(iter > iterMax)
                    break;
            }

            Matrix X(temp);
            X.print();
            cout << iter << endl;

            for(int i = 0; i < nRows; i++){
                eigneValues.push_back(temp[i][i]);
            }



        }

        void SVD(){
            Matrix A1(matMul(this->matrix, Matrix(this->matrix).T()));
            Matrix A2(matMul(Matrix(this->matrix).T(), this->matrix));

        }

        void updateRank(){
            // calculate the rank of nRows of the matrix
            mat C; // basis of columns
            mat R; // basis of nRows
            int r = 0; // rank

            // find basis of nRows
            mat rowElimination = gaussElimination(this->matrix);
            vector<int> basisInd;
            for(int i = 0; i < this->nRows; i++){
                vec row = rowElimination[i];
                bool isAllZero = all_of(row.begin(), row.end(), [](int x){return x==0;});
                if(!isAllZero){
                    basisInd.push_back(i);
                    r++;
                }
            }

            mat colElimination = gaussElimination(this->T());

        }

        void LUDecomposition(){
            auto LU = LinearAlgebra::LUDecomposition(this->matrix);
            this->L = get<0>(LU);
            this->U = get<1>(LU);
        }

        void solveLinearEquation(){
            if(this->U.empty()){
                this->LUDecomposition();
            }

            if(U.size() < U[0].size()){
                LinearAlgebra::solveLinearEquation(this->U);
            }
            else{
                cout << "### There is no RHS vector c ###" << endl;
                cout << "# of rows : " << this->U.size() << " # of cols : " << this->U[0].size() << endl;
            }


        }

        mat getMatrix(){
            return this->matrix;
        }

        mat T(){ //transpose
            mat temp = Matrix(this->nCols, this->nRows, 0).getMatrix();
            for(int i = 0; i < this->nCols; i++){
                for(int j = 0; j < this->nRows; j++){
                    temp[i][j] = this->matrix[j][i];
                }
            }
            return temp;

        }

        void print() {
            for (int r = 0; r < this->nRows; r++) {
                for (int c = 0; c < this->nCols; c++) {
                    cout << this->getElement(r, c) << "\t";
                }
                cout << endl;
            }
        }
    };

    mat gaussElimination(mat A){ // Only Upper Triangular Matrix
        // non-zero rows are basis

        int nRows = A.size();
        int nCols = A[0].size();

        for(int i = 0; i < nRows; i++){
            double div = A[i][i];
            if(div == 0.)
                continue;

            for(int j = i; j < nCols; j++){
                A[i][j] /= div;
            }

            for(int j = i+1; j< nRows; j++){
                double divRow = A[j][i] / A[i][i];

                for(int k = 0; k < nCols; k++){
                    A[j][k] -= A[i][k] * divRow;
                }
            }

        }

        Matrix B(A);
        B.print();

        return A;
    }



    vec projection(vec &v1, vec &v2){
        // project v2 onto v1
        double coef = innerProduct(v1, v2) / innerProduct(v1, v1);
        vec res;
        for(auto v : v1){
            res.push_back(coef * v);
        }
        return res;
    }

    vec vectorSum(vec v1, vec v2){
        assert(v1.size() == v2.size());
        int sz = v1.size();
        vec res(sz);
        for(int i = 0; i < sz; i++){
            res[i] = v1[i] + v2[i];
        }
        return res;
    }

    vec vectorMinus(vec v1, vec v2){
        assert(v1.size() == v2.size());
        int sz = v1.size();
        vec res(sz);
        for(int i = 0; i < sz; i++){
            res[i] = v1[i] - v2[i];
        }
        return res;
    }

    mat GramSchmidt(mat A){
        int nRows = A.size();
        int nCols = A[0].size();

        mat Q;
        Q.resize(nRows);
        for(int i = 0; i < nRows; i++){
            Q[i].resize(nCols);
        }

        for(int i = 0; i < nCols; i++){
            vec colVec(nRows);
            for(int j = 0; j < nRows; j++){
                colVec[j] = A[j][i];
            }

            if(i == 0){
                for(int j = 0; j < nRows; j++){
                    Q[j][i] = colVec[j];
                }
                continue;
            }
            else{
                for(int k = 0; k < i; k++){
                    vec u(nRows);
                    for(int j = 0; j < nRows; j++){
                        u[j] = Q[j][k];
                    }
                    colVec = vectorMinus(colVec, projection(u, colVec));
                }
                for(int j = 0; j < nRows; j++) {
                    Q[j][i] = colVec[j];
                }
            }

        }

        // normalize
        for(int i = 0; i < nCols; i++){
            double colSum = 0;
            for(int j =0; j < nRows; j++){
                colSum += Q[j][i] * Q[j][i];
            }
            double coef = 1./sqrt(colSum);
            for(int j =0; j < nRows; j++){
                Q[j][i] *= coef;
            }
        }

        return Q;
    }

    tuple<mat, mat> LUDecomposition(mat A){
        int nRows = A.size();
        int nCols = A[0].size();
//    assert(nRows == nCols);

        mat L;
        L.resize(nRows);
        for(int i = 0; i <nRows; i++){
            L[i].resize(nCols);
            L[i][i] = 1;
        }


        // Permutation
        if(A[0][0] == 0){
            int maxRowIdx = 0;
            double maxVal = A[0][0];
            for(int i = 0; i < nRows; i++){
                if(A[i][0] > maxVal) {
                    maxRowIdx = i;
                    maxVal = A[i][0];
                }
            }

            vec temp = A[0];
            A[0] = A[maxRowIdx];
            A[maxRowIdx] = temp;
        }

        for(int i = 0; i < nRows; i++){
            double base = A[i][i];
            for(int j = i+1; j < nRows; j++){
                double multiplier = A[j][i] / base;
                L[j][i] = multiplier;
                for(int k = 0; k < nCols; k++){
                    A[j][k] -= A[i][k] * multiplier;
                }
            }
        }

        cout << "L:Lower Triangular Matrix" << endl;
        Matrix C(L);
        C.print();
        cout << "U:Upper Triangular Matrix" <<endl;
        Matrix B(A);
        B.print();
//        cout << "LU" <<endl;
//        Matrix T(matMul(L, A));
//        T.print();

        tuple<mat, mat> res = make_tuple(L, A);

        return res;
    }


    double innerProduct(vec &vec1, vec &vec2) {
        assert(vec2.size() == vec1.size());
        double res = 0;

        for (int i = 0; i < vec1.size(); i++) {
            res += vec1[i] * vec2[i];
        }
        return res;
    }

    mat outerProduct(vec vec1, vec vec2){ // vec1 * vec2^T
        int sz1 = vec1.size();
        int sz2 = vec2.size();
        assert(sz1 == sz2);

        mat res;
        res.resize(sz1);
        for(int i = 0; i < sz1; i++){
            res[i].resize(sz1);
        }

        for(int i = 0; i < sz1; i++){
            double mul = vec1[i];
            std::transform(vec2.begin(), vec2.end(), res[i].begin(), [mul](double &c){return c*mul;});
        }

        return res;
    }

    mat matMul(mat A, mat B){
        int lRows = A.size();
        int lCols = A[0].size();
        int rRows = B.size();
        int rCols = B[0].size();
        assert(lCols == rRows);

        mat res;
        res.resize(lRows);
        for(int i = 0; i < lRows; i++){
            res[i].resize(rCols);
        }

        for(int i = 0; i < lRows; i++){
            for(int j = 0; j < rCols; j++){
                vec row = A[i];
                vec col;
                for(int k = 0; k < rRows; k++){
                    col.push_back(B[k][j]);
                }
                res[i][j] = innerProduct(row, col);
            }
        }

        return res;

    }

    void solveLinearEquation(mat Uc){ //Uc : Upper triangular Matrix(U) + RHS(c)
        // using back substitution
        int nRows = Uc.size(); // 2
        int nCols = Uc[0].size(); // 3

        vec res;
        res.resize(nCols-1, 0); // remove RHS column
        int iter = 1;
        for(int i = nRows-1; i >= 0; i--){
            double RHS = Uc[i][nCols-1];
            double sum = 0.;
            for(int j = 0; j < nCols - 1; j++){
                sum += res[j] * Uc[i][j];
            }
            res[i] = (RHS - sum)/Uc[i][nCols-1-iter];
            iter++;
        }
        cout << "Solution of Linear System" << endl;
        for(auto s:res){
            cout << s << " ";
        }
    }
}



#endif //IME824_LINALG_H
