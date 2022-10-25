# IME824



## Installation

```bash
g++ main.cpp -o prog -arch x86_64  -std=c++11 -O2 -larmadillo
./prog
```

## Usage
0. Load Data
```c++
matx temp;
readData<double>(temp, "wisc.dat", '\t');
Matrix A(temp);
```

1. Matrix Multiplication
```c++
   Matrix A;
   A.addRow({2./3,2./3,-1./3});
   A.addRow({2./3,-1./3,2./3});
   A.addRow({-1./3,2./3,2./3});

   Matrix B(A.T());

   Matrix C(matMul(A.getMatrix(), B.getMatrix()));
```
2. LU Decomposition
```c++
Matrix A;
A.addRow(vect{2, 3, 1, 5});
A.addRow(vect{6, 13, 5, 19});
A.addRow(vect{2, 19, 10, 23});
A.addRow(vect{4, 10, 11, 31});
A.LUDecomposition();

```
3. QR Decomposition
```c++
Matrix A;
A.addRow(vect{2, 3, 1, 5});
A.addRow(vect{6, 13, 5, 19});
A.addRow(vect{2, 19, 10, 23});
A.addRow(vect{4, 10, 11, 31});
A.QRDecomposition();
```
4. Eigen Value & Vector
```c++
Matrix A;
A.addRow(vect{8,7,9,5});
A.addRow(vect{6,7,9,8});
A.addRow(vect{4,3,3,1});
A.addRow(vect{2,1,1,4});
A.eigenValue();
A.eigenVector();
```
5. Singular Vector Decomposition
```c++
Matrix A;
A.addRow(vect{8,7,9,5});
A.addRow(vect{6,7,9,8});
A.addRow(vect{4,3,3,1});
A.addRow(vect{2,1,1,4});
A.SVD();
```

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License
[MIT](https://choosealicense.com/licenses/mit/)