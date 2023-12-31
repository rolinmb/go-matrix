package main

import (
  "fmt"
  "time"
  "math"
  "errors"
)

type Matrix struct {
  m int
  n int
  mat [][]float64
}

func printMatrix(mtrx Matrix) {
  fmt.Printf("%dx%d Matrix: [\n",mtrx.m, mtrx.n)
  for i := 0; i < mtrx.m; i++ {
    fmt.Println("\t",mtrx.mat[i])
  }
  fmt.Println("]")
}

func makeMatrix(data [][]float64) (Matrix, error) {
  var mRows int
  var nCols int
  if len(data) >= 1 {
    if len(data[0]) == 0 {
      return Matrix{}, errors.New("\n[ERROR] All input rows must be atleast length 1")
    }
    mRows = len(data)
    nCols = len(data[0])
    for i := 0; i < mRows; i++ {
      if len(data[i]) == 0 { 
        return Matrix{}, errors.New("\n[ERROR] All input rows must be atleast length 1")
      }
      if len(data[i]) != nCols {
	      return Matrix{}, errors.New("\n[ERROR] Input row width is not uniform")
      }
    }
    return Matrix{
      m: mRows,
      n: nCols,
      mat: data, 
    }, nil
  } else {
    return Matrix{}, errors.New("\n[ERROR] Matrix must have atleast one row (M >= 1)")
  }
}

func matrixTranspose(mtrx Matrix) Matrix {
  result := make([][]float64, mtrx.n)
  for i := 0; i < mtrx.n; i++ {
    col := make([]float64, mtrx.m)
    for j := 0; j < mtrx.m; j++ {
      col[j] = mtrx.mat[j][i]
    }
    result[i] = col
  }
  return Matrix{
    m: mtrx.n,
    n: mtrx.m,
    mat: result,
  }
}

func qrDecomp(mtrx Matrix) (Matrix, Matrix) {
  Q := make([][]float64, mtrx.m)
  R := make([][]float64, mtrx.n)
  for i := 0; i < mtrx.m; i++ {
    Q[i] = make([]float64, mtrx.n)
    R[i] = make([]float64, mtrx.n)
  }
  // Gram-Schmidt Process for QR-Decomposition
  for j := 0; j < mtrx.n; j++ {
    v := make([]float64, mtrx.m)
    for i := 0; i < mtrx.m; i++ {
      v[i] = mtrx.mat[i][j]
    }
    for i := 0; i < j; i++ {
      R[i][j] = 0.0
      dotp := 0.0
      for k := 0; k < mtrx.m; k++ {
        dotp += Q[k][i] * mtrx.mat[k][j]
      }
      for k := 0; k < mtrx.m; k++ {
        v[k] -= dotp * Q[k][i]
      }
    }
    vNorm := 0.0
    for _, val := range v {
      vNorm += val * val
    }
    vNorm = math.Sqrt(vNorm)
    for i := 0; i < mtrx.m; i++ {
      Q[i][j] = v[i] / vNorm
    }
    for i := 0; i < mtrx.n; i++ {
      dotp := 0.0
      for k := 0; k < mtrx.m; k++ {
        dotp += Q[k][j] * mtrx.mat[k][i]
      }
      R[j][i] = dotp
    }
  }
  return Matrix{
    m: mtrx.m,
    n: mtrx.n,
    mat: Q,
  }, Matrix{
    m: mtrx.n,
    n: mtrx.n,
    mat: R,
  }
}
// QR-Decomposition Algorithm for finding Eigenvalues
func qrAlgo(mtrx Matrix, iters int) []float64 {
  n0 := mtrx.m
  copyMtrx, err := makeMatrix(mtrx.mat)
  if err != nil {
    fmt.Println(err)
  }
  for it := 0; it < iters; it++ {
    Q, R := qrDecomp(copyMtrx)
    newMtrx, err := matrixMult(R, Q)
    if err != nil {
      fmt.Println(err)
    }
    copyMtrx = newMtrx
  }
  eigenvals := make([]float64, n0)
  for i := 0; i < n0; i++ {
    eigenvals[i] = copyMtrx.mat[i][i]
  }
  return eigenvals
}

func getEigenVecs(mtrx Matrix) ([]Matrix, error) {
  eigenVals := qrAlgo(mtrx, 100)
  eigenVecs := make([]Matrix, 0)
  for _, eigenVal := range eigenVals {
    id := make([][]float64, mtrx.m)
    for i := 0; i < mtrx.m; i++ {
      id[i] = make([]float64, mtrx.m)
      id[i][i] = 1
    }
    idMat := Matrix{
      m: len(id),
      n: len(id[0]),
      mat: id,
    }
    amli, err := matrixSub(mtrx, matrixScale(idMat, eigenVal))
    if err != nil {
      return nil, err
    }
    eigenVec, err := solveHmgSys(amli)
    if err != nil {
      return nil, err
    }
    eigenVecs = append(eigenVecs, eigenVec)
  }
  return eigenVecs, nil
}

func matrixAdd(a,b Matrix) (Matrix, error) {	
  if a.m != b.m || a.n != b.n {
    return Matrix{}, errors.New("\n[ERROR] Can only add Matrix objects of same MxN dimensions")
  } else {
    sum := make([][]float64, a.m)
    for i := 0; i < a.m; i++ {
      sum[i] = make([]float64, a.n)
      for j := 0; j < a.n; j++ {
        sum[i][j] = a.mat[i][j] + b.mat[i][j]
      }
    }
    return Matrix{
      m: a.m,
      n: a.n,
      mat: sum,
    }, nil
  }
}

func matrixSub(a,b Matrix) (Matrix, error) {
  if a.m != b.m || a.n != b.n {
    return Matrix{}, errors.New("\n[ERROR] Can only subtract Matrix objects of same MxN dimensions")
  } else {
    dif := make([][]float64, a.m)
    for i := 0; i < a.m; i++ {
      dif[i] = make([]float64, a.n)
      for j := 0; j < a.n; j++ {
        dif[i][j] = a.mat[i][j] - b.mat[i][j]
      }
    }
    return Matrix{
      m: a.m,
      n: a.n,
      mat: dif,
    }, nil
  }
}

func matrixScale(mtrx Matrix, scalar float64) Matrix {
  result := make([][]float64, mtrx.m)
  for i := 0; i < mtrx.m; i++ {
    result[i] = make([]float64, mtrx.n)
    for j := 0; j < mtrx.n; j++ {
      result[i][j] = mtrx.mat[i][j] * scalar
    }
  }
  return Matrix{
    m: mtrx.m,
    n: mtrx.n,
    mat: result,
  }
}

func matrixMult(a,b Matrix) (Matrix, error) {
  if a.n != b.m {
    return Matrix{}, errors.New("\n[ERROR] Incompatible matrix dimensions for multiplication")
  }
  result := make([][]float64, a.m)
  for i := 0; i < a.m; i++ {
    result[i] = make([]float64, b.n)
    for j := 0; j < b.n; j++ {
      result[i][j] = 0.0
      for k := 0; k < a.n; k++ {
	      result[i][j] += a.mat[i][k] * b.mat[k][j]
      }
      rounded := math.Round(result[i][j]*1e10)/1e10
      if math.Abs(rounded) < 1e-11 {
        result[i][j] = 0.0
      }
    }
  }
  return Matrix{
    m: a.m,
    n: b.n,
    mat: result,
  }, nil
}

func subMatrix(mtrx Matrix, i,j int) Matrix {
  sub := make([][]float64, mtrx.m-1)
	for k := range sub {
		sub[k] = make([]float64, mtrx.m-1)
	}
	row := 0
	for r := 0; r < mtrx.m; r++ {
		if r == i {
			continue
		}
		col := 0
		for c := 0; c < mtrx.m; c++ {
			if c == j {
				continue
			}
			sub[row][col] = mtrx.mat[r][c]
			col++
		}
		row++
	}
	return Matrix{
    m: mtrx.m-1,
    n: mtrx.n-1,
    mat: sub,
  }
}

func matrixDet(mtrx Matrix) (float64, error) {
  if mtrx.m != mtrx.n {
    return 0, errors.New("Cannot find determinant of a non-square Matrix")  
  }
  if mtrx.m <= 1 {
    return 0, errors.New("Determinant is not defined for a 1x1 matrix")
  }
  if mtrx.m == 2 {
    return mtrx.mat[0][0]*mtrx.mat[1][1] - mtrx.mat[0][1]*mtrx.mat[1][0], nil
  }
  if mtrx.m == 3 {
    return mtrx.mat[0][0]*mtrx.mat[1][1]*mtrx.mat[2][2] + mtrx.mat[0][1]*mtrx.mat[1][2]*mtrx.mat[2][0] + mtrx.mat[0][2]*mtrx.mat[1][0]*mtrx.mat[2][1] - mtrx.mat[0][2]*mtrx.mat[1][1]*mtrx.mat[2][0] - mtrx.mat[0][0]*mtrx.mat[1][2]*mtrx.mat[2][1] - mtrx.mat[0][1]*mtrx.mat[1][0]*mtrx.mat[2][2], nil
  }
  det := 0.0
  sign := 1.0
  for j := 0; j < mtrx.m; j++ {
    subMtrx := subMatrix(mtrx, 0, j)
    subDet, _ := matrixDet(subMtrx)
    det += sign * mtrx.mat[0][j] * subDet
    sign *= -1.0
  }
  return det, nil
}

func hasInverse(mtrx Matrix) bool {
  if mtrx.m != mtrx.n { // must be a square matrix
    return false 
  }
  det, err := matrixDet(mtrx)
  if err != nil {
    fmt.Println(err)
  }
  if det == 0 { // must have non-zero determinant
    return false
  }
  return true
}

func matrixInv(mtrx Matrix) (Matrix, error) {
  if !hasInverse(mtrx) {
    return Matrix{}, errors.New("Cannot find the inverse of supplied matrix")
  }
  id := make([][]float64, mtrx.m)
  copyMat := make([][]float64, mtrx.m)
  for i := 0; i < mtrx.m; i++ {
    id[i] = make([]float64, mtrx.m)
    id[i][i] = 1
    copyMat[i] = make([]float64, mtrx.n)
    copy(copyMat[i], mtrx.mat[i])
  }
  for i := 0; i < mtrx.m; i++ {
    if copyMat[i][i] == 0 {
      return Matrix{}, errors.New("Matrix is singular, cannot find its inverse")
    }
    div := copyMat[i][i]
    for j := 0; j < mtrx.m; j++ {
      copyMat[i][j] /= div
      id[i][j] /= div
    }
    for k := 0; k < mtrx.m; k ++ {
      if k != i {
        mul := copyMat[k][i]
        for l := 0; l < mtrx.m; l++ {
          copyMat[k][l] -= mul * copyMat[i][l]
          id[k][l] -= mul * id[i][l]
        }
      }
    }
  }
  return Matrix{
    m: mtrx.m,
    n: mtrx.n,
    mat: id,
  }, nil
}

func solveHmgSys(mtrx Matrix) (Matrix, error) {
  augMtrx := make([][]float64, mtrx.m)
  for i := 0; i < mtrx.m; i++ {
    augMtrx[i] = make([]float64, mtrx.n + 1)
    for j := 0; j < mtrx.n; j++ {
      augMtrx[i][j] = mtrx.mat[i][j]
    }
  }
  for i := 0; i < mtrx.m; i++ {
    augMtrx[i][mtrx.n] = 0
  }
  for i := 0; i < mtrx.m; i++ {
    if augMtrx[i][i] == 0 {
      for j := i + 1; j < mtrx.m; j++ {
        if augMtrx[j][i] != 0 {
          augMtrx[i], augMtrx[j] = augMtrx[j], augMtrx[i]
          break
        }
      }
    }
    factor := 1.0 / augMtrx[i][i]
    for j := 0; j < mtrx.m; j++ {
      augMtrx[i][j] *= factor
    }
    for j := 0; j < mtrx.m; j++ {
      if j != i {
        factor = -augMtrx[j][i]
        for k := 0; k <= mtrx.m; k++ {
          augMtrx[j][k] += factor * augMtrx[i][k]
        }
      }
    }
  }
  result := make([][]float64, mtrx.m)
  for i := 0; i < mtrx.m; i++ {
    result[i] = make([]float64, 1)
    result[i][0] = augMtrx[i][mtrx.n]
  }
  return Matrix{
    m: mtrx.m,
    n: 1,
    mat: result,
  }, nil
}

func dotProduct(a,b Matrix) (Matrix, error) {
  if a.n != b.m {
    return Matrix{}, errors.New("Incompatable dimensions for dot product.")
  }
  result := make([][]float64, a.m)
  for i := 0; i < a.m; i++ {
    result[i] = make([]float64, b.n)
    for j := 0; j < b.n; j++ {
      result[i][j] = 0.0
      for k := 0; k < a.n; k++ {
        result[i][j] += a.mat[i][k] * b.mat[k][j]
      }
    }
  }
  return Matrix{
    m: a.m,
    n: b.n,
    mat: result,
  }, nil
}

type Tensor struct {
  d int
  matrices [][]Matrix
}

func printTensor(tnsr Tensor) {
  fmt.Printf("\n\nTensor (DEPTH = %d) ((", tnsr.d)
  idx := 0
  for i, mtrxSlice := range tnsr.matrices {
    for j, mtrx := range mtrxSlice {
      fmt.Printf("\nMatrix idx = %d => {%d, %d}:", idx, i, j)
      printMatrix(mtrx)
      idx++
    }
  }
  fmt.Println("))")
}

func makeTensor(matrices [][]Matrix) (Tensor, error) {
  dDepth := len(matrices)
  if dDepth >= 1 {
    if len(matrices[0]) == 0 || len(matrices[0][0].mat) == 0 {
      return Tensor{}, errors.New("\n[ERROR] All input matrices must have atleast 1 row and 1 column")
    }
    return Tensor{
      d: dDepth,
      matrices: matrices,
    }, nil
  } else {
    return Tensor{}, errors.New("\n[ERROR] Tensor must have atleast 1 matrix (D >= 1)")
  }
}

func tensorProduct(t1,t2 Tensor) (Tensor, error) {
  if len(t1.matrices) == 0 || len(t2.matrices) == 0 {
    return Tensor{}, errors.New("Cannot take tensor product of empty tensor.")
  }
  if t1.d != t2.d {
    return Tensor{}, errors.New("Tensor product is only defined for tensors of the same dimension / depth")
  }
  var result [][]Matrix
  for _, sliceM1 := range t1.matrices {
    var mtrxSlice []Matrix
    for _, m1 := range sliceM1 {
      for _, sliceM2 := range t2.matrices {
        for _, m2 := range sliceM2 {
          product, err := matrixMult(m1, m2)
          if err != nil {
            return Tensor{}, err
          }
          mtrxSlice = append(mtrxSlice, product)
        }
      }
    }
    result = append(result, mtrxSlice)
  }
  return Tensor{
    d: t1.d,
    matrices: result,
  }, nil
}

func runMatrixTest() {
  testStart := time.Now()
  // TESTS ON INVALID STRUCTS
  data := [][]float64 {
    {0.0, 1.0},
    {0.0, 1.0},
    {0.0, 1.0},
    {0.0, 1.0},
    {0.0, 1.0},
    {0.0, 1.0, 2.0, 3.0},
  }
  badMat, err := makeMatrix(data)
  if err != nil {
    fmt.Println(err)
    fmt.Println("\nInvalid Matrix ", data)
    fmt.Println("\nReturned Matrix{}: ", badMat)
  }
  data = [][]float64 {
    {},
    {0.0},
    {1.0},
  }
  badMat, err = makeMatrix(data)
  if err != nil {
    fmt.Println(err)
    fmt.Println("\nInvalid Matrix ", data)
    fmt.Println("\nReturned Matrix{}: ", badMat)
  }
  data = [][]float64 {
    {0.0, 1.0},
    {},
    {0.0, 1.0},
  }
  badMat, err = makeMatrix(data)
  if err != nil {
    fmt.Println(err)
    fmt.Println("\nInvalid Matrix ", data)
    fmt.Println("\nReturned Matrix{}: ", badMat)
  }
  // TESTS ON VALID STRUCTS
  data0 := [][]float64 {
    {0.0},
  }
  goodMat0, err := makeMatrix(data0)
  if err != nil {
    fmt.Println(err)
  }
  printMatrix(goodMat0)  
  dataA := [][]float64 {
    {0.0, 1.0, 2.0, 3.0},
    {4.0, 5.0, 6.0, 7.0},
    {8.0, 9.0, 10.0, 11.0},
  }
  goodMatA, err := makeMatrix(dataA)
  if err != nil {
    fmt.Println(err)
  }
  printMatrix(goodMatA)
  dataB := [][]float64 {
    {-1.0, -1.0, -1.0, -1.0},
    {-1.0, -1.0, -1.0, -1.0},
    {-1.0, -1.0, -1.0, -1.0},
  }
  goodMatB, err := makeMatrix(dataB)
  if err != nil {
    fmt.Println(err)
  }
  printMatrix(goodMatB)
  sumMat, err := matrixAdd(goodMatA, goodMatB)
  if err != nil {
    fmt.Println(err)
  }
  fmt.Println("\nValid Sum Result 1 (A+B):")
  printMatrix(sumMat)
  sumMat, err = matrixAdd(goodMatB, goodMatA)
  if err != nil {
    fmt.Println(err)
  }
  fmt.Println("\nValid Sum Result 2 (B + A):")
  printMatrix(sumMat)
  difMatAB, err := matrixSub(goodMatA, goodMatB)
  if err != nil {
    fmt.Println(err)
  }
  fmt.Println("\nValid Difference Result 1 (A - B):")
  printMatrix(difMatAB)
  difMatBA, err := matrixSub(goodMatB, goodMatA)
  if err != nil {
    fmt.Println(err)
  }
  fmt.Println("\nValid Difference Result 2 (B - A):")
  printMatrix(difMatBA)
  dataC := [][]float64 {
    {0.0, 1.0, 2.0},
    {3.0, 4.0, 5.0},
  }
  dataD := [][]float64 {
    {0.0, 1.0},
    {2.0, 3.0},
    {4.0, 5.0},
  }
  goodMatC, err := makeMatrix(dataC)
  if err != nil {
    fmt.Println(err)
  }
  printMatrix(goodMatC)
  goodMatD, err := makeMatrix(dataD)
  if err != nil {
    fmt.Println(err)
  }
  printMatrix(goodMatD)
  multMatCD, err := matrixMult(goodMatC, goodMatD)
  if err != nil {
    fmt.Println(err)
  }
  fmt.Println("\nValid Multiplication Result (C * D):")
  printMatrix(multMatCD)
  multMatDC, err := matrixMult(goodMatD, goodMatC)
  if err != nil {
    fmt.Println(err)
  }
  fmt.Println("\nValid Multiplication Result (D * C):")
  printMatrix(multMatDC)
  dataE := [][]float64 {
    {4.0, 3.0},
    {3.0, 2.0},
  }
  goodMatE, err := makeMatrix(dataE)
  if err != nil {
    fmt.Println(err)
  }
  printMatrix(goodMatE)
  detTestE, err := matrixDet(goodMatE)
  if err != nil {
    fmt.Println(err)
  }
  fmt.Printf("\nMatrix E Determinant: %f", detTestE)
  invMatE, err := matrixInv(goodMatE)
  if err != nil {
    fmt.Println(err)
  }
  fmt.Println("\nInverse of Matrix E:")
  printMatrix(invMatE)
  shouldBeId, err := matrixMult(goodMatE, invMatE)
  if err != nil {
    fmt.Println(err)
  }
  fmt.Println("\nProduct of E * Inverse of E:")
  printMatrix(shouldBeId)
  dataF := [][]float64 {
    {2.0, 3.0, 1.0, 5.0, 4.0},
    {1.0, 2.0, 4.0, 2.0, 3.0},
    {3.0, 1.0, 2.0, 3.0, 1.0},
    {4.0, 2.0, 3.0, 1.0, 2.0},
    {1.0, 4.0, 2.0, 3.0, 5.0},
  }
  goodMatF, err := makeMatrix(dataF)
  if err != nil {
    fmt.Println(err)
  }
  printMatrix(goodMatF)
  detTestF, err := matrixDet(goodMatF)
  if err != nil {
    fmt.Println(err)
  }
  
  fmt.Printf("\nMatrix F Determinant: %f",detTestF)
  invMatF, err := matrixInv(goodMatF)
  if err != nil {
    fmt.Println(err)
  }
  fmt.Println("Inverse of Matrix F:\n")
  printMatrix(invMatF)
  shouldBeId, err = matrixMult(goodMatF, invMatF)
  if err != nil {
    fmt.Println(err)
  }
  fmt.Println("\nProduct of F * Inverse of F:")
  printMatrix(shouldBeId)
  dataG := [][]float64 {
    {1, 2, 3},
    {4, 5, 6},
  }
  dataH := [][]float64 {
    {7, 8},
    {9, 10},
    {11, 12},
  }
  goodMatG, err := makeMatrix(dataG)
  if err != nil {
    fmt.Println(err)
  }
  printMatrix(goodMatG)
  goodMatH, err := makeMatrix(dataH)
  if err != nil {
    fmt.Println(err)
  }
  printMatrix(goodMatH)
  dotTest, err := dotProduct(goodMatG, goodMatH)
  if err != nil {
    fmt.Println(err)
  }
  fmt.Println("\nDot Product of G & F:")
  printMatrix(dotTest)
  // TENSOR TESTS
  dataI := [][]float64 {
    {1.0, 2.0},
    {3.0, 4.0},
  }
  dataJ := [][]float64 {
    {5.0, 6.0},
    {7.0, 8.0},
  }
  dataK := [][]float64 {
    {9.0, 10.0},
    {11.0, 12.0},
  }
  dataL := [][]float64 {
    {13.0, 14.0},
    {15.0, 16.0},
  }
  goodMatI, err := makeMatrix(dataI)
  if err != nil {
    fmt.Println("Error: ", err)
  }
  goodMatJ, err := makeMatrix(dataJ)
  if err != nil {
    fmt.Println("Error: ", err)
  }
  goodMatK, err := makeMatrix(dataK)
  if err != nil {
    fmt.Println("Error: ", err)
  }
  goodMatL, err := makeMatrix(dataL)
  tensorMats1 := [][]Matrix {
    {goodMatI, goodMatJ},
    {goodMatK, goodMatL},
  }
  tensorMats2 := [][]Matrix {
    {goodMatJ, goodMatK},
    {goodMatL, goodMatI},
  }
  testTensor1, err := makeTensor(tensorMats1)
  if err != nil {
    fmt.Println(err)
  }
  fmt.Println("\nTest Tensor T1:")
  printTensor(testTensor1)
  testTensor2, err := makeTensor(tensorMats2)
  if err != nil {
    fmt.Println(err)
  }
  fmt.Println("\nTest Tensor T2:")
  printTensor(testTensor2)
  prodTensor, err := tensorProduct(testTensor1, testTensor2)
  if err != nil {
    fmt.Println(err)
  }
  fmt.Println("\nTensor Product of T1 and T2:")
  printTensor(prodTensor)
  // END OF TESTS
  testEnd := time.Now()
  testTime := testEnd.Sub(testStart).Nanoseconds()
  fmt.Printf("\nrunMatrixTest() Execution Time: %d nanoseconds", testTime)
}
// MAIN EXECUTION
func main() {
  runMatrixTest()
}
