package main

import (
  "fmt"
  "time"
  "errors"
)

type Matrix struct {
  m int
  n int
  mat [][]int
}

func printMatrix(mtrx Matrix) {
  fmt.Printf("%dx%d Matrix: [\n",mtrx.m, mtrx.n)
  for i := 0; i < mtrx.m; i++ {
    fmt.Println("\t",mtrx.mat[i])
  }
  fmt.Println("]")
}

func makeMatrix(data [][]int) (Matrix, error) {
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

func matrixAdd(a,b Matrix) (Matrix, error) {	
  if a.m != b.m || a.n != b.n {
    return Matrix{}, errors.New("\n[ERROR] Can only add Matrix objects of same MxN dimensions")
  } else {
    sum := make([][]int, a.m)
    for i := 0; i < a.m; i++ {
      sum[i] = make([]int, a.n)
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
    dif := make([][]int, a.m)
    for i := 0; i < a.m; i++ {
      dif[i] = make([]int, a.n)
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

func matrixMult(a,b Matrix) (Matrix, error) {
  if a.n != b.m {
    return Matrix{}, errors.New("\n[ERROR] Incompatible matrix dimensions for multiplication")
  }
  result := make([][]int, a.m)
  for i := 0; i < a.m; i++ {
    result[i] = make([]int, b.n)
    for j := 0; j < b.n; j++ {
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

func subMatrix(mtrx Matrix, i,j int) Matrix {
  sub := make([][]int, mtrx.m-1)
	for k := range sub {
		sub[k] = make([]int, mtrx.m-1)
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

func matrixDet(mtrx Matrix) (int, error) {
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
  det := 0
  sign := 1
  for j := 0; j < mtrx.m; j++ {
    subMtrx := subMatrix(mtrx, 0, j)
    subDet, _ := matrixDet(subMtrx)
    det += sign * mtrx.mat[0][j] * subDet
    sign = -sign
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
  id := make([][]int, mtrx.m)
  for i := 0; i < mtrx.m; i++ {
    id[i] = make([]int, mtrx.m)
    id[i][i] = 1
  }
  for i := 0; i < mtrx.m; i++ {
    if mtrx.mat[i][i] == 0 {
      return Matrix{}, errors.New("Matrix is singular, cannot find its inverse")
    }
    div := mtrx.mat[i][i]
    for j := 0; j < mtrx.m; j++ {
      mtrx.mat[i][j] /= div
      id[i][j] /= div
    }
    for k := 0; k < mtrx.m; k ++ {
      if k != i {
        mul := mtrx.mat[k][i]
        for l := 0; l < mtrx.m; l++ {
          mtrx.mat[k][l] -= mul * mtrx.mat[i][l]
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


func runMatrixTest() {
  testStart := time.Now()
  data := [][]int{
    {0, 1},
    {0, 1},
    {0, 1},
    {0, 1},
    {0, 1},
    {0, 1, 2, 3},
  }
  badMat, err := makeMatrix(data)
  if err != nil {
    fmt.Println(err)
    fmt.Println("\nInvalid Matrix ", data)
    fmt.Println("\nReturned Matrix{}: ", badMat)
  }
  data = [][]int {
    {},
    {0},
    {1},
  }
  badMat, err = makeMatrix(data)
  if err != nil {
    fmt.Println(err)
    fmt.Println("\nInvalid Matrix ", data)
    fmt.Println("\nReturned Matrix{}: ", badMat)
  }
  data = [][]int {
    {0, 1},
    {},
    {0, 1},
  }
  badMat, err = makeMatrix(data)
  if err != nil {
    fmt.Println(err)
    fmt.Println("\nInvalid Matrix ", data)
    fmt.Println("\nReturned Matrix{}: ", badMat)
  }
  data = [][]int {
    {0},
  }
  goodMat, err := makeMatrix(data)
  if err != nil {
    fmt.Println(err)
  }
  fmt.Println("\nValid Matrix{}: ", goodMat)
  printMatrix(goodMat)  
  data = [][]int{
    {0, 1, 2, 3},
    {4, 5, 6, 7},
    {8, 9, 10, 11},
  }
  goodMat, err = makeMatrix(data)
  if err != nil {
    fmt.Println(err)
  }
  fmt.Println("\nValid Matrix{} (A): ", goodMat)
  printMatrix(goodMat)
  data0 := [][]int {
    {-1, -1, -1, -1},
    {-1, -1, -1, -1},
    {-1, -1, -1, -1},
  }
  goodMat0, err := makeMatrix(data0)
  if err != nil {
    fmt.Println(err)
  }
  fmt.Println("\nValid Matrix{} (B): ", goodMat0)
  printMatrix(goodMat0)
  sumMat, err := matrixAdd(goodMat, goodMat0)
  if err != nil {
    fmt.Println(err)
  }
  fmt.Println("\nValid Sum Result 1 (A+B): ", sumMat)
  printMatrix(sumMat)
  sumMat, err = matrixAdd(goodMat0, goodMat)
  if err != nil {
    fmt.Println(err)
  }
  fmt.Println("\nValid Sum Result 2 (B+A): ", sumMat)
  printMatrix(sumMat)
  difMat, err := matrixSub(goodMat, goodMat0)
  if err != nil {
    fmt.Println(err)
  }
  fmt.Println("\nValid Difference Result 1 (A-B): ", difMat)
  printMatrix(difMat)
  difMat0, err := matrixSub(goodMat0, goodMat)
  if err != nil {
    fmt.Println(err)
  }
  fmt.Println("\nValid Difference Result 2 (B-A): ", difMat0)
  printMatrix(difMat0)
  data1 := [][]int {
    {0, 1, 2},
    {3, 4, 5},
  }
  data2 := [][]int {
    {0, 1},
    {2, 3},
    {4, 5},
  }
  goodMat1, err := makeMatrix(data1)
  if err != nil {
    fmt.Println(err)
  }
  fmt.Println("\nValid Matrix{} (C): ", goodMat1)
  printMatrix(goodMat1)
  goodMat2, err := makeMatrix(data2)
  if err != nil {
    fmt.Println(err)
  }
  fmt.Println("\nValid Matrix{} struct (D): ", goodMat2)
  multMat, err := matrixMult(goodMat1, goodMat2)
  if err != nil {
    fmt.Println(err)
  }
  fmt.Println("\nValid Multiplication Result (C*D): ", multMat)
  data3 := [][]int {
    {2, 4, 1, 3, 5},
		{5, 6, 2, 7, 2},
		{1, 2, 6, 3, 1},
		{3, 1, 3, 6, 2},
    {2, 2, 2, 2, 2},
  }
  goodMat3, err := makeMatrix(data3)
  if err != nil {
    fmt.Println(err)
  }
  fmt.Println("\nValid Matrix{} (E): ", goodMat3)
  printMatrix(goodMat3)
  detTest, err := matrixDet(goodMat3)
  if err != nil {
    fmt.Println(err)
  }
  fmt.Println("\nMatrix E Determinant: ", detTest)
  invMat, err := matrixInv(goodMat3)
  if err != nil {
    fmt.Println(err)
  }
  fmt.Println("\nInverse of Matrix E: ", invMat)
  printMatrix(invMat)
  testInvMult, err := matrixMult(goodMat3, invMat)
  if err != nil {
    fmt.Println(err)
  }
  fmt.Println("\nProduct of E * Inverse of E: ", testInvMult)
  printMatrix(testInvMult)
  testEnd := time.Now()
  testTime := testEnd.Sub(testStart).Nanoseconds()
  fmt.Printf("runMatrixTest() Execution Time: %d nanoseconds", testTime)
}

func main() {
  runMatrixTest()
}
