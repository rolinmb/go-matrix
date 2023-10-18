package main

import (
  "fmt"
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

func runMatrixTest() {
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
    fmt.Println("\nReturned Matrix{} struct: ", badMat)
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
    fmt.Println("\nReturned Matrix{} struct: ", badMat)
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
    fmt.Println("\nReturned Matrix{} struct: ", badMat)
  }
  data = [][]int {{0}}
  goodMat, err := makeMatrix(data)
  if err != nil {
    fmt.Println(err)
  } else {
    fmt.Println("\nValid Matrix{} struct: ", goodMat)
    printMatrix(goodMat)
  }  
  data = [][]int{
    {0, 1, 2, 3},
    {4, 5, 6, 7},
    {8, 9, 10, 11},
  }
  goodMat, err = makeMatrix(data)
  if err != nil {
    fmt.Println(err)
  } else {
    fmt.Println("\nValid Matrix{} struct (A): ", goodMat)
    printMatrix(goodMat)
  }
  data0 := [][]int {
    {-1, -1, -1, -1},
    {-1, -1, -1, -1},
    {-1, -1, -1, -1},
  }
  goodMat0, err := makeMatrix(data0)
  if err != nil {
    fmt.Println(err)
  } else {
    fmt.Println("\nValid Matrix{} struct (B): ", goodMat0)
    printMatrix(goodMat0)
  }
  sumMat, err := matrixAdd(goodMat, goodMat0)
  if err != nil {
    fmt.Println(err)
  } else {
    fmt.Println("\nValid Sum Result 1 (A+B): ", sumMat)
    printMatrix(sumMat)
  }
  sumMat, err = matrixAdd(goodMat0, goodMat)
  if err != nil {
    fmt.Println(err)
  } else {
    fmt.Println("\nValid Sum Result 2 (B+A): ", sumMat)
    printMatrix(sumMat)
  }
  difMat, err := matrixSub(goodMat, goodMat0)
  if err != nil {
    fmt.Println(err)
  } else {
    fmt.Println("\nValid Difference Result 1 (A-B): ", difMat)
    printMatrix(difMat)
  }
  difMat0, err := matrixSub(goodMat0, goodMat)
  if err != nil {
    fmt.Println(err)
  } else {
    fmt.Println("\nValid Difference Result 2 (B-A): ", difMat0)
    printMatrix(difMat)
  }
}

func main() {
  runMatrixTest()
}
