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

func makeMatrix(data [][]int) (Matrix, error) {
  var mRows int
  var nCols int
  if len(data) >= 1 {
    if len(data[0]) == 0 {
      return Matrix{}, errors.New("All input rows must be atleast length 1")
    }
    mRows = len(data)
    nCols = len(data[0])
    for i := 0; i < mRows; i++ {
      if len(data[i]) == 0 { 
        return Matrix{}, errors.New("All input rows must be atleast length 1")
      }
      if len(data[i]) != nCols {
	      return Matrix{}, errors.New("Input row width is not uniform")
      }
    }
    return Matrix{
      m: mRows,
      n: nCols,
      mat: data, 
    }, nil
  } else {
    return Matrix{}, errors.New("Matrix must have atleast one row (M >= 1)")
  }
}

func main() {
  data := [][]int{
    {0, 1},
    {0, 1},
    {0, 1},
    {0, 1},
    {0, 1},
    {},
  }
  badMat, err := makeMatrix(data)
  if err != nil {
    fmt.Println(err)
    fmt.Println("Invalid Matrix: \n",data)
    fmt.Println("Returned Matrix{} struct: ", badMat)
  }
  data = [][]int{
    {0, 1, 2, 3},
    {4, 5, 6, 7},
    {8, 9, 10, 11},
  }
  goodMat, err := makeMatrix(data)
  if err != nil {
    fmt.Println(err)
  } else {
    fmt.Println("Valid Matrix: \n",goodMat.mat)
  } 
}
