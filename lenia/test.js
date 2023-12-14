import { FFT } from "./fft_code/fftlib.js";
import { FFTUtils } from "./fft_code/FFTUtils.js";

const A = [[1,2,3,4],[1,2,3,4],[1,2,3,4],[1,2,3,4],[1,2,3,4],[1,2,3,4],[1,2,3,4],[1,2,3,4]];
var ker = [[1,2],[1,2]];

const aRows = A.length;
const aCols = A[0].length;
const kRows = ker.length;
const kCols = ker[0].length;
const output = new Array(aRows).fill(0).map(() => new Array(aCols).fill(0));

function mod(n, m) {
    return ((n % m) + m) % m;
}

for (let i = 0; i < aRows; i++) {
    for (let j = 0; j < aCols; j++) {
      for (let k = 0; k < kRows; k++) {
        for (let l = 0; l < kCols; l++) {
          const rowIdx = mod((i - math.floor(kRows / 2) + k), aRows);
          const colIdx = mod((j - math.floor(kCols / 2) + l), aCols);
          
          output[i][j] += A[rowIdx][colIdx] * ker[kRows-1-k][kCols-1-l];
        }
      }
    }
}

console.log(output);

// const test_signal = [1,2,3,4,5];
// const test_ker = [1,1];

// const out_arr = FFTUtils.convolute(A, ker, 4, 4);
// FFTUtils.convolute2DI(A, ker, 4, 4)
// console.log(A);

// test code

var nRows = A.length;
var nCols = A[0].length;
var data = new Array(nRows*nCols);
for(var i=0;i<nRows;i++){
    for(var j=0;j<nCols;j++){
        data[i*nCols+j]=A[i][j];
    }
}

// var kn = 2;
// var kernel = new Array(kn);
// for(var i=0;i<kn;i++){
//   kernel[i]=new Array(kn);
//   for(var j=0;j<kn;j++){
//       kernel[i][j]=ker[i][j];
//   }
// }

function return2Darray(array, n_row, n_col){
  let output_array = Array.from(new Array(n_row), _ => Array(n_col).fill(0));

  for (let i = 0; i < n_row; i++) {
    for (let j = 0; j < n_col; j++) {
      output_array[i][j] = array[i*n_col+j];
    }
  }
  return output_array;
}

var convolutedData = FFTUtils.convolute(data, ker, nRows, nCols);
console.log(return2Darray(convolutedData, nRows, nCols));







