import { FFTUtils } from "./fft_code/FFTUtils.js";

const canvas = document.getElementById('gridCanvas');
const context = canvas.getContext('2d');

// Create a color scale (adjust the color scheme as needed)
const colorScale = d3.scaleSequential(d3.interpolateInferno)
.domain([0, 1]);
canvas.height = window.innerHeight;
canvas.width = window.innerWidth;

const size_x = 256;
// const size_y = math.floor((canvas.height/canvas.width)*size_x);
const size_y = size_x;
const cellSize_x = canvas.width / size_x;
const cellSize_y = canvas.height / size_y;

// Setting up initial configuration

const R = 13;
const C = [[0,0,0,0,0,0,0.1,0.14,0.1,0,0,0.03,0.03,0,0,0.3,0,0,0,0],
[0,0,0,0,0,0.08,0.24,0.3,0.3,0.18,0.14,0.15,0.16,0.15,0.09,0.2,0,0,0,0],
[0,0,0,0,0,0.15,0.34,0.44,0.46,0.38,0.18,0.14,0.11,0.13,0.19,0.18,0.45,0,0,0],
[0,0,0,0,0.06,0.13,0.39,0.5,0.5,0.37,0.06,0,0,0,0.02,0.16,0.68,0,0,0],
[0,0,0,0.11,0.17,0.17,0.33,0.4,0.38,0.28,0.14,0,0,0,0,0,0.18,0.42,0,0],
[0,0,0.09,0.18,0.13,0.06,0.08,0.26,0.32,0.32,0.27,0,0,0,0,0,0,0.82,0,0],
[0.27,0,0.16,0.12,0,0,0,0.25,0.38,0.44,0.45,0.34,0,0,0,0,0,0.22,0.17,0],
[0,0.07,0.2,0.02,0,0,0,0.31,0.48,0.57,0.6,0.57,0,0,0,0,0,0,0.49,0],
[0,0.59,0.19,0,0,0,0,0.2,0.57,0.69,0.76,0.76,0.49,0,0,0,0,0,0.36,0],
[0,0.58,0.19,0,0,0,0,0,0.67,0.83,0.9,0.92,0.87,0.12,0,0,0,0,0.22,0.07],
[0,0,0.46,0,0,0,0,0,0.7,0.93,1,1,1,0.61,0,0,0,0,0.18,0.11],
[0,0,0.82,0,0,0,0,0,0.47,1,1,0.98,1,0.96,0.27,0,0,0,0.19,0.1],
[0,0,0.46,0,0,0,0,0,0.25,1,1,0.84,0.92,0.97,0.54,0.14,0.04,0.1,0.21,0.05],
[0,0,0,0.4,0,0,0,0,0.09,0.8,1,0.82,0.8,0.85,0.63,0.31,0.18,0.19,0.2,0.01],
[0,0,0,0.36,0.1,0,0,0,0.05,0.54,0.86,0.79,0.74,0.72,0.6,0.39,0.28,0.24,0.13,0],
[0,0,0,0.01,0.3,0.07,0,0,0.08,0.36,0.64,0.7,0.64,0.6,0.51,0.39,0.29,0.19,0.04,0],
[0,0,0,0,0.1,0.24,0.14,0.1,0.15,0.29,0.45,0.53,0.52,0.46,0.4,0.31,0.21,0.08,0,0],
[0,0,0,0,0,0.08,0.21,0.21,0.22,0.29,0.36,0.39,0.37,0.33,0.26,0.18,0.09,0,0,0],
[0,0,0,0,0,0,0.03,0.13,0.19,0.22,0.24,0.24,0.23,0.18,0.13,0.05,0,0,0,0],
[0,0,0,0,0,0,0,0,0.02,0.06,0.08,0.09,0.07,0.05,0.01,0,0,0,0,0]];

const cx = 20;
const cy = 20;
const scale = 1;

// Your JavaScript code goes here
let A = Array.from(new Array(size_y), _ => Array(size_x).fill(0));

// C = cells.zoom(scale, { order: 0 }); // Assuming C is your image data and scale is the zoom factor
// R *= scale;

// Assuming cx, cy are the coordinates where you want to place the zoomed image
// Ensure that cx + C.shape[0] and cy + C.shape[1] are within the bounds of A

// A[2*cx:2*cx + C.shape[0], 2*cy:2*cy + C.shape[1]] = np.fliplr(C)
// A[-2*cx:-2*cx + C.shape[0], -2*cy:-2*cy + C.shape[1]] = np.flipud(np.fliplr(C))
// A[10:10 + C.shape[0], 10:10 + C.shape[1]] = np.flipud(C)

const C_vertflip = C.map(row => row.slice().reverse());
const C_horflip = C.slice().reverse();

for (let i = 0; i < C.length; i++) {
  for (let j = 0; j < C[0].length; j++) {
    A[cx + i][cy + j] = C[i][j];
    A[3*cx + i][3*cy + j] = C_horflip[i][j];
  }
}

// Define convolution matrix K
function bell(x, m, s) {
    return math.exp(-math.pow((x - m) / s, 2) / 2);
}

// Create a meshgrid of indices
const X = math.range(-R+1, R + 1).toArray();
const Y = math.range(-R+1, R + 1).toArray();

const meshgrid = [];
for (let i = 0; i < X.length; i++) {
    const row = [];
    for (let j = 0; j < Y.length; j++) {
        row.push([X[i], Y[j]]);
    }
    meshgrid.push(row);
}

// Calculate the Euclidean norm of the indices
// const D = math.norm(meshgrid)/R;

const D = [];
for (let i = 0; i < X.length; i++) {
    const row = [];
    for (let j = 0; j < Y.length; j++) {
        row.push((math.norm(meshgrid[i][j]))/R);
    }
    D.push(row);
}
// console.log(D[25][25]);

// Create the kernel K
let K = math.map(D, value => (value < 1) ? bell(value, 0.5, 0.15) : 0);
K = math.divide(K, math.sum(K));
// console.log(K[12][12]);

// Function to map float values to RGB and draw the grid
function drawGrid(squareArray){
    // Render the array using colors on the canvas
    for (let i = 0; i < squareArray.length; i++) {
        for (let j = 0; j < squareArray[0].length; j++) {
            const color = d3.color(colorScale(squareArray[i][j]));
            context.fillStyle = color.toString();
            context.fillRect(j * cellSize_x, i * cellSize_y, cellSize_x, cellSize_y);
        }
    }
}

// Evolution
// const nj = require('numjs');
const steps = 500;
const m = 0.15;
const s = 0.015;
const T = 10;

function growth(U, m, s) {
    // return bell(U, m, s) * 2 - 1;
    return math.map(U, value => bell(value, m, s)*2-1);
}

function mod(n, m) {
    return ((n % m) + m) % m;
}

function convolution2D(A, B) {
    const aRows = A.length;
    const aCols = A[0].length;
    const bRows = B.length;
    const bCols = B[0].length;
    const output = new Array(aRows).fill(0).map(() => new Array(aCols).fill(0));

    const bRows_half = math.floor(bRows / 2);
    const bCols_half = math.floor(bCols / 2);
  
    for (let i = 0; i < aRows; i++) {
      for (let j = 0; j < aCols; j++) {
        for (let k = 0; k < bRows; k++) {
          for (let l = 0; l < bCols; l++) {
            const rowIdx = mod((i - bRows_half + k), aRows);
            const colIdx = mod((j - bCols_half + l), aCols);
            
            output[i][j] += A[rowIdx][colIdx] * B[bRows-1-k][bCols-1-l];
          }
        }
      }
    }
  
    return output;
}

function clip2DArray(arr) {
    const numRows = arr.length;
    const numCols = arr[0].length;
  
    // Create a new array to store the clipped values
    const clippedArray = new Array(numRows);
  
    for (let i = 0; i < numRows; i++) {
      clippedArray[i] = new Array(numCols);
  
      for (let j = 0; j < numCols; j++) {
        // Clip each value to be within the range [0, 1]
        clippedArray[i][j] = Math.min(Math.max(arr[i][j], 0), 1);
      }
    }
  
    return clippedArray;
}

function return2Darray(array, n_row, n_col){
  let output_array = Array.from(new Array(n_row), _ => Array(n_col).fill(0));

  for (let i = 0; i < n_row; i++) {
    for (let j = 0; j < n_col; j++) {
      output_array[i][j] = array[i*n_col+j];
    }
  }
  return output_array;
}

function simulate(step) {
    drawGrid(A);
    // Convolution
    // const U = convolution2D(A, K);

    const U = [];
    var nRows = size_y;
    var nCols = size_x;
    var data = new Array(nRows*nCols);
    for(var i=0;i<nRows;i++){
        for(var j=0;j<nCols;j++){
            data[i*nCols+j]=A[i][j];
        }
    }
    const conv_res = FFTUtils.convolute(data, K, size_y, size_x);
    while(conv_res.length) U.push(conv_res.splice(0,size_x));

    // Element-wise addition, scaling, and clipping
    const gr_mat = growth(U, m, s);
    const scaledGrowth = math.divide(gr_mat, T);
    const incrementedA = math.add(A, scaledGrowth);
    const clippedA = clip2DArray(incrementedA);

    // Update A with the clipped result
    A = clippedA;

    if (step < steps) {
        requestAnimationFrame(() => simulate(step + 1));
    }
}

// Draw the grid with the initial configuration
simulate(0);