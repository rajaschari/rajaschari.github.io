import { FFT } from './fftlib.js';

/**
 * Calculates the inverse of a 2D Fourier transform
 *
 * @param ft
 * @param ftRows
 * @param ftCols
 * @return
 */
function ifft2DArray(ft, ftRows, ftCols) {
  let tempTransform = new Array(ftRows * ftCols);
  let nRows = ftRows / 2;
  let nCols = (ftCols - 1) * 2;
  // reverse transform columns
  FFT.init(nRows);
  let tmpCols = { re: new Array(nRows), im: new Array(nRows) };
  let iRow, iCol;
  for (iCol = 0; iCol < ftCols; iCol++) {
    for (iRow = nRows - 1; iRow >= 0; iRow--) {
      tmpCols.re[iRow] = ft[iRow * 2 * ftCols + iCol];
      tmpCols.im[iRow] = ft[(iRow * 2 + 1) * ftCols + iCol];
    }
    //Unnormalized inverse transform
    FFT.bt(tmpCols.re, tmpCols.im);
    for (iRow = nRows - 1; iRow >= 0; iRow--) {
      tempTransform[iRow * 2 * ftCols + iCol] = tmpCols.re[iRow];
      tempTransform[(iRow * 2 + 1) * ftCols + iCol] = tmpCols.im[iRow];
    }
  }

  // reverse row transform
  let finalTransform = new Array(nRows * nCols);
  FFT.init(nCols);
  let tmpRows = { re: new Array(nCols), im: new Array(nCols) };
  let scale = nCols * nRows;
  for (iRow = 0; iRow < ftRows; iRow += 2) {
    tmpRows.re[0] = tempTransform[iRow * ftCols];
    tmpRows.im[0] = tempTransform[(iRow + 1) * ftCols];
    for (iCol = 1; iCol < ftCols; iCol++) {
      tmpRows.re[iCol] = tempTransform[iRow * ftCols + iCol];
      tmpRows.im[iCol] = tempTransform[(iRow + 1) * ftCols + iCol];
      tmpRows.re[nCols - iCol] = tempTransform[iRow * ftCols + iCol];
      tmpRows.im[nCols - iCol] = -tempTransform[(iRow + 1) * ftCols + iCol];
    }
    //Unnormalized inverse transform
    FFT.bt(tmpRows.re, tmpRows.im);

    let indexB = (iRow / 2) * nCols;
    for (iCol = nCols - 1; iCol >= 0; iCol--) {
      finalTransform[indexB + iCol] = tmpRows.re[iCol] / scale;
    }
  }
  return finalTransform;
}

/**
 * Calculates the fourier transform of a matrix of size (nRows,nCols) It is
 * assumed that both nRows and nCols are a power of two
 *
 * On exit the matrix has dimensions (nRows * 2, nCols / 2 + 1) where the
 * even rows contain the real part and the odd rows the imaginary part of the
 * transform
 * @param data
 * @param nRows
 * @param nCols
 * @return
 */
function fft2DArray(data, nRows, nCols) {
  let iRow, iCol;
  let ftCols = nCols / 2 + 1;
  let ftRows = nRows * 2;
  let tempTransform = new Array(ftRows * ftCols);
  FFT.init(nCols);
  // transform rows
  let tmpRows = { re: new Array(nCols), im: new Array(nCols) };
  let row1 = { re: new Array(nCols), im: new Array(nCols) };
  let row2 = { re: new Array(nCols), im: new Array(nCols) };
  let index, iRow0, iRow1, iRow2, iRow3;
  for (iRow = 0; iRow < nRows / 2; iRow++) {
    index = iRow * 2 * nCols;
    tmpRows.re = data.slice(index, index + nCols);

    index = (iRow * 2 + 1) * nCols;
    tmpRows.im = data.slice(index, index + nCols);

    FFT.fft1d(tmpRows.re, tmpRows.im);

    reconstructTwoRealFFT(tmpRows, row1, row2);
    //Now lets put back the result into the output array
    iRow0 = iRow * 4 * ftCols;
    iRow1 = (iRow * 4 + 1) * ftCols;
    iRow2 = (iRow * 4 + 2) * ftCols;
    iRow3 = (iRow * 4 + 3) * ftCols;
    for (let k = ftCols - 1; k >= 0; k--) {
      tempTransform[iRow0 + k] = row1.re[k];
      tempTransform[iRow1 + k] = row1.im[k];
      tempTransform[iRow2 + k] = row2.re[k];
      tempTransform[iRow3 + k] = row2.im[k];
    }
  }

  //console.log(tempTransform);
  row1 = null;
  row2 = null;
  // transform columns
  let finalTransform = new Array(ftRows * ftCols);

  FFT.init(nRows);
  let tmpCols = { re: new Array(nRows), im: new Array(nRows) };
  for (iCol = ftCols - 1; iCol >= 0; iCol--) {
    for (iRow = nRows - 1; iRow >= 0; iRow--) {
      tmpCols.re[iRow] = tempTransform[iRow * 2 * ftCols + iCol];
      tmpCols.im[iRow] = tempTransform[(iRow * 2 + 1) * ftCols + iCol];
      //TODO Check why this happens
      if (isNaN(tmpCols.re[iRow])) {
        tmpCols.re[iRow] = 0;
      }
      if (isNaN(tmpCols.im[iRow])) {
        tmpCols.im[iRow] = 0;
      }
    }
    FFT.fft1d(tmpCols.re, tmpCols.im);
    for (iRow = nRows - 1; iRow >= 0; iRow--) {
      finalTransform[iRow * 2 * ftCols + iCol] = tmpCols.re[iRow];
      finalTransform[(iRow * 2 + 1) * ftCols + iCol] = tmpCols.im[iRow];
    }
  }
  return finalTransform;
}

/**
 *
 * @param fourierTransform
 * @param realTransform1
 * @param realTransform2
 *
 * Reconstructs the individual Fourier transforms of two simultaneously
 * transformed series. Based on the Symmetry relationships (the asterisk
 * denotes the complex conjugate)
 *
 * F_{N-n} = F_n^{*} for a purely real f transformed to F
 *
 * G_{N-n} = G_n^{*} for a purely imaginary g transformed to G
 *
 */
function reconstructTwoRealFFT(
  fourierTransform,
  realTransform1,
  realTransform2,
) {
  let length = fourierTransform.re.length;

  // the components n=0 are trivial
  realTransform1.re[0] = fourierTransform.re[0];
  realTransform1.im[0] = 0.0;
  realTransform2.re[0] = fourierTransform.im[0];
  realTransform2.im[0] = 0.0;
  let rm, rp, im, ip, j;
  for (let i = length / 2; i > 0; i--) {
    j = length - i;
    rm = 0.5 * (fourierTransform.re[i] - fourierTransform.re[j]);
    rp = 0.5 * (fourierTransform.re[i] + fourierTransform.re[j]);
    im = 0.5 * (fourierTransform.im[i] - fourierTransform.im[j]);
    ip = 0.5 * (fourierTransform.im[i] + fourierTransform.im[j]);
    realTransform1.re[i] = rp;
    realTransform1.im[i] = im;
    realTransform1.re[j] = rp;
    realTransform1.im[j] = -im;
    realTransform2.re[i] = ip;
    realTransform2.im[i] = -rm;
    realTransform2.re[j] = ip;
    realTransform2.im[j] = rm;
  }
}

/**
 * In place version of convolute 2D
 *
 * @param ftSignal
 * @param ftFilter
 * @param ftRows
 * @param ftCols
 * @return
 */
function convolute2DI(ftSignal, ftFilter, ftRows, ftCols) {
  let re, im;
  for (let iRow = 0; iRow < ftRows / 2; iRow++) {
    for (let iCol = 0; iCol < ftCols; iCol++) {
      //
      re =
        ftSignal[iRow * 2 * ftCols + iCol] *
          ftFilter[iRow * 2 * ftCols + iCol] -
        ftSignal[(iRow * 2 + 1) * ftCols + iCol] *
          ftFilter[(iRow * 2 + 1) * ftCols + iCol];
      im =
        ftSignal[iRow * 2 * ftCols + iCol] *
          ftFilter[(iRow * 2 + 1) * ftCols + iCol] +
        ftSignal[(iRow * 2 + 1) * ftCols + iCol] *
          ftFilter[iRow * 2 * ftCols + iCol];
      //
      ftSignal[iRow * 2 * ftCols + iCol] = re;
      ftSignal[(iRow * 2 + 1) * ftCols + iCol] = im;
    }
  }
}

/**
 *
 * @param data
 * @param kernel
 * @param nRows
 * @param nCols
 * @returns {*}
 */
function convolute(data, kernel, nRows, nCols) {
  let i;
  let ftSpectrum = new Array(nCols * nRows);
  for (i = 0; i < nRows * nCols; i++) {
    ftSpectrum[i] = data[i];
  }

  ftSpectrum = fft2DArray(ftSpectrum, nRows, nCols);

  let dimR = kernel.length;
  let dimC = kernel[0].length;
  let ftFilterData = new Array(nCols * nRows);
  for (i = 0; i < nCols * nRows; i++) {
    ftFilterData[i] = 0;
  }

  let iRow, iCol;
  let shiftR = Math.floor((dimR - 1) / 2);
  let shiftC = Math.floor((dimC - 1) / 2);
  for (let ir = 0; ir < dimR; ir++) {
    iRow = (ir - shiftR + nRows) % nRows;
    for (let ic = 0; ic < dimC; ic++) {
      iCol = (ic - shiftC + nCols) % nCols;
      ftFilterData[iRow * nCols + iCol] = kernel[ir][ic];
    }
  }
  ftFilterData = fft2DArray(ftFilterData, nRows, nCols);

  let ftRows = nRows * 2;
  let ftCols = nCols / 2 + 1;
  convolute2DI(ftSpectrum, ftFilterData, ftRows, ftCols);

  return ifft2DArray(ftSpectrum, ftRows, ftCols);
}

function toRadix2(data, nRows, nCols) {
  let i, j, irow, icol;
  let cols = nCols;
  let rows = nRows;
  //var prows = 0;
  //var pcols = 0;
  if (!(nCols !== 0 && (nCols & (nCols - 1)) === 0)) {
    //Then we have to make a pading to next radix2
    cols = 0;
    while (nCols >> ++cols !== 0);
    cols = 1 << cols;
    //pcols = cols - nCols;
  }
  if (!(nRows !== 0 && (nRows & (nRows - 1)) === 0)) {
    //Then we have to make a pading to next radix2
    rows = 0;
    while (nRows >> ++rows !== 0);
    rows = 1 << rows;
    //prows = (rows - nRows) * cols;
  }
  if (rows === nRows && cols === nCols) {
    //Do nothing. Returns the same input!!! Be careful
    return { data: data, rows: nRows, cols: nCols };
  }

  let output = new Array(rows * cols);
  let shiftR = Math.floor((rows - nRows) / 2) - nRows;
  let shiftC = Math.floor((cols - nCols) / 2) - nCols;

  for (i = 0; i < rows; i++) {
    irow = i * cols;
    icol = ((i - shiftR) % nRows) * nCols;
    for (j = 0; j < cols; j++) {
      output[irow + j] = data[icol + ((j - shiftC) % nCols)];
    }
  }
  return { data: output, rows: rows, cols: cols };
}

/**
 * Crop the given matrix to fit the corresponding number of rows and columns
 */
function crop(data, rows, cols, nRows, nCols) {
  if (rows === nRows && cols === nCols) {
    //Do nothing. Returns the same input!!! Be careful
    return data;
  }

  let output = new Array(nCols * nRows);

  let shiftR = Math.floor((rows - nRows) / 2);
  let shiftC = Math.floor((cols - nCols) / 2);
  let destinyRow, sourceRow, i, j;
  for (i = 0; i < nRows; i++) {
    destinyRow = i * nCols;
    sourceRow = (i + shiftR) * cols;
    for (j = 0; j < nCols; j++) {
      output[destinyRow + j] = data[sourceRow + (j + shiftC)];
    }
  }

  return output;
}

export const FFTUtils = {
  crop,
  convolute,
  toRadix2,
  convolute2DI,
  reconstructTwoRealFFT,
  fft2DArray,
  ifft2DArray,
};