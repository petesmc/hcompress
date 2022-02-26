use std::io::{BufRead, Cursor, Read};

use bytes::Buf;
//use bytes::{BytesMut, BufMut, Buf, Bytes};

fn ffpmsg(_m: &str) {}

/*  #########################################################################
These routines to apply the H-compress decompression algorithm to a 2-D Fits
image were written by R. White at the STScI and were obtained from the STScI at
http://www.stsci.edu/software/hcompress.html

This source file is a concatination of the following sources files in the
original distribution
  hinv.c
  hsmooth.c
  undigitize.c
  decode.c
  dodecode.c
  qtree_decode.c
  qread.c
  bit_input.c


The following modifications have been made to the original code:

  - commented out redundant "include" statements
  - added the nextchar global variable
  - changed all the 'extern' declarations to 'static', since all the routines are in
    the same source file
  - changed the first parameter in decode (and in lower level routines from a file stream
    to a char array
  - modified the myread routine, and lower level byte reading routines,  to copy
    the input bytes to a char array, instead of reading them from a file stream
  - changed the function declarations to the more modern ANSI C style
  - changed calls to printf and perror to call the CFITSIO ffpmsg routine
  - replace "exit" statements with "return" statements

 ############################################################################  */

// static long nextchar;
#[must_use]
pub fn min(a: i32, b: i32) -> i32 {
    if a < b {
        a
    } else {
        b
    }
}

#[must_use]
pub fn max(a: i32, b: i32) -> i32 {
    if a > b {
        a
    } else {
        b
    }
}

#[must_use]
pub fn min64(a: i64, b: i64) -> i64 {
    if a < b {
        a
    } else {
        b
    }
}

#[must_use]
pub fn max64(a: i64, b: i64) -> i64 {
    if a > b {
        a
    } else {
        b
    }
}

/* ---------------------------------------------------------------------- */
pub fn fits_hdecompress(inputz: &[u8], smooth: i32, a: &mut [i32]) -> (i32, usize, usize, i32) {
    /*
        decompress the input byte stream using the H-compress algorithm

      input  - input array of compressed bytes
      a - pre-allocated array to hold the output uncompressed image
      nx - returned X axis size
      ny - returned Y axis size

    NOTE: the nx and ny dimensions as defined within this code are reversed from
    the usual FITS notation.  ny is the fastest varying dimension, which is
    usually considered the X axis in the FITS image display

     */

    let mut input = Cursor::new(inputz);

    /* decode the input array */
    let res = decode(&mut input, a);

    if res.0 != 0 {
        return (res.0, res.1, res.2, res.3)
    }

    let nx = res.1;
    let ny = res.2;
    let scale = res.3;
    /*
     * Un-Digitize
     */
    undigitize(a, nx, ny, scale);

    /*
     * Inverse H-transform
     */
    let stat = hinv(a, nx, ny, smooth, scale);

    (stat, nx, ny, scale)
}

/*  ############################################################################  */
/*  ############################################################################  */

/* Copyright (c) 1993 Association of Universities for Research
 * in Astronomy. All rights reserved. Produced under National
 * Aeronautics and Space Administration Contract No. NAS5-26555.
 */
/* hinv.c   Inverse H-transform of NX x NY integer image
 *
 * Programmer: R. White		Date: 23 July 1993
 */

/*  ############################################################################  */
pub fn hinv(a: &mut [i32], nx: usize, ny: usize, smooth: i32, scale: i32) -> i32
/*
 int smooth;    0 for no smoothing, else smooth during inversion 
 int scale;     used if smoothing is specified 
 */ {
    // int nmax, log2n, i, j, k;
    // int nxtop,nytop,nxf,nyf,c;
    // int oddx,oddy;
    // int shift, bit0, bit1, bit2, mask0, mask1, mask2,
    //     prnd0, prnd1, prnd2, nrnd0, nrnd1, nrnd2, lowbit0, lowbit1;
    let mut lowbit0;
    let mut lowbit1;
    // int h0, hx, hy, hc;
    // int s10, s00;
    // int *tmp;

    let mut h0: i32;
    let mut hx: i32;
    let mut hy: i32;
    let mut hc: i32;

    let mut oddx: usize;
    let mut oddy: usize;

    let mut s10: usize;
    let mut s00: usize;

    /*
     * log2n is log2 of max(nx,ny) rounded up to next power of 2
     */
    let nmax: usize = if nx > ny { nx } else { ny };
    let mut log2n: usize = ((nmax as f64).ln() / 2.0_f64.ln() + 0.5) as usize;

    if nmax > (1 << log2n) {
        log2n += 1;
    }
    /*
     * get temporary storage for shuffling elements
     */
    let mut tmp: Vec<i32> = vec![0; (nmax + 1) / 2];

    /*
     * set up masks, rounding parameters
     */
    let mut shift = 1;
    let mut bit0 = 1 << (log2n - 1);
    let mut bit1 = bit0 << 1;
    let mut bit2 = bit0 << 2;
    let mut mask0 = -bit0;
    let mut mask1 = mask0 << 1;
    let mask2 = mask0 << 2;
    let mut prnd0 = bit0 >> 1;
    let mut prnd1 = bit1 >> 1;
    let prnd2 = bit2 >> 1;
    let mut nrnd0 = prnd0 - 1;
    let mut nrnd1 = prnd1 - 1;
    let nrnd2 = prnd2 - 1;
    /*
     * round h0 to multiple of bit2
     */
    a[0] = (a[0] + (if a[0] >= 0 { prnd2 } else { nrnd2 })) & mask2;
    /*
     * do log2n expansions
     *
     * We're indexing a as a 2-D array with dimensions (nx,ny).
     */
    let mut nxtop = 1;
    let mut nytop = 1;
    let mut nxf = nx;
    let mut nyf = ny;
    let mut c = 1 << log2n;
    for k in (0..log2n).rev() {
        //for (k = log2n-1; k>=0; k--) {
        /*
         * this somewhat cryptic code generates the sequence
         * ntop[k-1] = (ntop[k]+1)/2, where ntop[log2n] = n
         */
        c >>= 1;
        nxtop <<= 1;
        nytop <<= 1;
        if nxf <= c {
            nxtop -= 1;
        } else {
            nxf -= c;
        }
        if nyf <= c {
            nytop -= 1;
        } else {
            nyf -= c;
        }
        /*
         * double shift and fix nrnd0 (because prnd0=0) on last pass
         */
        if k == 0 {
            nrnd0 = 0;
            shift = 2;
        }
        /*
         * unshuffle in each dimension to interleave coefficients
         */
        for i in 0..nxtop {
            //for (i = 0; i<nxtop; i++) {
            unshuffle(&mut a[ny * i..], nytop, 1, &mut tmp);
        }

        for j in 0..nytop {
            // for (j = 0; j<nytop; j++) {
            unshuffle(&mut a[j..], nxtop, ny, &mut tmp);
        }
        /*
         * smooth by interpolating coefficients if SMOOTH != 0
         */
        if smooth > 0 {
            hsmooth(a, nxtop, nytop, ny, scale);
        }

        oddx = nxtop % 2;
        oddy = nytop % 2;

        for i in (0..(nxtop - oddx)).step_by(2) {
            // for (i = 0; i<nxtop-oddx; i += 2) {
            s00 = ny * i; /* s00 is index of a[i,j]	*/
            s10 = s00 + ny; /* s10 is index of a[i+1,j]	*/
            for _j in (0..(nytop - oddy)).step_by(2) {
                // for (j = 0; j<nytop-oddy; j += 2) {
                h0 = a[s00];
                hx = a[s10];
                hy = a[s00 + 1];
                hc = a[s10 + 1];
                /*
                 * round hx and hy to multiple of bit1, hc to multiple of bit0
                 * h0 is already a multiple of bit2
                 */
                hx = (hx + (if hx >= 0 { prnd1 } else { nrnd1 })) & mask1;
                hy = (hy + (if hy >= 0 { prnd1 } else { nrnd1 })) & mask1;
                hc = (hc + (if hc >= 0 { prnd0 } else { nrnd0 })) & mask0;
                /*
                 * propagate bit0 of hc to hx,hy
                 */
                lowbit0 = hc & bit0;
                hx = if hx >= 0 { hx - lowbit0 } else { hx + lowbit0 };
                hy = if hy >= 0 { hy - lowbit0 } else { hy + lowbit0 };
                /*
                 * Propagate bits 0 and 1 of hc,hx,hy to h0.
                 * This could be simplified if we assume h0>0, but then
                 * the inversion would not be lossless for images with
                 * negative pixels.
                 */
                lowbit1 = (hc ^ hx ^ hy) & bit1;
                h0 = if h0 >= 0 {
                    h0 + lowbit0 - lowbit1
                } else {
                    h0 + (if lowbit0 == 0 {
                        lowbit1
                    } else {
                        lowbit0 - lowbit1
                    })
                };
                /*
                 * Divide sums by 2 (4 last time)
                 */
                a[s10 + 1] = (h0 + hx + hy + hc) >> shift;
                a[s10] = (h0 + hx - hy - hc) >> shift;
                a[s00 + 1] = (h0 - hx + hy - hc) >> shift;
                a[s00] = (h0 - hx - hy + hc) >> shift;
                s00 += 2;
                s10 += 2;
            }
            if oddy > 0 {
                /*
                 * do last element in row if row length is odd
                 * s00+1, s10+1 are off edge
                 */
                h0 = a[s00];
                hx = a[s10];
                hx = (if hx >= 0 { hx + prnd1 } else { hx + nrnd1 }) & mask1;
                lowbit1 = hx & bit1;
                h0 = if h0 >= 0 { h0 - lowbit1 } else { h0 + lowbit1 };
                a[s10] = (h0 + hx) >> shift;
                a[s00] = (h0 - hx) >> shift;
            }
        }
        if oddx > 0 {
            /*
             * do last row if column length is odd
             * s10, s10+1 are off edge
             */
            s00 = ny * (nxtop-1); //TODO doublec heck
            for _j in (0..(nytop - oddy)).step_by(2) {
                //for (j = 0; j<nytop-oddy; j += 2) {
                h0 = a[s00];
                hy = a[s00 + 1];
                hy = (if hy >= 0 { hy + prnd1 } else { hy + nrnd1 }) & mask1;
                lowbit1 = hy & bit1;
                h0 = if h0 >= 0 { h0 - lowbit1 } else { h0 + lowbit1 };
                a[s00 + 1] = (h0 + hy) >> shift;
                a[s00] = (h0 - hy) >> shift;
                s00 += 2;
            }
            if oddy > 0 {
                /*
                 * do corner element if both row and column lengths are odd
                 * s00+1, s10, s10+1 are off edge
                 */
                h0 = a[s00];
                a[s00] = h0 >> shift;
            }
        }
        /*
         * divide all the masks and rounding values by 2
         */
        bit2 = bit1;
        bit1 = bit0;
        bit0 >>= 1;
        mask1 = mask0;
        mask0 >>= 1;
        prnd1 = prnd0;
        prnd0 >>= 1;
        nrnd1 = nrnd0;
        nrnd0 = prnd0 - 1;
    }

    0
}
/*  ############################################################################  */

/*  ############################################################################  */
pub fn unshuffle(a: &mut [i32], n: usize, n2: usize, tmp: &mut [i32])
/*
int a[];	 array to shuffle
int n;		 number of elements to shuffle
int n2;		 second dimension
int tmp[];	 scratch storage
*/
{
    // int i;
    // int nhalf;
    // int *p1, *p2, *pt;

    /*
     * copy 2nd half of array to tmp
     */
    let nhalf = (n + 1) >> 1;
    let mut pt = 0;
    let mut p1 = n2 * nhalf; /* pointer to a[i]			*/
    for _i in nhalf..n {
        //for (i=nhalf; i<n; i++) {
        tmp[pt] = a[p1];
        p1 += n2;
        pt += 1;
    }
    /*
     * distribute 1st half of array to even elements
     */
    let mut p2 = n2 * (nhalf - 1); /* pointer to a[i]			*/
    p1 = (n2 * (nhalf - 1)) << 1; /* pointer to a[2*i]		*/
    for i in (0..nhalf).rev() {
        //for (i=nhalf-1; i >= 0; i--) {
        a[p1] = a[p2];
        if i == 0 {
            break;
        }
        p2 -= n2;
        p1 -= n2 + n2;
    }

    /*
    n2 = 1
    nhalf = 1
    p2 = 0
    p1 = 0
    i = 0
    a[0] = a[0]
    p2 = 0 - 1   -> !!!!! overflow
    p1 = 0 - (1+1) = -2 -> !!!! overflow

    */

    /*
     * now distribute 2nd half of array (in tmp) to odd elements
     */
    pt = 0;
    p1 = n2; /* pointer to a[i]			*/
    for _i in (1..n).step_by(2) {
        //for (i=1; i<n; i += 2) {
        a[p1] = tmp[pt];
        p1 += n2 + n2;
        pt += 1;
    }
}
/*  ############################################################################  */
pub fn unshuffle64(a: &mut [i64], n: usize, n2: usize, tmp: &mut [i64])
/*
LONGLONG a[];	 array to shuffle
int n;		 number of elements to shuffle
int n2;		 second dimension
LONGLONG tmp[];	 scratch storage
*/
{
    /*
     * copy 2nd half of array to tmp
     */
    let nhalf = (n + 1) >> 1;
    let mut pt = 0;
    let mut p1 = n2 * nhalf; /* pointer to a[i]			*/
    for _i in nhalf..n {
        //for (i=nhalf; i<n; i++) {
        tmp[pt] = a[p1];
        p1 += n2;
        pt += 1;
    }
    /*
     * distribute 1st half of array to even elements
     */
    let mut p2 = n2 * (nhalf - 1); /* pointer to a[i]			*/
    p1 = (n2 * (nhalf - 1)) << 1; /* pointer to a[2*i]		*/
    for _i in (0..nhalf).rev() {
        //for (i=nhalf-1; i >= 0; i--) {
        a[p1] = a[p2];
        p2 -= n2;
        p1 -= n2 + n2;
    }
    /*
     * now distribute 2nd half of array (in tmp) to odd elements
     */
    pt = 0;
    p1 = n2; /* pointer to a[i]			*/
    for _i in (1..n).step_by(2) {
        //for (i=1; i<n; i += 2) {
        a[p1] = tmp[pt];
        p1 += n2 + n2;
        pt += 1;
    }
}

/*  ############################################################################  */
/*  ############################################################################  */

/* Copyright (c) 1993 Association of Universities for Research
 * in Astronomy. All rights reserved. Produced under National
 * Aeronautics and Space Administration Contract No. NAS5-26555.
 */
/* hsmooth.c	Smooth H-transform image by adjusting coefficients toward
 *				interpolated values
 *
 * Programmer: R. White		Date: 13 April 1992
 */

/*  ############################################################################  */
pub fn hsmooth(a: &mut [i32], nxtop: usize, nytop: usize, ny: usize, scale: i32)
/*
int a[];			 array of H-transform coefficients
int nxtop,nytop;	 size of coefficient block to use
int ny;				 actual 1st dimension of array
int scale;			 truncation scale factor that was used
*/
{
    // int i, j;
    // int ny2, s10, s00, diff, dmax, dmin, s, smax;
    // int hm, h0, hp, hmm, hpm, hmp, hpp, hx2, hy2;
    // int m1,m2;
    let mut m1;
    let mut m2;
    let mut hm;
    let mut h0;
    let mut hp;
    let mut hmm;
    let mut hpm;
    let mut hmp;
    let mut hpp;
    let mut hx2;
    let mut hy2;
    let mut diff;
    let mut dmax;
    let mut dmin;
    let mut s;
    let mut s00: usize;
    let mut s10: usize;

    /*
     * Maximum change in coefficients is determined by scale factor.
     * Since we rounded during division (see digitize.c), the biggest
     * permitted change is scale/2.
     */
    let smax = scale >> 1;
    if smax <= 0 {
        return;
    }
    let ny2 = ny << 1;
    /*
     * We're indexing a as a 2-D array with dimensions (nxtop,ny) of which
     * only (nxtop,nytop) are used.  The coefficients on the edge of the
     * array are not adjusted (which is why the loops below start at 2
     * instead of 0 and end at nxtop-2 instead of nxtop.)
     */
    /*
     * Adjust x difference hx
     */
    for i in (2..(nxtop - 2)).step_by(2) {
        // for (i = 2; i<nxtop-2; i += 2) {
        s00 = ny * i; /* s00 is index of a[i,j]	*/
        s10 = s00 + ny; /* s10 is index of a[i+1,j]	*/

        for _j in (0..nytop).step_by(2) {
            // for (j = 0; j<nytop; j += 2) {
            /*
             * hp is h0 (mean value) in next x zone, hm is h0 in previous x zone
             */
            hm = a[s00 - ny2];
            h0 = a[s00];
            hp = a[s00 + ny2];
            /*
             * diff = 8 * hx slope that would match h0 in neighboring zones
             */
            diff = hp - hm;
            /*
             * monotonicity constraints on diff
             */
            dmax = max(min(hp - h0, h0 - hm), 0) << 2;
            dmin = min(max(hp - h0, h0 - hm), 0) << 2;
            /*
             * if monotonicity would set slope = 0 then don't change hx.
             * note dmax>=0, dmin<=0.
             */
            if dmin < dmax {
                diff = max(min(diff, dmax), dmin);
                /*
                 * Compute change in slope limited to range +/- smax.
                 * Careful with rounding negative numbers when using
                 * shift for divide by 8.
                 */
                s = diff - (a[s10] << 3);
                s = if s >= 0 { s >> 3 } else { (s + 7) >> 3 };
                s = max(min(s, smax), -smax);
                a[s10] += s;
            }
            s00 += 2;
            s10 += 2;
        }
    }
    /*
     * Adjust y difference hy
     */
    for i in (0..nxtop).step_by(2) {
        // for (i = 0; i<nxtop; i += 2) {
        s00 = ny * i + 2;
        s10 = s00 + ny;
        for _j in (2..(nytop - 2)).step_by(2) {
            // for (j = 2; j<nytop-2; j += 2) {
            hm = a[s00 - 2];
            h0 = a[s00];
            hp = a[s00 + 2];
            diff = hp - hm;
            dmax = max(min(hp - h0, h0 - hm), 0) << 2;
            dmin = min(max(hp - h0, h0 - hm), 0) << 2;
            if dmin < dmax {
                diff = max(min(diff, dmax), dmin);
                s = diff - (a[s00 + 1] << 3);
                s = if s >= 0 { s >> 3 } else { (s + 7) >> 3 };
                s = max(min(s, smax), -smax);
                a[s00 + 1] += s;
            }
            s00 += 2;
            s10 += 2;
        }
    }
    /*
     * Adjust curvature difference hc
     */
    for i in (2..(nxtop - 2)).step_by(2) {
        //for (i = 2; i<nxtop-2; i += 2) {
        s00 = ny * i + 2;
        s10 = s00 + ny;
        for _j in (2..(nytop - 2)).step_by(2) {
            //for (j = 2; j<nytop-2; j += 2) {
            /*
             * ------------------    y
             * | hmp |    | hpp |    |
             * ------------------    |
             * |     | h0 |     |    |
             * ------------------    -------x
             * | hmm |    | hpm |
             * ------------------
             */
            hmm = a[s00 - ny2 - 2];
            hpm = a[s00 + ny2 - 2];
            hmp = a[s00 - ny2 + 2];
            hpp = a[s00 + ny2 + 2];
            h0 = a[s00];
            /*
             * diff = 64 * hc value that would match h0 in neighboring zones
             */
            diff = hpp + hmm - hmp - hpm;
            /*
             * 2 times x,y slopes in this zone
             */
            hx2 = a[s10] << 1;
            hy2 = a[s00 + 1] << 1;
            /*
             * monotonicity constraints on diff
             */
            m1 = min(max(hpp - h0, 0) - hx2 - hy2, max(h0 - hpm, 0) + hx2 - hy2);
            m2 = min(max(h0 - hmp, 0) - hx2 + hy2, max(hmm - h0, 0) + hx2 + hy2);
            dmax = min(m1, m2) << 4;
            m1 = max(min(hpp - h0, 0) - hx2 - hy2, min(h0 - hpm, 0) + hx2 - hy2);
            m2 = max(min(h0 - hmp, 0) - hx2 + hy2, min(hmm - h0, 0) + hx2 + hy2);
            dmin = max(m1, m2) << 4;
            /*
             * if monotonicity would set slope = 0 then don't change hc.
             * note dmax>=0, dmin<=0.
             */
            if dmin < dmax {
                diff = max(min(diff, dmax), dmin);
                /*
                 * Compute change in slope limited to range +/- smax.
                 * Careful with rounding negative numbers when using
                 * shift for divide by 64.
                 */
                s = diff - (a[s10 + 1] << 6);
                s = if s >= 0 { s >> 6 } else { (s + 63) >> 6 };
                s = max(min(s, smax), -smax);
                a[s10 + 1] += s;
            }
            s00 += 2;
            s10 += 2;
        }
    }
}
/*  ############################################################################  */
pub fn hsmooth64(a: &mut [i64], nxtop: usize, nytop: usize, ny: usize, scale: i32)
/*
LONGLONG a[];			 array of H-transform coefficients
int nxtop,nytop;	 size of coefficient block to use
int ny;				 actual 1st dimension of array
int scale;			 truncation scale factor that was used
*/
{
    //  int i, j;
    //  int ny2, s10, s00;
    //  LONGLONG hm, h0, hp, hmm, hpm, hmp, hpp, hx2, hy2, diff, dmax, dmin, s, smax, m1, m2;

    /*
     * Maximum change in coefficients is determined by scale factor.
     * Since we rounded during division (see digitize.c), the biggest
     * permitted change is scale/2.
     */
    let mut hm;
    let mut h0;
    let mut hp;
    let mut hmm;
    let mut hpm;
    let mut hmp;
    let mut hpp;
    let mut hx2;
    let mut hy2;
    let mut diff;
    let mut dmax;
    let mut dmin;
    let mut s;
    let mut s00: usize;
    let mut s10: usize;
    let mut m1;
    let mut m2;

    let smax = i64::from(scale >> 1);
    if smax <= 0 {
        return;
    }
    let ny2 = ny << 1;
    /*
     * We're indexing a as a 2-D array with dimensions (nxtop,ny) of which
     * only (nxtop,nytop) are used.  The coefficients on the edge of the
     * array are not adjusted (which is why the loops below start at 2
     * instead of 0 and end at nxtop-2 instead of nxtop.)
     */
    /*
     * Adjust x difference hx
     */
    for i in (2..(nxtop - 2)).step_by(2) {
        // for (i = 2; i<nxtop-2; i += 2) {
        s00 = ny * i; /* s00 is index of a[i,j]	*/
        s10 = s00 + ny; /* s10 is index of a[i+1,j]	*/

        for _j in (0..nytop).step_by(2) {
            // for (j = 0; j<nytop; j += 2) {
            /*
             * hp is h0 (mean value) in next x zone, hm is h0 in previous x zone
             */
            hm = a[s00 - ny2];
            h0 = a[s00];
            hp = a[s00 + ny2];
            /*
             * diff = 8 * hx slope that would match h0 in neighboring zones
             */
            diff = hp - hm;
            /*
             * monotonicity constraints on diff
             */
            dmax = max64(min64(hp - h0, h0 - hm), 0) << 2;
            dmin = min64(max64(hp - h0, h0 - hm), 0) << 2;
            /*
             * if monotonicity would set slope = 0 then don't change hx.
             * note dmax>=0, dmin<=0.
             */
            if dmin < dmax {
                diff = max64(min64(diff, dmax), dmin);
                /*
                 * Compute change in slope limited to range +/- smax.
                 * Careful with rounding negative numbers when using
                 * shift for divide by 8.
                 */
                s = diff - (a[s10] << 3);
                s = if s >= 0 { s >> 3 } else { (s + 7) >> 3 };
                s = max64(min64(s, smax), -smax);
                a[s10] += s;
            }
            s00 += 2;
            s10 += 2;
        }
    }
    /*
     * Adjust y difference hy
     */
    for i in (0..nxtop).step_by(2) {
        // for (i = 0; i<nxtop; i += 2) {
        s00 = ny * i + 2;
        s10 = s00 + ny;
        for _j in (2..(nytop - 2)).step_by(2) {
            // for (j = 2; j<nytop-2; j += 2) {
            hm = a[s00 - 2];
            h0 = a[s00];
            hp = a[s00 + 2];
            diff = hp - hm;
            dmax = max64(min64(hp - h0, h0 - hm), 0) << 2;
            dmin = min64(max64(hp - h0, h0 - hm), 0) << 2;
            if dmin < dmax {
                diff = max64(min64(diff, dmax), dmin);
                s = diff - (a[s00 + 1] << 3);
                s = if s >= 0 { s >> 3 } else { (s + 7) >> 3 };
                s = max64(min64(s, smax), -smax);
                a[s00 + 1] += s;
            }
            s00 += 2;
            s10 += 2;
        }
    }
    /*
     * Adjust curvature difference hc
     */
    for i in (2..(nxtop - 2)).step_by(2) {
        //for (i = 2; i<nxtop-2; i += 2) {
        s00 = ny * i + 2;
        s10 = s00 + ny;
        for _j in (2..(nytop - 2)).step_by(2) {
            //for (j = 2; j<nytop-2; j += 2) {
            /*
             * ------------------    y
             * | hmp |    | hpp |    |
             * ------------------    |
             * |     | h0 |     |    |
             * ------------------    -------x
             * | hmm |    | hpm |
             * ------------------
             */
            hmm = a[s00 - ny2 - 2];
            hpm = a[s00 + ny2 - 2];
            hmp = a[s00 - ny2 + 2];
            hpp = a[s00 + ny2 + 2];
            h0 = a[s00];
            /*
             * diff = 64 * hc value that would match h0 in neighboring zones
             */
            diff = hpp + hmm - hmp - hpm;
            /*
             * 2 times x,y slopes in this zone
             */
            hx2 = a[s10] << 1;
            hy2 = a[s00 + 1] << 1;
            /*
             * monotonicity constraints on diff
             */
            m1 = min64(
                max64(hpp - h0, 0) - hx2 - hy2,
                max64(h0 - hpm, 0) + hx2 - hy2,
            );
            m2 = min64(
                max64(h0 - hmp, 0) - hx2 + hy2,
                max64(hmm - h0, 0) + hx2 + hy2,
            );
            dmax = min64(m1, m2) << 4;
            m1 = max64(
                min64(hpp - h0, 0) - hx2 - hy2,
                min64(h0 - hpm, 0) + hx2 - hy2,
            );
            m2 = max64(
                min64(h0 - hmp, 0) - hx2 + hy2,
                min64(hmm - h0, 0) + hx2 + hy2,
            );
            dmin = max64(m1, m2) << 4;
            /*
             * if monotonicity would set slope = 0 then don't change hc.
             * note dmax>=0, dmin<=0.
             */
            if dmin < dmax {
                diff = max64(min64(diff, dmax), dmin);
                /*
                 * Compute change in slope limited to range +/- smax.
                 * Careful with rounding negative numbers when using
                 * shift for divide by 64.
                 */
                s = diff - (a[s10 + 1] << 6);
                s = if s >= 0 { s >> 6 } else { (s + 63) >> 6 };
                s = max64(min64(s, smax), -smax);
                a[s10 + 1] += s;
            }
            s00 += 2;
            s10 += 2;
        }
    }
}

/*  ############################################################################  */
/*  ############################################################################  */
/* Copyright (c) 1993 Association of Universities for Research
 * in Astronomy. All rights reserved. Produced under National
 * Aeronautics and Space Administration Contract No. NAS5-26555.
 */
/* undigitize.c		undigitize H-transform
 *
 * Programmer: R. White		Date: 9 May 1991
 */

/*  ############################################################################  */
pub fn undigitize(a: &mut [i32], nx: usize, ny: usize, scale: i32) {
    let _p = 0;

    /*
     * multiply by scale
     */
    if scale <= 1 {
        return;
    }

    if nx == 0 || ny == 0 {
        return;
    }

    for item in a.iter_mut().take(nx * ny - 1) {
        *item *= scale;
    }
    //for (p=a; p <= &a[nx*ny-1]; p++) *p = (*p)*scale;
}
/*  ############################################################################  */
pub fn undigitize64(a: &mut [i64], nx: usize, ny: usize, scale: i32) {
    /*
     * multiply by scale
     */
    if scale <= 1 {
        return;
    }
    let scale64: i64 = i64::from(scale); /* use a 64-bit int for efficiency in the big loop */

    if nx == 0 || ny == 0 {
        return;
    }

    a.iter_mut().take(nx * ny - 1).for_each(|x| *x *= scale64);

    //for (p=a; p <= &a[nx*ny-1]; p++) *p = (*p)*scale64;
}

/*  ############################################################################  */
/*  ############################################################################  */
/* Copyright (c) 1993 Association of Universities for Research
 * in Astronomy. All rights reserved. Produced under National
 * Aeronautics and Space Administration Contract No. NAS5-26555.
 */
/* decode.c		read codes from infile and construct array
 *
 * Programmer: R. White		Date: 2 February 1994
 */

pub const CODE_MAGIC: [u8; 2] = [0xDD, 0x99];

/*  ############################################################################  */
pub fn decode(infile: &mut Cursor<&[u8]>, a: &mut [i32]) -> (i32, usize, usize, i32) {
    //nx: usize, ny: usize, scale: i32
    /*
    char *infile;				 input file
    int  *a;				 address of output array [nx][ny]
    int  *nx,*ny;				 size of output array
    int  *scale;				 scale factor for digitization
    */

    // LONGLONG sumall;
    // int stat;
    // unsigned char nbitplanes[3];
    let mut nbitplanes: [u8; 3] = [0; 3];
    let mut tmagic: [u8; 2] = [0; 2];
    /* initialize the byte read position to the beginning of the array */

    /*
     * File starts either with special 2-byte magic code or with
     * FITS keyword "SIMPLE  ="
     */
    qread(infile, &mut tmagic, 2);
    /*
     * check for correct magic code value
     */
    if !tmagic.eq(&CODE_MAGIC) {
        ffpmsg("bad file format");
        return (-1, 0, 0, 0); // TODO!!!!
    }

    let nx = readint(infile) as usize; /* x size of image			*/
    let ny = readint(infile) as usize; /* y size of image			*/
    let scale = readint(infile); /* scale factor for digitization	*/

    /* sum of all pixels	*/
    let sumall = readlonglong(infile);
    /* # bits in quadrants	*/

    let len = nbitplanes.len();
    qread(infile, &mut nbitplanes, len);

    let stat = dodecode(infile, a, nx, ny, nbitplanes);
    /*
     * put sum of all pixels back into pixel 0
     */
    a[0] = sumall as i32;
    (stat, nx, ny, scale)
}

/*  ############################################################################  */
/*  ############################################################################  */
/* Copyright (c) 1993 Association of Universities for Research
 * in Astronomy. All rights reserved. Produced under National
 * Aeronautics and Space Administration Contract No. NAS5-26555.
 */
/* dodecode.c	Decode stream of characters on infile and return array
 *
 * This version encodes the different quadrants separately
 *
 * Programmer: R. White		Date: 9 May 1991
 */

/*  ############################################################################  */
pub fn dodecode(
    infile: &mut Cursor<&[u8]>,
    a: &mut [i32],
    nx: usize,
    ny: usize,
    nbitplanes: [u8; 3],
) -> i32 {
    /* int a[];
       int nx,ny;					 Array dimensions are [nx][ny]
       unsigned char nbitplanes[3];		 Number of bit planes in quadrants
    */

    // int i, nel, nx2, ny2, stat;
    let nel = nx * ny;
    let nx2 = (nx + 1) / 2;
    let ny2 = (ny + 1) / 2;

    /*
     * initialize a to zero
     */
    a.iter_mut().take(nel).for_each(|x| *x = 0);

    /*
     * Initialize bit input
     */
    let mut b2 = start_inputing_bits();
    /*
     * read bit planes for each quadrant
     */
    let mut stat = qtree_decode(infile, a, ny, nx2, ny2, i32::from(nbitplanes[0]), &mut b2);
    if stat != 0 {
        return stat;
    }

    stat = qtree_decode(
        infile,
        &mut a[ny2..],
        ny,
        nx2,
        ny / 2,
        i32::from(nbitplanes[1]),
        &mut b2,
    );
    if stat != 0 {
        return stat;
    }

    stat = qtree_decode(
        infile,
        &mut a[(ny * nx2)..],
        ny,
        nx / 2,
        ny2,
        i32::from(nbitplanes[1]),
        &mut b2,
    );
    if stat != 0 {
        return stat;
    }

    stat = qtree_decode(
        infile,
        &mut a[(ny * nx2 + ny2)..],
        ny,
        nx / 2,
        ny / 2,
        i32::from(nbitplanes[2]),
        &mut b2,
    );
    if stat != 0 {
        return stat;
    }

    /*
     * make sure there is an EOF symbol (nybble=0) at end
     */
    if input_nybble(infile, &mut b2) != 0 {
        ffpmsg("dodecode: bad bit plane values");
        return -1;
    }
    /*
     * now get the sign bits
     * Re-initialize bit input
     */
    let mut b2 = start_inputing_bits();
    for item in a.iter_mut().take(nel) {
        //for (i=0; i<nel; i++) {
        if *item > 0 {
            /* tried putting the input_bit code in-line here, instead of */
            /* calling the function, but it made no difference in the speed */
            if input_bit(infile, &mut b2) > 0 {
                *item = -*item;
            }
        }
    }
    0
}
/*  ############################################################################  */

/*  ############################################################################  */
/*  ############################################################################  */
/* Copyright (c) 1993 Association of Universities for Research
 * in Astronomy. All rights reserved. Produced under National
 * Aeronautics and Space Administration Contract No. NAS5-26555.
 */
/* qtree_decode.c	Read stream of codes from infile and construct bit planes
 *					in quadrant of 2-D array using binary quadtree coding
 *
 * Programmer: R. White		Date: 7 May 1991
 */

/*  ############################################################################  */
pub fn qtree_decode(
    infile: &mut Cursor<&[u8]>,
    a: &mut [i32],
    n: usize,
    nqx: usize,
    nqy: usize,
    nbitplanes: i32,
    b2: &mut Buffer2,
) -> i32 {
    /*
    char *infile;
    int a[];				 a is 2-D array with dimensions (n,n)
    int n;					 length of full row in a
    int nqx;				 partial length of row to decode
    int nqy;				 partial length of column (<=n)
    int nbitplanes;				 number of bitplanes to decode
    */

    // int log2n, k, bit, b, nqmax;
    // int nx,ny,nfx,nfy,c;
    // int nqx2, nqy2;
    // unsigned char *scratch;
    let mut b;
    let mut nfx;
    let mut nfy;
    let mut c;

    let mut nx: usize;
    let mut ny: usize;

    /*
     * log2n is log2 of max(nqx,nqy) rounded up to next power of 2
     */
    let nqmax: usize = if nqx > nqy { nqx } else { nqy };
    let mut log2n: usize = ((nqmax as f32).ln() / 2.0_f32.ln() + 0.5) as usize;
    if nqmax > (1 << log2n) {
        log2n += 1;
    }
    /*
     * allocate scratch array for working space
     */
    let nqx2 = (nqx + 1) / 2;
    let nqy2 = (nqy + 1) / 2;

    let mut scratch: Vec<u8> = vec![0; nqx2 * nqy2+20];
    let mut scratch2: Vec<u8> = vec![0; nqx2 * nqy2+20];
    /*
     * now decode each bit plane, starting at the top
     * A is assumed to be initialized to zero
     */
    for bit in (0..nbitplanes).rev() {
        //for (bit = nbitplanes-1; bit >= 0; bit--) {
        /*
         * Was bitplane was quadtree-coded or written directly?
         */
        b = input_nybble(infile, b2);

        if b == 0 {
            /*
             * bit map was written directly
             */
            read_bdirect(infile, a, n, nqx, nqy, &mut scratch, bit, b2);
        } else if b != 0xf {
            ffpmsg("qtree_decode: bad format code");
            return -1;
        } else {
            /*
             * bitmap was quadtree-coded, do log2n expansions
             *
             * read first code
             */
            scratch[0] = input_huffman(infile, b2) as u8;
            /*
             * now do log2n expansions, reading codes from file as necessary
             */
            nx = 1;
            ny = 1;
            nfx = nqx;
            nfy = nqy;
            c = 1 << log2n;
            for _k in 1..log2n {
                //for (k = 1; k<log2n; k++) {
                /*
                 * this somewhat cryptic code generates the sequence
                 * n[k-1] = (n[k]+1)/2 where n[log2n]=nqx or nqy
                 */
                c >>= 1;
                nx <<= 1;
                ny <<= 1;
                if nfx <= c {
                    nx -= 1;
                } else {
                    nfx -= c;
                }
                if nfy <= c {
                    ny -= 1;
                } else {
                    nfy -= c;
                }
                qtree_expand(infile, &scratch, nx, ny, &mut scratch2, b2);
                scratch.copy_from_slice(&scratch2);
            }
            /*
             * now copy last set of 4-bit codes to bitplane bit of array a
             */
            qtree_bitins(&mut scratch, nqx, nqy, a, n, bit);
        }
    }
    0
}
/*  ############################################################################  */

/*  ############################################################################  */
/*
 * do one quadtree expansion step on array a[(nqx+1)/2,(nqy+1)/2]
 * results put into b[nqx,nqy] (which may be the same as a)
 */
pub fn qtree_expand(
    infile: &mut Cursor<&[u8]>,
    a: &[u8],
    nx: usize,
    ny: usize,
    b: &mut [u8],
    b2: &mut Buffer2,
) {
    /*
     * first copy a to b, expanding each 4-bit value
     */
    qtree_copy(a, nx, ny, b, ny);
    /*
     * now read new 4-bit values into b for each non-zero element
     */
    for i in (0..(nx * ny)).rev() {
        //for (i = nx*ny-1; i >= 0; i--) {
        if b[i] > 0 {
            b[i] = input_huffman(infile, b2) as u8;
        }
    }
}

/*  ############################################################################  */
/*
 * copy 4-bit values from a[(nx+1)/2,(ny+1)/2] to b[nx,ny], expanding
 * each value to 2x2 pixels
 * a,b may be same array
 */
pub fn qtree_copy(a: &mut [u8], nx: usize, ny: usize, n: usize)
/*   int n;		declared y dimension of b */
{
    // int i, j, k, nx2, ny2;
    // int s00, s10;

    let mut s00: usize;
    let mut s10: usize;
    /*
     * first copy 4-bit values to b
     * start at end in case a,b are same array
     */
    let nx2 = (nx + 1) / 2;
    let ny2 = (ny + 1) / 2;

    let mut k = ny2 * (nx2 - 1) + ny2;// was ny2 * (nx2 - 1) + ny2 -1 ;  /* k   is index of a[i,j]		*/
    for i in (0..nx2).rev() {
        //for (i = nx2-1; i >= 0; i--) {
            //TODO added +2 below to re-order the -2 to prevent underflow
        s00 = 2 * (n * i + ny2 - 1) + 2; /* s00 is index of b[2*i,2*j]		*/
        for _j in (0..ny2).rev() {
            //for (j = ny2-1; j >= 0; j--) {
            k -= 1;
            s00 -= 2;
            a[s00] = a[k];
            
        }
    }
    /*
     * now expand each 2x2 block
     */
    let oddx = nx % 2;
    let oddy = ny % 2;

    if nx == 0 {
        return;
    } // Return early and prevent underflow
    for i in (0..(nx - 1)).step_by(2) {
        //for (i = 0; i<nx-1; i += 2) {

        /* Note:
           Unlike the case in qtree_bitins, this code runs faster on a 32-bit linux
           machine using the s10 intermediate variable, rather that using s00+n.
           Go figure!
        */
        s00 = n * i; /* s00 is index of b[i,j]	*/
        s10 = s00 + n; /* s10 is index of b[i+1,j]	*/

        if ny == 0 {
            continue;
        } // Continue early and prevent underflow
        for _j in (0..(ny - 1)).step_by(2) {
            //for (j = 0; j<ny-1; j += 2) {

            match b[s00] {
                0 => {
                    b[s10 + 1] = 0;
                    b[s10] = 0;
                    b[s00 + 1] = 0;
                    b[s00] = 0;
                }

                1 => {
                    b[s10 + 1] = 1;
                    b[s10] = 0;
                    b[s00 + 1] = 0;
                    b[s00] = 0;
                }
                2 => {
                    b[s10 + 1] = 0;
                    b[s10] = 1;
                    b[s00 + 1] = 0;
                    b[s00] = 0;
                }
                3 => {
                    b[s10 + 1] = 1;
                    b[s10] = 1;
                    b[s00 + 1] = 0;
                    b[s00] = 0;
                }
                4 => {
                    b[s10 + 1] = 0;
                    b[s10] = 0;
                    b[s00 + 1] = 1;
                    b[s00] = 0;
                }
                5 => {
                    b[s10 + 1] = 1;
                    b[s10] = 0;
                    b[s00 + 1] = 1;
                    b[s00] = 0;
                }
                6 => {
                    b[s10 + 1] = 0;
                    b[s10] = 1;
                    b[s00 + 1] = 1;
                    b[s00] = 0;
                }
                7 => {
                    b[s10 + 1] = 1;
                    b[s10] = 1;
                    b[s00 + 1] = 1;
                    b[s00] = 0;
                }
                8 => {
                    b[s10 + 1] = 0;
                    b[s10] = 0;
                    b[s00 + 1] = 0;
                    b[s00] = 1;
                }
                9 => {
                    b[s10 + 1] = 1;
                    b[s10] = 0;
                    b[s00 + 1] = 0;
                    b[s00] = 1;
                }
                10 => {
                    b[s10 + 1] = 0;
                    b[s10] = 1;
                    b[s00 + 1] = 0;
                    b[s00] = 1;
                }
                11 => {
                    b[s10 + 1] = 1;
                    b[s10] = 1;
                    b[s00 + 1] = 0;
                    b[s00] = 1;
                }
                12 => {
                    b[s10 + 1] = 0;
                    b[s10] = 0;
                    b[s00 + 1] = 1;
                    b[s00] = 1;
                }
                13 => {
                    b[s10 + 1] = 1;
                    b[s10] = 0;
                    b[s00 + 1] = 1;
                    b[s00] = 1;
                }
                14 => {
                    b[s10 + 1] = 0;
                    b[s10] = 1;
                    b[s00 + 1] = 1;
                    b[s00] = 1;
                }
                15 => {
                    b[s10 + 1] = 1;
                    b[s10] = 1;
                    b[s00 + 1] = 1;
                    b[s00] = 1;
                }
                _ => (),
            }
            /*
                        b[s10+1] =  b[s00]     & 1;
                        b[s10  ] = (b[s00]>>1) & 1;
                        b[s00+1] = (b[s00]>>2) & 1;
                        b[s00  ] = (b[s00]>>3) & 1;
            */

            s00 += 2;
            s10 += 2;
        }

        if oddy > 0 {
            /*
             * row size is odd, do last element in row
             * s00+1, s10+1 are off edge
             */
            /* not worth converting this to use 16 case statements */
            b[s10] = (b[s00] >> 1) & 1;
            b[s00] = (b[s00] >> 3) & 1;
        }
    }
    if oddx > 0 {
        /*
         * column size is odd, do last row
         * s10, s10+1 are off edge
         */
        s00 = n * (nx - 1);
        if ny == 0 {
            return;
        } // Return early and prevent underflow
        for _j in (0..(ny - 1)).step_by(2) {
            //for (j = 0; j<ny-1; j += 2) {
            /* not worth converting this to use 16 case statements */
            b[s00 + 1] = (b[s00] >> 2) & 1;
            b[s00] = (b[s00] >> 3) & 1;
            s00 += 2;
        }
        if oddy > 0 {
            /*
             * both row and column size are odd, do corner element
             * s00+1, s10, s10+1 are off edge
             */
            /* not worth converting this to use 16 case statements */
            b[s00] = (b[s00] >> 3) & 1;
        }
    }
}

/*  ############################################################################  */
/*
 * Copy 4-bit values from a[(nx+1)/2,(ny+1)/2] to b[nx,ny], expanding
 * each value to 2x2 pixels and inserting into bitplane BIT of B.
 * A,B may NOT be same array (it wouldn't make sense to be inserting
 * bits into the same array anyway.)
 */
pub fn qtree_bitins(a: &mut [u8], nx: usize, ny: usize, b: &mut [i32], n: usize, bit: i32)
/*
   int n;		declared y dimension of b
*/
{
    // int i, j, k;

    let mut s00: usize;

    let plane_val = 1 << bit;

    /*
     * expand each 2x2 block
     */
    let mut k: usize = 0; /* k   is index of a[i/2,j/2]	*/
    let oddx = nx % 2;
    let oddy = nx % 2;

    if nx == 0 {
        return;
    } // Return early and prevent underflow
    for i in (0..(nx - 1)).step_by(2) {
        //for (i = 0; i<nx-1; i += 2) {
        s00 = n * i; /* s00 is index of b[i,j]	*/

        /* Note:
           this code appears to run very slightly faster on a 32-bit linux
           machine using s00+n rather than the s10 intermediate variable
        */
        /*		s10 = s00+n;	*/
        /* s10 is index of b[i+1,j]	*/
        if ny == 0 {
            continue;
        } // Continue early and prevent underflow
        for _j in (0..(ny - 1)).step_by(2) {
            //for (j = 0; j<ny-1; j += 2) {

            match a[k] {
                0 => (),
                1 => {
                    b[s00 + n + 1] |= plane_val;
                }
                2 => {
                    b[s00 + n] |= plane_val;
                }
                3 => {
                    b[s00 + n + 1] |= plane_val;
                    b[s00 + n] |= plane_val;
                }
                4 => {
                    b[s00 + 1] |= plane_val;
                }
                5 => {
                    b[s00 + n + 1] |= plane_val;
                    b[s00 + 1] |= plane_val;
                }
                6 => {
                    b[s00 + n] |= plane_val;
                    b[s00 + 1] |= plane_val;
                }
                7 => {
                    b[s00 + n + 1] |= plane_val;
                    b[s00 + n] |= plane_val;
                    b[s00 + 1] |= plane_val;
                }
                8 => {
                    b[s00] |= plane_val;
                }
                9 => {
                    b[s00 + n + 1] |= plane_val;
                    b[s00] |= plane_val;
                }
                10 => {
                    b[s00 + n] |= plane_val;
                    b[s00] |= plane_val;
                }
                11 => {
                    b[s00 + n + 1] |= plane_val;
                    b[s00 + n] |= plane_val;
                    b[s00] |= plane_val;
                }
                12 => {
                    b[s00 + 1] |= plane_val;
                    b[s00] |= plane_val;
                }
                13 => {
                    b[s00 + n + 1] |= plane_val;
                    b[s00 + 1] |= plane_val;
                    b[s00] |= plane_val;
                }
                14 => {
                    b[s00 + n] |= plane_val;
                    b[s00 + 1] |= plane_val;
                    b[s00] |= plane_val;
                }
                15 => {
                    b[s00 + n + 1] |= plane_val;
                    b[s00 + n] |= plane_val;
                    b[s00 + 1] |= plane_val;
                    b[s00] |= plane_val;
                }
                _ => (),
            }

            /*
                        b[s10+1] |= ( a[k]     & 1) << bit;
                        b[s10  ] |= ((a[k]>>1) & 1) << bit;
                        b[s00+1] |= ((a[k]>>2) & 1) << bit;
                        b[s00  ] |= ((a[k]>>3) & 1) << bit;
            */
            s00 += 2;
            /*			s10 += 2; */
            k += 1;
        }
        if oddy > 0 {
            /*
             * row size is odd, do last element in row
             * s00+1, s10+1 are off edge
             */

            match a[k] {
                2 => {
                    b[s00 + n] |= plane_val;
                }
                3 => {
                    b[s00 + n] |= plane_val;
                }
                6 => {
                    b[s00 + n] |= plane_val;
                }
                7 => {
                    b[s00 + n] |= plane_val;
                }
                8 => {
                    b[s00] |= plane_val;
                }
                9 => {
                    b[s00] |= plane_val;
                }
                10 => {
                    b[s00 + n] |= plane_val;
                    b[s00] |= plane_val;
                }
                11 => {
                    b[s00 + n] |= plane_val;
                    b[s00] |= plane_val;
                }
                12 => {
                    b[s00] |= plane_val;
                }
                13 => {
                    b[s00] |= plane_val;
                }
                14 => {
                    b[s00 + n] |= plane_val;
                    b[s00] |= plane_val;
                }
                15 => {
                    b[s00 + n] |= plane_val;
                    b[s00] |= plane_val;
                }
                _ => (),
            }

            /*
                        b[s10  ] |= ((a[k]>>1) & 1) << bit;
                        b[s00  ] |= ((a[k]>>3) & 1) << bit;
            */
            k += 1;
        }
    }
    if oddx > 0 {
        /*
         * column size is odd, do last row
         * s10, s10+1 are off edge
         */
        s00 = n * (nx - 1);
        if ny == 0 {
            return;
        } // Return early and prevent underflow
        for _j in (0..(ny - 1)).step_by(2) {
            //for (j = 0; j<ny-1; j += 2) {

            match a[k] {
                4 => {
                    b[s00 + 1] |= plane_val;
                }
                5 => {
                    b[s00 + 1] |= plane_val;
                }
                6 => {
                    b[s00 + 1] |= plane_val;
                }
                7 => {
                    b[s00 + 1] |= plane_val;
                }
                8 => {
                    b[s00] |= plane_val;
                }
                9 => {
                    b[s00] |= plane_val;
                }
                10 => {
                    b[s00] |= plane_val;
                }
                11 => {
                    b[s00] |= plane_val;
                }
                12 => {
                    b[s00 + 1] |= plane_val;
                    b[s00] |= plane_val;
                }
                13 => {
                    b[s00 + 1] |= plane_val;
                    b[s00] |= plane_val;
                }
                14 => {
                    b[s00 + 1] |= plane_val;
                    b[s00] |= plane_val;
                }
                15 => {
                    b[s00 + 1] |= plane_val;
                    b[s00] |= plane_val;
                }
                _ => (),
            }

            /*
                        b[s00+1] |= ((a[k]>>2) & 1) << bit;
                        b[s00  ] |= ((a[k]>>3) & 1) << bit;
            */

            s00 += 2;
            k += 1;
        }
        if oddy > 0 {
            /*
             * both row and column size are odd, do corner element
             * s00+1, s10, s10+1 are off edge
             */

            match a[k] {
                8 => {
                    b[s00] |= plane_val;
                }
                9 => {
                    b[s00] |= plane_val;
                }
                10 => {
                    b[s00] |= plane_val;
                }
                11 => {
                    b[s00] |= plane_val;
                }
                12 => {
                    b[s00] |= plane_val;
                }
                13 => {
                    b[s00] |= plane_val;
                }
                14 => {
                    b[s00] |= plane_val;
                }
                15 => {
                    b[s00] |= plane_val;
                }
                _ => (),
            }

            /*
                        b[s00  ] |= ((a[k]>>3) & 1) << bit;
            */
            k += 1;
        }
    }
}

/*  ############################################################################  */
pub fn read_bdirect(
    infile: &mut Cursor<&[u8]>,
    a: &mut [i32],
    n: usize,
    nqx: usize,
    nqy: usize,
    scratch: &mut [u8],
    bit: i32,
    b2: &mut Buffer2,
) {
    /*
     * read bit image packed 4 pixels/nybble
     */
    input_nnybble(infile, ((nqx + 1) / 2) * ((nqy + 1) / 2), scratch, b2);

    /*
     * insert in bitplane BIT of image A
     */
    qtree_bitins(scratch, nqx, nqy, a, n, bit);
}
/*  ############################################################################  */
pub fn read_bdirect64(
    infile: &mut Cursor<&[u8]>,
    _a: &[i64],
    _n: usize,
    nqx: usize,
    nqy: usize,
    scratch: &mut [u8],
    _bit: i32,
    b2: &mut Buffer2,
) {
    /*
     * read bit image packed 4 pixels/nybble
     */

    input_nnybble(infile, ((nqx + 1) / 2) * ((nqy + 1) / 2), scratch, b2);

    /*
     * insert in bitplane BIT of image A
     */
    //qtree_bitins64(scratch, nqx, nqy, a, n, bit);
}

/*  ############################################################################  */
/*
 * Huffman decoding for fixed codes
 *
 * Coded values range from 0-15
 *
 * Huffman code values (hex):
 *
 *	3e, 00, 01, 08, 02, 09, 1a, 1b,
 *	03, 1c, 0a, 1d, 0b, 1e, 3f, 0c
 *
 * and number of bits in each code:
 *
 *	6,  3,  3,  4,  3,  4,  5,  5,
 *	3,  5,  4,  5,  4,  5,  6,  4
 */
pub fn input_huffman(infile: &mut Cursor<&[u8]>, b2: &mut Buffer2) -> i32 {
    /*
     * get first 3 bits to start
     */
    let mut c: i32 = input_nbits(infile, 3, b2);
    if c < 4 {
        /*
         * this is all we need
         * return 1,2,4,8 for c=0,1,2,3
         */
        return 1 << c;
    }
    /*
     * get the next bit
     */
    c = input_bit(infile, b2) | (c << 1);
    if c < 13 {
        /*
         * OK, 4 bits is enough
         */
        match c {
            8 => return 3,
            9 => return 5,
            10 => return 10,
            11 => return 12,
            12 => return 15,
            _ => (),
        }
    }
    /*
     * get yet another bit
     */
    c = input_bit(infile, b2) | (c << 1);
    if c < 31 {
        /*
         * OK, 5 bits is enough
         */
        match c {
            26 => return 6,
            27 => return 7,
            28 => return 9,
            29 => return 11,
            30 => return 13,
            _ => (),
        }
    }
    /*
     * need the 6th bit
     */
    c = input_bit(infile, b2) | (c << 1);
    if c == 62 {
        0
    } else {
        14
    }
}

/*  ############################################################################  */
/*  ############################################################################  */
/* Copyright (c) 1993 Association of Universities for Research
 * in Astronomy. All rights reserved. Produced under National
 * Aeronautics and Space Administration Contract No. NAS5-26555.
 */
/* qread.c	Read binary data
 *
 * Programmer: R. White		Date: 11 March 1991
 */

#[must_use]
pub fn readint(infile: &mut Cursor<&[u8]>) -> i32 {
    let mut b: [u8; 4] = [0; 4];

    /* Read integer A one byte at a time from infile.
     *
     * This is portable from Vax to Sun since it eliminates the
     * need for byte-swapping.
     *
     *  This routine is only called to read the first 3 values
     *  in the compressed file, so it doesn't have to be
     *  super-efficient
     */
    /*
    for i in 0..4 {
        //for (i=0; i<4; i++) {
        qread(infile, &mut b[i..i], 1);
    }
    */
    qread(infile, &mut b, 4);

    let mut a: i32 = i32::from(b[0]);
    for i in 1..4 {
        //for (i=1; i<4; i++) {
        a = (a << 8) + i32::from(b[i]);
    }
    a
}

/*  ############################################################################  */
#[must_use]
pub fn readlonglong(infile: &mut Cursor<&[u8]>) -> i64 {
    let mut b: [u8; 8] = [0; 8];

    /* Read integer A one byte at a time from infile.
     *
     * This is portable from Vax to Sun since it eliminates the
     * need for byte-swapping.
     *
     *  This routine is only called to read the first 3 values
     *  in the compressed file, so it doesn't have to be
     *  super-efficient
     */
    for i in 0..8 {
        //for (i=0; i<8; i++) {
        qread(infile, &mut b[i..], 1);
    }
    let mut a: i64 = i64::from(b[0]);
    for i in 1..8 {
        //for (i=1; i<8; i++) {
        a = (a << 8) + i64::from(b[i]);
    }
    a
}

/*  ############################################################################  */
pub fn qread(file: &mut Cursor<&[u8]>, buffer: &mut [u8], n: usize) {
    /*
     * read n bytes from file into buffer
     *
     */
    file.copy_to_slice(&mut buffer[0..n]);
    //buffer[0..n].copy_from_slice(&file[0..n]);
    //file.advance(n);

    //nextchar += n;

    //memcpy(buffer, &file[nextchar], n);
    //nextchar += n;
}

/*  ############################################################################  */
/*  ############################################################################  */
/* Copyright (c) 1993 Association of Universities for Research
 * in Astronomy. All rights reserved. Produced under National
 * Aeronautics and Space Administration Contract No. NAS5-26555.
 */

/* BIT INPUT ROUTINES */

/* THE BIT BUFFER */
pub struct Buffer2 {
    pub buffer2: usize, /* Bits waiting to be input	*/
    pub bits_to_go: i32, /* Number of bits still in buffer */
                        //pub nextchar: usize,
}

/* INITIALIZE BIT INPUT */

/*  ############################################################################  */
#[must_use]
pub fn start_inputing_bits() -> Buffer2 {
    /*
     * Buffer starts out with no bits in it
     */
    Buffer2 {
        buffer2: 0, /* Buffer is empty to start	*/
        bits_to_go: 0, /* with				*/
                    // nextchar: 0,
    }
}

/*  ############################################################################  */
/* INPUT A BIT */

pub fn input_bit(infile: &mut Cursor<&[u8]>, b2: &mut Buffer2) -> i32 {
    if b2.bits_to_go == 0 {
        /* Read the next byte if no	*/

        b2.buffer2 = infile.get_u8() as usize;
        //b2.nextchar += 1;

        b2.bits_to_go = 8;
    }
    /*
     * Return the next bit
     */
    b2.bits_to_go -= 1;
    ((b2.buffer2 >> b2.bits_to_go) & 1) as i32
}

/*  ############################################################################  */
/* INPUT N BITS (N must be <= 8) */

pub fn input_nbits(infile: &mut Cursor<&[u8]>, n: usize, b2: &mut Buffer2) -> i32 {
    /* AND mask for retreiving the right-most n bits */
    let mask: [i32; 9] = [0, 1, 3, 7, 15, 31, 63, 127, 255];

    if b2.bits_to_go < n as i32 {
        /*
         * need another byte's worth of bits
         */

        b2.buffer2 = (b2.buffer2 << 8) | infile.get_u8() as usize;
        //b2.nextchar += 1;
        b2.bits_to_go += 8;
    }
    /*
     * now pick off the first n bits
     */
    b2.bits_to_go -= n as i32;

    /* there was a slight gain in speed by replacing the following line */
    /*	return( (buffer2>>bits_to_go) & ((1<<n)-1) ); */
    (b2.buffer2 >> b2.bits_to_go) as i32 & (mask[n])
}
/*  ############################################################################  */
/* INPUT 4 BITS  */

pub fn input_nybble(infile: &mut Cursor<&[u8]>, b2: &mut Buffer2) -> i32 {
    if b2.bits_to_go < 4 {
        /*
         * need another byte's worth of bits
         */

        b2.buffer2 = (b2.buffer2 << 8) | infile.get_u8() as usize;
        //b2.nextchar += 1;
        b2.bits_to_go += 8;
    }
    /*
     * now pick off the first 4 bits
     */
    b2.bits_to_go -= 4;

    ((b2.buffer2 >> b2.bits_to_go) & 15) as i32
}
/*  ############################################################################  */
/* INPUT array of 4 BITS  */

pub fn input_nnybble(
    infile: &mut Cursor<&[u8]>,
    n: usize,
    array: &mut [u8],
    b2: &mut Buffer2,
) -> i32 {
    /* copy n 4-bit nybbles from infile to the lower 4 bits of array */

    let mut _ii: usize;

    /*  forcing byte alignment doesn;t help, and even makes it go slightly slower
    if (bits_to_go != 8) input_nbits(infile, bits_to_go);
    */
    if n == 1 {
        array[0] = input_nybble(infile, b2) as u8;
        return 0;
    }

    if b2.bits_to_go == 8 {
        /*
           already have 2 full nybbles in buffer2, so
           backspace the infile array to reuse last char
        */
        // TODO
        infile.set_position(infile.position() - 1);
        //b2.nextchar -= 1;
        b2.bits_to_go = 0;
    }

    /* bits_to_go now has a value in the range 0 - 7.  After adding  */
    /* another byte, bits_to_go effectively will be in range 8 - 15 */

    let shift1: i32 = b2.bits_to_go + 4; /* shift1 will be in range 4 - 11 */
    let shift2: i32 = b2.bits_to_go; /* shift2 will be in range 0 -  7 */
    let mut kk: usize = 0;

    /* special case */
    if b2.bits_to_go == 0 {
        for _ii in 0..(n / 2) {
            // for (ii = 0; ii < n/2; ii++) {
            /*
             * refill the buffer with next byte
             */

            b2.buffer2 = (b2.buffer2 << 8) | infile.get_u8() as usize;
            // b2.nextchar += 1;
            array[kk] = ((b2.buffer2 >> 4) & 15) as u8;
            array[kk + 1] = ((b2.buffer2) & 15) as u8; /* no shift required */
            kk += 2;
        }
    } else {
        for _ii in 0..(n / 2) {
            //for (ii = 0; ii < n/2; ii++) {
            /*
             * refill the buffer with next byte
             */
            b2.buffer2 = (b2.buffer2 << 8) | infile.get_u8() as usize;
            //b2.nextchar += 1;
            array[kk] = ((b2.buffer2 >> shift1) & 15) as u8;
            array[kk + 1] = ((b2.buffer2 >> shift2) & 15) as u8;
            kk += 2;
        }
    }

    let ii = (n / 2) - 1;
    if ii * 2 != n {
        /* have to read last odd byte */
        array[n - 1] = input_nybble(infile, b2) as u8;
    }

    ((b2.buffer2 >> b2.bits_to_go) & 15) as i32
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_fits_decompress() {
        let input: [u8; 48] = [
            221, 153, 0, 0, 0, 4, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 32, 5, 5, 5, 245,
            231, 227, 199, 253, 227, 199, 253, 247, 255, 120, 249, 245, 239, 254, 241, 255, 124,
            120, 251, 0, 68, 200,
        ];
        let mut output: Vec<i32> = vec![0; 16];
        let res = fits_hdecompress(&input, 0, &mut output);

        assert_eq!(output.len(), 16);
        assert_eq!(res.1, 4);
        assert_eq!(res.2, 4);

        assert_eq!(output, [2, 2, 1, 2, 3, 2, 7, 7, 4, 2, 2, 1, 2, 4, 25, 2]);
        assert_eq!(
            input,
            [
                221, 153, 0, 0, 0, 4, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 32, 5, 5, 5,
                245, 231, 227, 199, 253, 227, 199, 253, 247, 255, 120, 249, 245, 239, 254, 241,
                255, 124, 120, 251, 0, 68, 200
            ]
        );
    }

    #[test]
    fn test_fits_decompress_strange() {
        let input: [u8; 101] = [
            221,153,0,0,0,10,0,0,0,1,0,0,0,0,255,255,255,255,255,255,109,64,16,17,0,2,136,255,191,224,40,143,251,254,246,207,253,238,168,251,53,238,168,255,223,247,253,255,127,223,247,253,255,127,223,247,253,255,127,223,247,253,255,127,223,247,253,255,127,223,247,178,255,239,251,217,127,247,178,253,151,255,127,222,234,143,186,163,255,123,170,62,234,143,186,163,238,168,255,192,100];
        let mut output: Vec<i32> = vec![0; 10];
        let res = fits_hdecompress(&input, 0, &mut output);

        assert_eq!(res.0, 0);
        assert_eq!(res.1, 10);
        assert_eq!(res.2, 1);
        assert_eq!(res.3, 0);

        assert_eq!(output.len(), 10);

        assert_eq!(output, [-1, -1, -112, -1, 9983, -28528, -112, -1, -1, -1]);
        assert_eq!(
            input,
            [
                221,153,0,0,0,10,0,0,0,1,0,0,0,0,255,255,255,255,255,255,109,64,16,17,0,2,136,255,191,224,40,143,251,254,246,207,253,238,168,251,53,238,168,255,223,247,253,255,127,223,247,253,255,127,223,247,253,255,127,223,247,253,255,127,223,247,253,255,127,223,247,178,255,239,251,217,127,247,178,253,151,255,127,222,234,143,186,163,255,123,170,62,234,143,186,163,238,168,255,192,100
            ]
        );
    }

    #[test]
    fn test_fits_decompress_strange2() {
        let input: [u8; 106] = [
            221,153,0,0,0,12,0,0,0,1,0,0,0,0,255,255,255,255,255,251,193,64,17,16,0,246,79,151,94,233,144,42,143,46,128,170,0,130,128,42,0,32,128,130,0,130,0,162,191,239,251,254,255,191,239,251,254,255,191,239,251,254,255,191,239,251,254,255,191,239,251,254,255,191,224,40,47,46,189,150,5,4,1,20,5,5,229,215,178,252,182,4,21,236,191,45,249,111,203,96,81,95,240,0,175,0];
        let mut output: Vec<i32> = vec![0; 12];
        let res = fits_hdecompress(&input, 0, &mut output);

        assert_eq!(res.0, 0);
        assert_eq!(res.1, 12);
        assert_eq!(res.2, 1);
        assert_eq!(res.3, 0);

        assert_eq!(output.len(), 12);

        assert_eq!(output, [-1,
            -1,
            -9584,
            -28561,
            -112,
            -24321,
            -1,
            -1,
            -1,
            -9584,
            -28561,
            -112]);
        assert_eq!(
            input,
            [
                221,153,0,0,0,12,0,0,0,1,0,0,0,0,255,255,255,255,255,251,193,64,17,16,0,246,79,151,94,233,144,42,143,46,128,170,0,130,128,42,0,32,128,130,0,130,0,162,191,239,251,254,255,191,239,251,254,255,191,239,251,254,255,191,239,251,254,255,191,239,251,254,255,191,224,40,47,46,189,150,5,4,1,20,5,5,229,215,178,252,182,4,21,236,191,45,249,111,203,96,81,95,240,0,175,0]
        
        );
    }


    #[test]
    fn test_qread() {
        let input: [u8; 4] = [1, 2, 3, 4];
        let mut input_c = Cursor::new(&input[..]);
        let mut buffer = [99, 99, 99, 99];

        qread(&mut input_c, &mut buffer, 2);

        assert_eq!(input_c.position(), 2);
        assert_eq!(buffer, [1, 2, 99, 99]);

        // Try in the middle now
        let mut input_c = Cursor::new(&input[..]);
        let mut buffer = [99, 99, 99, 99];

        qread(&mut input_c, &mut buffer[2..], 2);

        assert_eq!(input_c.position(), 2);
        assert_eq!(buffer, [99, 99, 1, 2]);
    }

    #[test]
    fn test_hinv() {
        let mut a = [32, 16, -2, 2, 12, 6, 0, -24, 2, 12, -1, -1, 0, 24, 4, -22];
        let nx = 4;
        let ny = 4;
        let smooth = 0;
        let scale = 0;

        let res = hinv(&mut a[0..], nx, ny, smooth, scale);

        assert_eq!(a, [2, 2, 1, 2, 3, 2, 7, 7, 4, 2, 2, 1, 2, 4, 25, 2]);
    }
}
