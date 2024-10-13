use std::io::Cursor;

use bytemuck::cast_slice_mut;
use bytes::Buf;

use crate::CODE_MAGIC;

fn ffpmsg(_m: &str) {}

#[derive(Debug, PartialEq)]
pub enum DecodeError {
    BadBitPlaneValues,
    MemoryAllocationError,
    BadFormatCode,
    BadFileFormat,
}

/// The bit buffer
struct Buffer2 {
    pub buffer2: i32,    // Bits waiting to be input
    pub bits_to_go: i32, // Number of bits still in buffer
}

pub struct HCDecoder {}

impl Default for HCDecoder {
    fn default() -> Self {
        Self::new()
    }
}

impl HCDecoder {
    pub fn new() -> Self {
        Self {}
    }

    /* ---------------------------------------------------------------------- */
    /// Decompress the input byte stream using the H-compress algorithm.
    ///
    /// # Arguments
    ///
    ///  `input`  - input array of compressed bytes
    ///  `a` - pre-allocated array to hold the output uncompressed image
    ///  `nx` - returned X axis size
    ///  `ny` - returned Y axis size
    ///
    /// NOTE: the nx and ny dimensions as defined within this code are reversed from
    /// the usual FITS notation.  ny is the fastest varying dimension, which is
    /// usually considered the X axis in the FITS image display
    pub fn read(
        &mut self,
        inputz: &[u8],
        smooth: i32,
        a: &mut [i32],
    ) -> Result<(usize, usize, i32), DecodeError> {
        let mut input = Cursor::new(inputz);

        // decode the input array
        let (nx, ny, scale) = decode(&mut input, a)?;

        // Un-Digitize
        undigitize(a, nx, ny, scale);

        // Inverse H-transform
        let stat = hinv(a, nx, ny, smooth, scale);

        if stat.is_err() {
            return Err(stat.unwrap_err());
        }

        Ok((nx, ny, scale))
    }

    /* ---------------------------------------------------------------------- */
    /// Decompress the input byte stream using the H-compress algorithm.
    ///
    /// # Arguments
    ///
    ///  `input`  - input array of compressed bytes
    ///  `a` - pre-allocated array to hold the output uncompressed image
    ///  `nx` - returned X axis size
    ///  `ny` - returned Y axis size
    ///
    /// NOTE: the nx and ny dimensions as defined within this code are reversed from
    /// the usual FITS notation.  ny is the fastest varying dimension, which is
    /// usually considered the X axis in the FITS image display
    pub fn read64(
        &mut self,
        inputz: &[u8],
        smooth: i32,
        a: &mut [i64],
    ) -> Result<(usize, usize, i32), DecodeError> {
        let mut input = Cursor::new(inputz);

        // decode the input array
        let (nx, ny, scale) = decode64(&mut input, a)?;

        // Un-Digitize
        undigitize64(a, nx, ny, scale);

        // Inverse H-transform
        let stat = hinv64(a, nx, ny, smooth, scale);

        if stat.is_err() {
            return Err(stat.unwrap_err());
        }

        // pack the I*8 values back into an I*4 array
        // Original code it does the conversion in place, but rust struggle with the borrow
        let nval = (nx) * (ny);
        let mut iarray = Vec::with_capacity(nval);

        for ii in 0..nval {
            iarray.push(a[ii] as i32);
        }

        // Put back into A
        let a_int: &mut [i32] = cast_slice_mut(a);
        a_int[..nval].copy_from_slice(&iarray[..nval]);

        Ok((nx, ny, scale))
    }
}

/*  ############################################################################  */
/// Inverse H-transform of NX x NY integer image
///
/// # Arguments
///
/// a - array of H-transform coefficients
/// nx - size of coefficient block to use
/// ny - actual 1st dimension of array
/// smooth - 0 for no smoothing, else smooth during inversion
/// scale - used if smoothing is specified
fn hinv(a: &mut [i32], nx: usize, ny: usize, smooth: i32, scale: i32) -> Result<(), DecodeError> {
    let mut lowbit0;
    let mut lowbit1;

    let mut h0: i32;
    let mut hx: i32;
    let mut hy: i32;
    let mut hc: i32;

    let mut oddx: usize;
    let mut oddy: usize;

    let mut s10: usize;
    let mut s00: usize;

    // log2n is log2 of i32::max(nx,ny) rounded up to next power of 2
    let nmax: usize = if nx > ny { nx } else { ny };
    let mut log2n: usize = ((nmax as f64).ln() / 2.0_f64.ln() + 0.5) as usize;

    if nmax > (1 << log2n) {
        log2n += 1;
    }

    // get temporary storage for shuffling elements
    let mut tmp: Vec<i32> = Vec::new();
    if tmp.try_reserve_exact((nmax + 1) / 2).is_err() {
        ffpmsg("hinv: insufficient memory");
        return Err(DecodeError::MemoryAllocationError);
    } else {
        tmp.resize((nmax + 1) / 2, 0);
    }

    // set up masks, rounding parameters
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

    // round h0 to multiple of bit2
    a[0] = (a[0] + (if a[0] >= 0 { prnd2 } else { nrnd2 })) & mask2;

    // do log2n expansions
    // We're indexing a as a 2-D array with dimensions (nx,ny).
    let mut nxtop = 1;
    let mut nytop = 1;
    let mut nxf = nx;
    let mut nyf = ny;
    let mut c = 1 << log2n;
    for k in (0..log2n).rev() {
        // this somewhat cryptic code generates the sequence
        // ntop[k-1] = (ntop[k]+1)/2, where ntop[log2n] = n
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

        // double shift and fix nrnd0 (because prnd0=0) on last pass
        if k == 0 {
            nrnd0 = 0;
            shift = 2;
        }

        // unshuffle in each dimension to interleave coefficients
        for i in 0..nxtop {
            unshuffle(&mut a[ny * i..], nytop, 1, &mut tmp);
        }

        for j in 0..nytop {
            unshuffle(&mut a[j..], nxtop, ny, &mut tmp);
        }

        // smooth by interpolating coefficients if SMOOTH != 0
        if smooth > 0 {
            hsmooth(a, nxtop, nytop, ny, scale);
        }

        oddx = nxtop % 2;
        oddy = nytop % 2;

        for i in (0..(nxtop - oddx)).step_by(2) {
            s00 = ny * i; // s00 is index of a[i,j]
            s10 = s00 + ny; // s10 is index of a[i+1,j]

            for _j in (0..(nytop - oddy)).step_by(2) {
                h0 = a[s00];
                hx = a[s10];
                hy = a[s00 + 1];
                hc = a[s10 + 1];

                // round hx and hy to multiple of bit1, hc to multiple of bit0
                // h0 is already a multiple of bit2
                hx = (hx + (if hx >= 0 { prnd1 } else { nrnd1 })) & mask1;
                hy = (hy + (if hy >= 0 { prnd1 } else { nrnd1 })) & mask1;
                hc = (hc + (if hc >= 0 { prnd0 } else { nrnd0 })) & mask0;

                // propagate bit0 of hc to hx,hy
                lowbit0 = hc & bit0;
                hx = if hx >= 0 { hx - lowbit0 } else { hx + lowbit0 };
                hy = if hy >= 0 { hy - lowbit0 } else { hy + lowbit0 };

                // Propagate bits 0 and 1 of hc,hx,hy to h0.
                // This could be simplified if we assume h0>0, but then
                // the inversion would not be lossless for images with
                // negative pixels.
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

                // Divide sums by 2 (4 last time)
                a[s10 + 1] = (h0 + hx + hy + hc) >> shift;
                a[s10] = (h0 + hx - hy - hc) >> shift;
                a[s00 + 1] = (h0 - hx + hy - hc) >> shift;
                a[s00] = (h0 - hx - hy + hc) >> shift;
                s00 += 2;
                s10 += 2;
            }

            if oddy > 0 {
                // do last element in row if row length is odd
                // s00+1, s10+1 are off edge
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
            // do last row if column length is odd
            // s10, s10+1 are off edge
            s00 = ny * (nxtop - 1);

            for _j in (0..(nytop - oddy)).step_by(2) {
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
                // do corner element if both row and column lengths are odd
                // s00+1, s10, s10+1 are off edge
                h0 = a[s00];
                a[s00] = h0 >> shift;
            }
        }

        // divide all the masks and rounding values by 2
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

    Ok(())
}

/*  ############################################################################  */
/// Inverse H-transform of NX x NY integer image
///
/// # Arguments
///
/// a - array of H-transform coefficients
/// nx - size of coefficient block to use
/// ny - actual 1st dimension of array
/// smooth - 0 for no smoothing, else smooth during inversion
/// scale - used if smoothing is specified
fn hinv64(a: &mut [i64], nx: usize, ny: usize, smooth: i32, scale: i32) -> Result<(), DecodeError> {
    let mut lowbit0;
    let mut lowbit1;

    let mut h0: i64;
    let mut hx: i64;
    let mut hy: i64;
    let mut hc: i64;

    let mut oddx: usize;
    let mut oddy: usize;

    let mut s10: usize;
    let mut s00: usize;

    // log2n is log2 of i32::max(nx,ny) rounded up to next power of 2
    let nmax: usize = if nx > ny { nx } else { ny };
    let mut log2n: usize = ((nmax as f64).ln() / 2.0_f64.ln() + 0.5) as usize;

    if nmax > (1 << log2n) {
        log2n += 1;
    }

    // get temporary storage for shuffling elements
    let mut tmp: Vec<i64> = Vec::new();
    if tmp.try_reserve_exact((nmax + 1) / 2).is_err() {
        ffpmsg("hinv64: insufficient memory");
        return Err(DecodeError::MemoryAllocationError);
    } else {
        tmp.resize((nmax + 1) / 2, 0);
    }

    // set up masks, rounding parameters
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

    // round h0 to multiple of bit2
    a[0] = (a[0] + (if a[0] >= 0 { prnd2 } else { nrnd2 })) & mask2;

    // do log2n expansions
    // We're indexing a as a 2-D array with dimensions (nx,ny).
    let mut nxtop = 1;
    let mut nytop = 1;
    let mut nxf = nx;
    let mut nyf = ny;
    let mut c = 1 << log2n;
    for k in (0..log2n).rev() {
        // this somewhat cryptic code generates the sequence
        // ntop[k-1] = (ntop[k]+1)/2, where ntop[log2n] = n
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

        // double shift and fix nrnd0 (because prnd0=0) on last pass
        if k == 0 {
            nrnd0 = 0;
            shift = 2;
        }

        // unshuffle in each dimension to interleave coefficients
        for i in 0..nxtop {
            unshuffle64(&mut a[ny * i..], nytop, 1, &mut tmp);
        }

        for j in 0..nytop {
            unshuffle64(&mut a[j..], nxtop, ny, &mut tmp);
        }

        // smooth by interpolating coefficients if SMOOTH != 0
        if smooth > 0 {
            hsmooth64(a, nxtop, nytop, ny, scale);
        }

        oddx = nxtop % 2;
        oddy = nytop % 2;

        for i in (0..(nxtop - oddx)).step_by(2) {
            s00 = ny * i; // s00 is index of a[i,j]
            s10 = s00 + ny; // s10 is index of a[i+1,j]

            for _j in (0..(nytop - oddy)).step_by(2) {
                h0 = a[s00];
                hx = a[s10];
                hy = a[s00 + 1];
                hc = a[s10 + 1];

                // round hx and hy to multiple of bit1, hc to multiple of bit0
                // h0 is already a multiple of bit2
                hx = (hx + (if hx >= 0 { prnd1 } else { nrnd1 })) & mask1;
                hy = (hy + (if hy >= 0 { prnd1 } else { nrnd1 })) & mask1;
                hc = (hc + (if hc >= 0 { prnd0 } else { nrnd0 })) & mask0;

                // propagate bit0 of hc to hx,hy
                lowbit0 = hc & bit0;
                hx = if hx >= 0 { hx - lowbit0 } else { hx + lowbit0 };
                hy = if hy >= 0 { hy - lowbit0 } else { hy + lowbit0 };

                // Propagate bits 0 and 1 of hc,hx,hy to h0.
                // This could be simplified if we assume h0>0, but then
                // the inversion would not be lossless for images with
                // negative pixels.
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

                // Divide sums by 2 (4 last time)
                a[s10 + 1] = (h0 + hx + hy + hc) >> shift;
                a[s10] = (h0 + hx - hy - hc) >> shift;
                a[s00 + 1] = (h0 - hx + hy - hc) >> shift;
                a[s00] = (h0 - hx - hy + hc) >> shift;
                s00 += 2;
                s10 += 2;
            }

            if oddy > 0 {
                // do last element in row if row length is odd
                // s00+1, s10+1 are off edge
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
            // do last row if column length is odd
            // s10, s10+1 are off edge
            s00 = ny * (nxtop - 1);

            for _j in (0..(nytop - oddy)).step_by(2) {
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
                // do corner element if both row and column lengths are odd
                // s00+1, s10, s10+1 are off edge
                h0 = a[s00];
                a[s00] = h0 >> shift;
            }
        }

        // divide all the masks and rounding values by 2
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

    Ok(())
}

/*  ############################################################################  */
/// Unshuffle an array of integers
///
/// # Arguments
///
/// `a` - array to shuffle
/// `n` - number of elements to shuffle
/// `n2` - second dimension
/// `tmp` - scratch storage
fn unshuffle(a: &mut [i32], n: usize, n2: usize, tmp: &mut [i32]) {
    // copy 2nd half of array to tmp
    let nhalf = (n + 1) >> 1;
    let mut pt = 0;
    let mut p1 = n2 * nhalf; // pointer to a[i]

    for _i in nhalf..n {
        tmp[pt] = a[p1];
        p1 += n2;
        pt += 1;
    }

    // distribute 1st half of array to even elements
    let mut p2 = n2 * (nhalf - 1); // pointer to a[i]
    p1 = (n2 * (nhalf - 1)) << 1; // pointer to a[2*i]

    for i in (0..nhalf).rev() {
        a[p1] = a[p2];
        if i == 0 {
            break;
        }
        p2 -= n2;
        p1 -= n2 + n2;
    }

    // now distribute 2nd half of array (in tmp) to odd elements
    pt = 0;
    p1 = n2; // pointer to a[i]
    for _i in (1..n).step_by(2) {
        a[p1] = tmp[pt];
        p1 += n2 + n2;
        pt += 1;
    }
}

/*  ############################################################################  */
/// Unshuffle an array of longs
///
/// # Arguments
///
/// `a` - array to shuffle
/// `n` - number of elements to shuffle
/// `n2` - second dimension
/// `tmp` - scratch storage
fn unshuffle64(a: &mut [i64], n: usize, n2: usize, tmp: &mut [i64]) {
    // copy 2nd half of array to tmp
    let nhalf = (n + 1) >> 1;
    let mut pt = 0;
    let mut p1 = n2 * nhalf; // pointer to a[i]

    for _i in nhalf..n {
        tmp[pt] = a[p1];
        p1 += n2;
        pt += 1;
    }

    // distribute 1st half of array to even elements
    let mut p2 = n2 * (nhalf - 1); // pointer to a[i]
    p1 = (n2 * (nhalf - 1)) << 1; // pointer to a[2*i]
    for _i in (0..nhalf).rev() {
        a[p1] = a[p2];
        if _i > 0 {
            p2 -= n2;
            p1 -= n2 + n2;
        }
    }

    // now distribute 2nd half of array (in tmp) to odd elements
    pt = 0;
    p1 = n2; // pointer to a[i]
    for _i in (1..n).step_by(2) {
        a[p1] = tmp[pt];
        p1 += n2 + n2;
        pt += 1;
    }
}

/*  ############################################################################  */
/// Smooth H-transform image by adjusting coefficients toward interpolated values
///
fn hsmooth(a: &mut [i32], nxtop: usize, nytop: usize, ny: usize, scale: i32)
/*
int a[];			 array of H-transform coefficients
int nxtop,nytop;	 size of coefficient block to use
int ny;				 actual 1st dimension of array
int scale;			 truncation scale factor that was used
*/
{
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

    // Maximum change in coefficients is determined by scale factor.
    // Since we rounded during division (see digitize.c), the biggest
    // permitted change is scale/2.
    let smax = scale >> 1;
    if smax <= 0 {
        return;
    }

    let ny2 = ny << 1;

    // We're indexing a as a 2-D array with dimensions (nxtop,ny) of which
    // only (nxtop,nytop) are used.  The coefficients on the edge of the
    // array are not adjusted (which is why the loops below start at 2
    // instead of 0 and end at nxtop-2 instead of nxtop.)

    // Adjust x difference hx
    for i in (2..(nxtop - 2)).step_by(2) {
        s00 = ny * i; // s00 is index of a[i,j]
        s10 = s00 + ny; // s10 is index of a[i+1,j]

        for _j in (0..nytop).step_by(2) {
            // hp is h0 (mean value) in next x zone, hm is h0 in previous x zone
            hm = a[s00 - ny2];
            h0 = a[s00];
            hp = a[s00 + ny2];

            // diff = 8 * hx slope that would match h0 in neighboring zones
            diff = hp - hm;

            // monotonicity constraints on diff
            dmax = i32::max(i32::min(hp - h0, h0 - hm), 0) << 2;
            dmin = i32::min(i32::max(hp - h0, h0 - hm), 0) << 2;

            // if monotonicity would set slope = 0 then don't change hx.
            // note dmax>=0, dmin<=0.
            if dmin < dmax {
                diff = i32::max(i32::min(diff, dmax), dmin);

                // Compute change in slope limited to range +/- smax.
                // Careful with rounding negative numbers when using
                // shift for divide by 8.

                s = diff - (a[s10] << 3);
                s = if s >= 0 { s >> 3 } else { (s + 7) >> 3 };
                s = i32::max(i32::min(s, smax), -smax);
                a[s10] += s;
            }
            s00 += 2;
            s10 += 2;
        }
    }

    // Adjust y difference hy
    for i in (0..nxtop).step_by(2) {
        s00 = ny * i + 2;
        s10 = s00 + ny;

        for _j in (2..(nytop - 2)).step_by(2) {
            hm = a[s00 - 2];
            h0 = a[s00];
            hp = a[s00 + 2];
            diff = hp - hm;
            dmax = i32::max(i32::min(hp - h0, h0 - hm), 0) << 2;
            dmin = i32::min(i32::max(hp - h0, h0 - hm), 0) << 2;
            if dmin < dmax {
                diff = i32::max(i32::min(diff, dmax), dmin);
                s = diff - (a[s00 + 1] << 3);
                s = if s >= 0 { s >> 3 } else { (s + 7) >> 3 };
                s = i32::max(i32::min(s, smax), -smax);
                a[s00 + 1] += s;
            }
            s00 += 2;
            s10 += 2;
        }
    }

    // Adjust curvature difference hc
    for i in (2..(nxtop - 2)).step_by(2) {
        s00 = ny * i + 2;
        s10 = s00 + ny;

        for _j in (2..(nytop - 2)).step_by(2) {
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

            // diff = 64 * hc value that would match h0 in neighboring zones
            diff = hpp + hmm - hmp - hpm;

            // 2 times x,y slopes in this zone
            hx2 = a[s10] << 1;
            hy2 = a[s00 + 1] << 1;

            // monotonicity constraints on diff
            m1 = i32::min(
                i32::max(hpp - h0, 0) - hx2 - hy2,
                i32::max(h0 - hpm, 0) + hx2 - hy2,
            );
            m2 = i32::min(
                i32::max(h0 - hmp, 0) - hx2 + hy2,
                i32::max(hmm - h0, 0) + hx2 + hy2,
            );
            dmax = i32::min(m1, m2) << 4;
            m1 = i32::max(
                i32::min(hpp - h0, 0) - hx2 - hy2,
                i32::min(h0 - hpm, 0) + hx2 - hy2,
            );
            m2 = i32::max(
                i32::min(h0 - hmp, 0) - hx2 + hy2,
                i32::min(hmm - h0, 0) + hx2 + hy2,
            );
            dmin = i32::max(m1, m2) << 4;

            // if monotonicity would set slope = 0 then don't change hc.
            // note dmax>=0, dmin<=0.
            if dmin < dmax {
                diff = i32::max(i32::min(diff, dmax), dmin);

                // Compute change in slope limited to range +/- smax.
                // Careful with rounding negative numbers when using
                // shift for divide by 64.

                s = diff - (a[s10 + 1] << 6);
                s = if s >= 0 { s >> 6 } else { (s + 63) >> 6 };
                s = i32::max(i32::min(s, smax), -smax);
                a[s10 + 1] += s;
            }
            s00 += 2;
            s10 += 2;
        }
    }
}

/*  ############################################################################  */
/// Smooth H-transform image by adjusting coefficients toward interpolated values
///
fn hsmooth64(a: &mut [i64], nxtop: usize, nytop: usize, ny: usize, scale: i32)
/*
LONGLONG a[];			 array of H-transform coefficients
int nxtop,nytop;	 size of coefficient block to use
int ny;				 actual 1st dimension of array
int scale;			 truncation scale factor that was used
*/
{
    // Maximum change in coefficients is determined by scale factor.
    // Since we rounded during division (see digitize.c), the biggest
    // permitted change is scale/2.
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

    // We're indexing a as a 2-D array with dimensions (nxtop,ny) of which
    // only (nxtop,nytop) are used.  The coefficients on the edge of the
    // array are not adjusted (which is why the loops below start at 2
    // instead of 0 and end at nxtop-2 instead of nxtop.)

    // Adjust x difference hx
    for i in (2..(nxtop - 2)).step_by(2) {
        s00 = ny * i; // s00 is index of a[i,j]
        s10 = s00 + ny; // s10 is index of a[i+1,j]

        for _j in (0..nytop).step_by(2) {
            // hp is h0 (mean value) in next x zone, hm is h0 in previous x zone
            hm = a[s00 - ny2];
            h0 = a[s00];
            hp = a[s00 + ny2];

            // diff = 8 * hx slope that would match h0 in neighboring zones
            diff = hp - hm;

            // monotonicity constraints on diff
            dmax = i64::max(i64::min(hp - h0, h0 - hm), 0) << 2;
            dmin = i64::min(i64::max(hp - h0, h0 - hm), 0) << 2;

            // if monotonicity would set slope = 0 then don't change hx.
            // note dmax>=0, dmin<=0.
            if dmin < dmax {
                diff = i64::max(i64::min(diff, dmax), dmin);

                // Compute change in slope limited to range +/- smax.
                // Careful with rounding negative numbers when using
                // shift for divide by 8.

                s = diff - (a[s10] << 3);
                s = if s >= 0 { s >> 3 } else { (s + 7) >> 3 };
                s = i64::max(i64::min(s, smax), -smax);
                a[s10] += s;
            }
            s00 += 2;
            s10 += 2;
        }
    }

    // Adjust y difference hy
    for i in (0..nxtop).step_by(2) {
        s00 = ny * i + 2;
        s10 = s00 + ny;

        for _j in (2..(nytop - 2)).step_by(2) {
            hm = a[s00 - 2];
            h0 = a[s00];
            hp = a[s00 + 2];
            diff = hp - hm;
            dmax = i64::max(i64::min(hp - h0, h0 - hm), 0) << 2;
            dmin = i64::min(i64::max(hp - h0, h0 - hm), 0) << 2;
            if dmin < dmax {
                diff = i64::max(i64::min(diff, dmax), dmin);
                s = diff - (a[s00 + 1] << 3);
                s = if s >= 0 { s >> 3 } else { (s + 7) >> 3 };
                s = i64::max(i64::min(s, smax), -smax);
                a[s00 + 1] += s;
            }
            s00 += 2;
            s10 += 2;
        }
    }

    // Adjust curvature difference hc
    for i in (2..(nxtop - 2)).step_by(2) {
        s00 = ny * i + 2;
        s10 = s00 + ny;

        for _j in (2..(nytop - 2)).step_by(2) {
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

            // diff = 64 * hc value that would match h0 in neighboring zones
            diff = hpp + hmm - hmp - hpm;

            // 2 times x,y slopes in this zone
            hx2 = a[s10] << 1;
            hy2 = a[s00 + 1] << 1;

            // monotonicity constraints on diff
            m1 = i64::min(
                i64::max(hpp - h0, 0) - hx2 - hy2,
                i64::max(h0 - hpm, 0) + hx2 - hy2,
            );
            m2 = i64::min(
                i64::max(h0 - hmp, 0) - hx2 + hy2,
                i64::max(hmm - h0, 0) + hx2 + hy2,
            );
            dmax = i64::min(m1, m2) << 4;
            m1 = i64::max(
                i64::min(hpp - h0, 0) - hx2 - hy2,
                i64::min(h0 - hpm, 0) + hx2 - hy2,
            );
            m2 = i64::max(
                i64::min(h0 - hmp, 0) - hx2 + hy2,
                i64::min(hmm - h0, 0) + hx2 + hy2,
            );
            dmin = i64::max(m1, m2) << 4;

            // if monotonicity would set slope = 0 then don't change hc.
            // note dmax>=0, dmin<=0.
            if dmin < dmax {
                diff = i64::max(i64::min(diff, dmax), dmin);

                // Compute change in slope limited to range +/- smax.
                // Careful with rounding negative numbers when using
                // shift for divide by 64.

                s = diff - (a[s10 + 1] << 6);
                s = if s >= 0 { s >> 6 } else { (s + 63) >> 6 };
                s = i64::max(i64::min(s, smax), -smax);
                a[s10 + 1] += s;
            }
            s00 += 2;
            s10 += 2;
        }
    }
}

/*  ############################################################################  */
/// undigitize H-transform
fn undigitize(a: &mut [i32], nx: usize, ny: usize, scale: i32) {
    // multiply by scale
    if scale <= 1 {
        return;
    }

    if nx == 0 || ny == 0 {
        return;
    }

    for item in a.iter_mut().take(nx * ny - 1) {
        *item *= scale;
    }
}

/*  ############################################################################  */
/// undigitize H-transform
fn undigitize64(a: &mut [i64], nx: usize, ny: usize, scale: i32) {
    // multiply by scale
    if scale <= 1 {
        return;
    }

    let scale64: i64 = i64::from(scale); // use a 64-bit int for efficiency in the big loop

    if nx == 0 || ny == 0 {
        return;
    }

    a.iter_mut().take(nx * ny - 1).for_each(|x| *x *= scale64);
}

/*  ############################################################################  */
/// read codes from infile and construct array
fn decode(infile: &mut Cursor<&[u8]>, a: &mut [i32]) -> Result<(usize, usize, i32), DecodeError> {
    //nx: usize, ny: usize, scale: i32
    /*
    char *infile;				 input file
    int  *a;				 address of output array [nx][ny]
    int  *nx,*ny;				 size of output array
    int  *scale;				 scale factor for digitization
    */

    let mut nbitplanes: [u8; 3] = [0; 3];
    let mut tmagic: [u8; 2] = [0; 2];
    // initialize the byte read position to the beginning of the array

    // File starts either with special 2-byte magic code or with
    // FITS keyword "SIMPLE  ="
    qread(infile, &mut tmagic, 2);

    // check for correct magic code value
    if !tmagic.eq(&CODE_MAGIC) {
        ffpmsg("bad file format");
        return Err(DecodeError::BadFileFormat);
    }

    let nx = readint(infile) as usize; // x size of image
    let ny = readint(infile) as usize; // y size of image
    let scale = readint(infile); // scale factor for digitization

    // sum of all pixels
    let sumall = readlonglong(infile);

    // # bits in quadrants
    let len = nbitplanes.len();
    qread(infile, &mut nbitplanes, len);

    let stat = dodecode(infile, a, nx, ny, nbitplanes);

    // put sum of all pixels back into pixel 0
    a[0] = sumall as i32;

    match stat {
        Ok(_) => Ok((nx, ny, scale)),
        Err(e) => Err(e),
    }
}

/*  ############################################################################  */
/// read codes from infile and construct array
fn decode64(infile: &mut Cursor<&[u8]>, a: &mut [i64]) -> Result<(usize, usize, i32), DecodeError> {
    //nx: usize, ny: usize, scale: i32
    /*
    char *infile;				 input file
    i64  *a;				 address of output array [nx][ny]
    int  *nx,*ny;				 size of output array
    int  *scale;				 scale factor for digitization
    */

    let mut nbitplanes: [u8; 3] = [0; 3];
    let mut tmagic: [u8; 2] = [0; 2];
    // initialize the byte read position to the beginning of the array

    // File starts either with special 2-byte magic code or with
    // FITS keyword "SIMPLE  ="
    qread(infile, &mut tmagic, 2);

    // check for correct magic code value
    if !tmagic.eq(&CODE_MAGIC) {
        ffpmsg("bad file format");
        return Err(DecodeError::BadFileFormat);
    }

    let nx = readint(infile) as usize; // x size of image
    let ny = readint(infile) as usize; // y size of image
    let scale = readint(infile); // scale factor for digitization

    // sum of all pixels
    let sumall = readlonglong(infile);

    // # bits in quadrants
    let len = nbitplanes.len();
    qread(infile, &mut nbitplanes, len);

    let stat = dodecode64(infile, a, nx, ny, nbitplanes);

    // put sum of all pixels back into pixel 0
    a[0] = sumall;

    match stat {
        Ok(_) => Ok((nx, ny, scale)),
        Err(e) => Err(e),
    }
}

/*  ############################################################################  */
/// Decode stream of characters on infile and return array
///
/// # Arguments
///
/// `infile` - input file
/// `a` - output array
/// `nx` - x dimension
/// `ny` - y dimension
/// `nbitplanes` - number of bit planes in quadrants
fn dodecode(
    infile: &mut Cursor<&[u8]>,
    a: &mut [i32],
    nx: usize,
    ny: usize,
    nbitplanes: [u8; 3],
) -> Result<(), DecodeError> {
    let nel = nx * ny;
    let nx2 = (nx + 1) / 2;
    let ny2 = (ny + 1) / 2;

    // initialize a to zero
    a.iter_mut().take(nel).for_each(|x| *x = 0);

    // Initialize bit input
    let mut b2 = start_inputing_bits();

    // read bit planes for each quadrant
    qtree_decode(infile, a, ny, nx2, ny2, i32::from(nbitplanes[0]), &mut b2)?;

    qtree_decode(
        infile,
        &mut a[ny2..],
        ny,
        nx2,
        ny / 2,
        i32::from(nbitplanes[1]),
        &mut b2,
    )?;

    qtree_decode(
        infile,
        &mut a[(ny * nx2)..],
        ny,
        nx / 2,
        ny2,
        i32::from(nbitplanes[1]),
        &mut b2,
    )?;

    if ny * nx2 + ny2 < a.len() {
        qtree_decode(
            infile,
            &mut a[(ny * nx2 + ny2)..],
            ny,
            nx / 2,
            ny / 2,
            i32::from(nbitplanes[2]),
            &mut b2,
        )?;
    }

    // make sure there is an EOF symbol (nybble=0) at end
    if input_nybble(infile, &mut b2) != 0 {
        ffpmsg("dodecode: bad bit plane values");
        return Err(DecodeError::BadBitPlaneValues);
    }

    // now get the sign bits
    // Re-initialize bit input
    let mut b2 = start_inputing_bits();
    for item in a.iter_mut().take(nel) {
        if *item > 0 {
            // tried putting the input_bit code in-line here, instead of
            // calling the function, but it made no difference in the speed
            if input_bit(infile, &mut b2) > 0 {
                *item = -*item;
            }
        }
    }
    Ok(())
}

/*  ############################################################################  */
/// Decode stream of characters on infile and return array
///
/// # Arguments
///
/// `infile` - input file
/// `a` - output array
/// `nx` - x dimension
/// `ny` - y dimension
/// `nbitplanes` - number of bit planes in quadrants
fn dodecode64(
    infile: &mut Cursor<&[u8]>,
    a: &mut [i64],
    nx: usize,
    ny: usize,
    nbitplanes: [u8; 3],
) -> Result<(), DecodeError> {
    let nel = nx * ny;
    let nx2 = (nx + 1) / 2;
    let ny2 = (ny + 1) / 2;

    // initialize a to zero
    a.iter_mut().take(nel).for_each(|x| *x = 0);

    // Initialize bit input
    let mut b2 = start_inputing_bits();

    // read bit planes for each quadrant
    qtree_decode64(infile, a, ny, nx2, ny2, i32::from(nbitplanes[0]), &mut b2)?;

    qtree_decode64(
        infile,
        &mut a[ny2..],
        ny,
        nx2,
        ny / 2,
        i32::from(nbitplanes[1]),
        &mut b2,
    )?;

    qtree_decode64(
        infile,
        &mut a[(ny * nx2)..],
        ny,
        nx / 2,
        ny2,
        i32::from(nbitplanes[1]),
        &mut b2,
    )?;

    if ny * nx2 + ny2 < a.len() {
        qtree_decode64(
            infile,
            &mut a[(ny * nx2 + ny2)..],
            ny,
            nx / 2,
            ny / 2,
            i32::from(nbitplanes[2]),
            &mut b2,
        )?;
    }

    // make sure there is an EOF symbol (nybble=0) at end
    if input_nybble(infile, &mut b2) != 0 {
        ffpmsg("dodecode64: bad bit plane values");
        return Err(DecodeError::BadBitPlaneValues);
    }

    // now get the sign bits
    // Re-initialize bit input
    let mut b2 = start_inputing_bits();
    for item in a.iter_mut().take(nel) {
        if *item > 0 && input_bit(infile, &mut b2) > 0 {
            *item = -*item;
        }
    }
    Ok(())
}

/*  ############################################################################  */
/// Read stream of codes from infile and construct bit planes in quadrant of 2-D array using binary quadtree coding
fn qtree_decode(
    infile: &mut Cursor<&[u8]>,
    a: &mut [i32],
    n: usize,
    nqx: usize,
    nqy: usize,
    nbitplanes: i32,
    b2: &mut Buffer2,
) -> Result<(), DecodeError> {
    /*
    char *infile;
    int a[];				 a is 2-D array with dimensions (n,n)
    int n;					 length of full row in a
    int nqx;				 partial length of row to decode
    int nqy;				 partial length of column (<=n)
    int nbitplanes;				 number of bitplanes to decode
    */

    let mut b;
    let mut nfx;
    let mut nfy;
    let mut c;

    let mut nx: usize;
    let mut ny: usize;

    // log2n is log2 of i32::max(nqx,nqy) rounded up to next power of 2
    let nqmax: usize = if nqx > nqy { nqx } else { nqy };
    let mut log2n: usize = ((nqmax as f32).ln() / 2.0_f32.ln() + 0.5) as usize;
    if nqmax > (1 << log2n) {
        log2n += 1;
    }

    // allocate scratch array for working space
    let nqx2 = (nqx + 1) / 2;
    let nqy2 = (nqy + 1) / 2;

    let mut scratch: Vec<u8> = Vec::new();
    let scratch_len = (nqx2 + 1) * (nqy2 + 1);
    if scratch.try_reserve_exact(scratch_len).is_err() {
        ffpmsg("qtree_decode: memory allocation error");
        return Err(DecodeError::MemoryAllocationError);
    } else {
        scratch.resize(scratch_len, 0);
    }

    //let mut scratch2: Vec<u8> = vec![0; nqx2 * nqy2+20];

    // now decode each bit plane, starting at the top
    // A is assumed to be initialized to zero
    for bit in (0..nbitplanes).rev() {
        // Was bitplane was quadtree-coded or written directly?
        b = input_nybble(infile, b2);

        if b == 0 {
            // bit map was written directly
            read_bdirect(infile, a, n, nqx, nqy, &mut scratch, bit, b2);
        } else if b != 0xf {
            ffpmsg("qtree_decode: bad format code");
            return Err(DecodeError::BadFormatCode);
        } else {
            // bitmap was quadtree-coded, do log2n expansions

            // read first code
            scratch[0] = input_huffman(infile, b2) as u8;

            // now do log2n expansions, reading codes from file as necessary
            nx = 1;
            ny = 1;
            nfx = nqx;
            nfy = nqy;
            c = 1 << log2n;
            for _k in 1..log2n {
                // this somewhat cryptic code generates the sequence
                // n[k-1] = (n[k]+1)/2 where n[log2n]=nqx or nqy
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
                qtree_expand(infile, &mut scratch, nx, ny, b2);
                // scratch.copy_from_slice(&scratch2);
            }

            // now copy last set of 4-bit codes to bitplane bit of array a
            qtree_bitins(&mut scratch, nqx, nqy, a, n, bit);
        }
    }
    Ok(())
}

/*  ############################################################################  */
/// Read stream of codes from infile and construct bit planes in quadrant of 2-D array using binary quadtree coding
fn qtree_decode64(
    infile: &mut Cursor<&[u8]>,
    a: &mut [i64],
    n: usize,
    nqx: usize,
    nqy: usize,
    nbitplanes: i32,
    b2: &mut Buffer2,
) -> Result<(), DecodeError> {
    /*
    char *infile;
    int a[];				 a is 2-D array with dimensions (n,n)
    int n;					 length of full row in a
    int nqx;				 partial length of row to decode
    int nqy;				 partial length of column (<=n)
    int nbitplanes;				 number of bitplanes to decode
    */

    let mut b;
    let mut nfx;
    let mut nfy;
    let mut c;

    let mut nx: usize;
    let mut ny: usize;

    // log2n is log2 of i32::max(nqx,nqy) rounded up to next power of 2
    let nqmax: usize = if nqx > nqy { nqx } else { nqy };
    let mut log2n: usize = ((nqmax as f32).ln() / 2.0_f32.ln() + 0.5) as usize;
    if nqmax > (1 << log2n) {
        log2n += 1;
    }

    // allocate scratch array for working space
    let nqx2 = (nqx + 1) / 2;
    let nqy2 = (nqy + 1) / 2;

    let mut scratch: Vec<u8> = Vec::new();
    let scratch_len = (nqx2 + 1) * (nqy2 + 1);
    if scratch.try_reserve_exact(scratch_len).is_err() {
        ffpmsg("qtree_decode64: memory allocation error");
        return Err(DecodeError::MemoryAllocationError);
    } else {
        scratch.resize(scratch_len, 0);
    }

    //let mut scratch2: Vec<u8> = vec![0; nqx2 * nqy2+20];

    // now decode each bit plane, starting at the top
    // A is assumed to be initialized to zero
    for bit in (0..nbitplanes).rev() {
        // Was bitplane was quadtree-coded or written directly?
        b = input_nybble(infile, b2);

        if b == 0 {
            // bit map was written directly
            read_bdirect64(infile, a, n, nqx, nqy, &mut scratch, bit, b2);
        } else if b != 0xf {
            ffpmsg("qtree_decode64: bad format code");
            return Err(DecodeError::BadFormatCode);
        } else {
            // bitmap was quadtree-coded, do log2n expansions

            // read first code
            scratch[0] = input_huffman(infile, b2) as u8;

            // now do log2n expansions, reading codes from file as necessary
            nx = 1;
            ny = 1;
            nfx = nqx;
            nfy = nqy;
            c = 1 << log2n;
            for _k in 1..log2n {
                // this somewhat cryptic code generates the sequence
                // n[k-1] = (n[k]+1)/2 where n[log2n]=nqx or nqy
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
                qtree_expand(infile, &mut scratch, nx, ny, b2);
                // scratch.copy_from_slice(&scratch2);
            }

            // now copy last set of 4-bit codes to bitplane bit of array a
            qtree_bitins64(&mut scratch, nqx, nqy, a, n, bit);
        }
    }
    Ok(())
}

/*  ############################################################################  */
/// Do one quadtree expansion step on array a[(nqx+1)/2,(nqy+1)/2]
///
/// Results put into b[nqx,nqy] (which may be the same as a)
fn qtree_expand(infile: &mut Cursor<&[u8]>, a: &mut [u8], nx: usize, ny: usize, b2: &mut Buffer2) {
    // first copy a to b, expanding each 4-bit value
    qtree_copy(a, nx, ny, ny);

    // now read new 4-bit values into b for each non-zero element
    for i in (0..(nx * ny)).rev() {
        if a[i] > 0 {
            a[i] = input_huffman(infile, b2) as u8;
        }
    }
}

/*  ############################################################################  */
/// copy 4-bit values from a[(nx+1)/2,(ny+1)/2] to b[nx,ny], expanding each value to 2x2 pixels
///
/// a,b may be same array
///
///
fn qtree_copy(a: &mut [u8], nx: usize, ny: usize, n: usize)
/*   int n;		declared y dimension of b */
{
    let mut s00: usize;
    let mut s10: usize;

    // first copy 4-bit values to b
    // start at end in case a,b are same array
    let nx2 = (nx + 1) / 2;
    let ny2 = (ny + 1) / 2;

    let mut k = ny2 * (nx2 - 1) + ny2; // k   is index of a[i,j]
    for i in (0..nx2).rev() {
        //TODO added +2 below to re-order the -2 to prevent underflow
        s00 = 2 * (n * i + ny2 - 1) + 2; // s00 is index of b[2*i,2*j]
        for _j in (0..ny2).rev() {
            k -= 1;
            s00 -= 2;
            a[s00] = a[k];
        }
    }

    // now expand each 2x2 block
    let oddx = nx % 2;
    let oddy = ny % 2;

    if nx == 0 {
        return;
    } // Return early and prevent underflow

    for i in (0..(nx - 1)).step_by(2) {
        // Note:
        // Unlike the case in qtree_bitins, this code runs faster on a 32-bit linux
        // machine using the s10 intermediate variable, rather that using s00+n.
        // Go figure!

        s00 = n * i; // s00 is index of b[i,j]
        s10 = s00 + n; // s10 is index of b[i+1,j]

        if ny == 0 {
            continue;
        } // Continue early and prevent underflow

        for _j in (0..(ny - 1)).step_by(2) {
            match a[s00] {
                0 => {
                    a[s10 + 1] = 0;
                    a[s10] = 0;
                    a[s00 + 1] = 0;
                    a[s00] = 0;
                }

                1 => {
                    a[s10 + 1] = 1;
                    a[s10] = 0;
                    a[s00 + 1] = 0;
                    a[s00] = 0;
                }
                2 => {
                    a[s10 + 1] = 0;
                    a[s10] = 1;
                    a[s00 + 1] = 0;
                    a[s00] = 0;
                }
                3 => {
                    a[s10 + 1] = 1;
                    a[s10] = 1;
                    a[s00 + 1] = 0;
                    a[s00] = 0;
                }
                4 => {
                    a[s10 + 1] = 0;
                    a[s10] = 0;
                    a[s00 + 1] = 1;
                    a[s00] = 0;
                }
                5 => {
                    a[s10 + 1] = 1;
                    a[s10] = 0;
                    a[s00 + 1] = 1;
                    a[s00] = 0;
                }
                6 => {
                    a[s10 + 1] = 0;
                    a[s10] = 1;
                    a[s00 + 1] = 1;
                    a[s00] = 0;
                }
                7 => {
                    a[s10 + 1] = 1;
                    a[s10] = 1;
                    a[s00 + 1] = 1;
                    a[s00] = 0;
                }
                8 => {
                    a[s10 + 1] = 0;
                    a[s10] = 0;
                    a[s00 + 1] = 0;
                    a[s00] = 1;
                }
                9 => {
                    a[s10 + 1] = 1;
                    a[s10] = 0;
                    a[s00 + 1] = 0;
                    a[s00] = 1;
                }
                10 => {
                    a[s10 + 1] = 0;
                    a[s10] = 1;
                    a[s00 + 1] = 0;
                    a[s00] = 1;
                }
                11 => {
                    a[s10 + 1] = 1;
                    a[s10] = 1;
                    a[s00 + 1] = 0;
                    a[s00] = 1;
                }
                12 => {
                    a[s10 + 1] = 0;
                    a[s10] = 0;
                    a[s00 + 1] = 1;
                    a[s00] = 1;
                }
                13 => {
                    a[s10 + 1] = 1;
                    a[s10] = 0;
                    a[s00 + 1] = 1;
                    a[s00] = 1;
                }
                14 => {
                    a[s10 + 1] = 0;
                    a[s10] = 1;
                    a[s00 + 1] = 1;
                    a[s00] = 1;
                }
                15 => {
                    a[s10 + 1] = 1;
                    a[s10] = 1;
                    a[s00 + 1] = 1;
                    a[s00] = 1;
                }
                _ => (),
            }

            /*
                        a[s10+1] =  a[s00]     & 1;
                        a[s10  ] = (a[s00]>>1) & 1;
                        a[s00+1] = (a[s00]>>2) & 1;
                        a[s00  ] = (a[s00]>>3) & 1;
            */

            s00 += 2;
            s10 += 2;
        }

        if oddy > 0 {
            // row size is odd, do last element in row
            // s00+1, s10+1 are off edge
            // not worth converting this to use 16 case statements
            a[s10] = (a[s00] >> 1) & 1;
            a[s00] = (a[s00] >> 3) & 1;
        }
    }

    if oddx > 0 {
        // column size is odd, do last row
        // s10, s10+1 are off edge
        s00 = n * (nx - 1);

        if ny == 0 {
            return;
        } // Return early and prevent underflow

        for _j in (0..(ny - 1)).step_by(2) {
            // not worth converting this to use 16 case statements
            a[s00 + 1] = (a[s00] >> 2) & 1;
            a[s00] = (a[s00] >> 3) & 1;
            s00 += 2;
        }

        if oddy > 0 {
            // both row and column size are odd, do corner element
            // s00+1, s10, s10+1 are off edge
            // not worth converting this to use 16 case statements
            a[s00] = (a[s00] >> 3) & 1;
        }
    }
}

/*  ############################################################################  */
/// Copy 4-bit values from a[(nx+1)/2,(ny+1)/2] to b[nx,ny], expanding each value to 2x2 pixels and inserting into bitplane BIT of B.
///
/// A,B may NOT be same array (it wouldn't make sense to be inserting
/// bits into the same array anyway.)
///
/// # Arguments
///
/// * `a` - 1-D array of 4-bit values
/// * `nx` - x dimension of a
/// * `ny` - y dimension of a
/// * `b` - 1-D array of 32-bit values
/// * `n` - y dimension of b
/// * `bit` - bitplane to insert into
///
fn qtree_bitins(a: &mut [u8], nx: usize, ny: usize, b: &mut [i32], n: usize, bit: i32) {
    let mut s00: usize;

    let plane_val = 1 << bit;

    // expand each 2x2 block
    let mut k: usize = 0; // k   is index of a[i/2,j/2]
    let oddx = nx % 2;
    let oddy = ny % 2;

    if nx == 0 {
        return;
    } // Return early and prevent underflow

    for i in (0..(nx - 1)).step_by(2) {
        s00 = n * i; // s00 is index of b[i,j]

        // Note:
        // this code appears to run very slightly faster on a 32-bit linux
        // machine using s00+n rather than the s10 intermediate variable

        //		s10 = s00+n;
        // s10 is index of b[i+1,j]
        if ny == 0 {
            // TODO
            continue;
        } // Continue early and prevent underflow

        for _j in (0..(ny - 1)).step_by(2) {
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
            //			s10 += 2;
            k += 1;
        }

        if oddy > 0 {
            // row size is odd, do last element in row
            // s00+1, s10+1 are off edge

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
        // column size is odd, do last row
        // s10, s10+1 are off edge
        s00 = n * (nx - 1);

        if ny == 0 {
            return;
        } // Return early and prevent underflow

        for _j in (0..(ny - 1)).step_by(2) {
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
            // both row and column size are odd, do corner element
            // s00+1, s10, s10+1 are off edge

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
/// Copy 4-bit values from a[(nx+1)/2,(ny+1)/2] to b[nx,ny], expanding each value to 2x2 pixels and inserting into bitplane BIT of B.
///
/// A,B may NOT be same array (it wouldn't make sense to be inserting
/// bits into the same array anyway.)
///
/// # Arguments
///
/// * `a` - 1-D array of 4-bit values
/// * `nx` - x dimension of a
/// * `ny` - y dimension of a
/// * `b` - 1-D array of 32-bit values
/// * `n` - y dimension of b
/// * `bit` - bitplane to insert into
///
fn qtree_bitins64(a: &mut [u8], nx: usize, ny: usize, b: &mut [i64], n: usize, bit: i32) {
    let mut s00: usize;

    let plane_val: i64 = 1 << bit;

    // expand each 2x2 block
    let mut k: usize = 0; // k   is index of a[i/2,j/2]
    let oddx = nx % 2;
    let oddy = ny % 2;

    if nx == 0 {
        return;
    } // Return early and prevent underflow

    for i in (0..(nx - 1)).step_by(2) {
        s00 = n * i; // s00 is index of b[i,j]

        // Note:
        // this code appears to run very slightly faster on a 32-bit linux
        // machine using s00+n rather than the s10 intermediate variable

        //		s10 = s00+n;
        // s10 is index of b[i+1,j]
        if ny == 0 {
            // TODO
            continue;
        } // Continue early and prevent underflow

        for _j in (0..(ny - 1)).step_by(2) {
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
            //			s10 += 2;
            k += 1;
        }

        if oddy > 0 {
            // row size is odd, do last element in row
            // s00+1, s10+1 are off edge

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
        // column size is odd, do last row
        // s10, s10+1 are off edge
        s00 = n * (nx - 1);

        if ny == 0 {
            return;
        } // Return early and prevent underflow

        for _j in (0..(ny - 1)).step_by(2) {
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
            // both row and column size are odd, do corner element
            // s00+1, s10, s10+1 are off edge

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
fn read_bdirect(
    infile: &mut Cursor<&[u8]>,
    a: &mut [i32],
    n: usize,
    nqx: usize,
    nqy: usize,
    scratch: &mut [u8],
    bit: i32,
    b2: &mut Buffer2,
) {
    // read bit image packed 4 pixels/nybble
    input_nnybble(infile, ((nqx + 1) / 2) * ((nqy + 1) / 2), scratch, b2);

    // insert in bitplane BIT of image A
    qtree_bitins(scratch, nqx, nqy, a, n, bit);
}

/*  ############################################################################  */
fn read_bdirect64(
    infile: &mut Cursor<&[u8]>,
    a: &mut [i64],
    n: usize,
    nqx: usize,
    nqy: usize,
    scratch: &mut [u8],
    bit: i32,
    b2: &mut Buffer2,
) {
    // read bit image packed 4 pixels/nybble
    input_nnybble(infile, ((nqx + 1) / 2) * ((nqy + 1) / 2), scratch, b2);

    // insert in bitplane BIT of image A
    qtree_bitins64(scratch, nqx, nqy, a, n, bit);
}

/*  ############################################################################  */
/// Huffman decoding for fixed codes
///
/// Coded values range from 0-15
///
/// Huffman code values (hex):
///
///	3e, 00, 01, 08, 02, 09, 1a, 1b,
///	03, 1c, 0a, 1d, 0b, 1e, 3f, 0c
///
/// and number of bits in each code:
///
///	6,  3,  3,  4,  3,  4,  5,  5,
///	3,  5,  4,  5,  4,  5,  6,  4
///
fn input_huffman(infile: &mut Cursor<&[u8]>, b2: &mut Buffer2) -> i32 {
    // get first 3 bits to start
    let mut c: i32 = input_nbits(infile, 3, b2);

    if c < 4 {
        // this is all we need
        // return 1,2,4,8 for c=0,1,2,3
        return 1 << c;
    }

    // get the next bit
    c = input_bit(infile, b2) | (c << 1);
    if c < 13 {
        // OK, 4 bits is enough
        match c {
            8 => return 3,
            9 => return 5,
            10 => return 10,
            11 => return 12,
            12 => return 15,
            _ => (),
        }
    }

    // get yet another bit
    c = input_bit(infile, b2) | (c << 1);
    if c < 31 {
        // OK, 5 bits is enough
        match c {
            26 => return 6,
            27 => return 7,
            28 => return 9,
            29 => return 11,
            30 => return 13,
            _ => (),
        }
    }

    // need the 6th bit
    c = input_bit(infile, b2) | (c << 1);
    if c == 62 {
        0
    } else {
        14
    }
}

/*  ############################################################################  */
/// Read integer A one byte at a time from infile.
///
/// This is portable from Vax to Sun since it eliminates the need for byte-swapping.
///
/// This routine is only called to read the first 3 values
/// in the compressed file, so it doesn't have to be
/// super-efficient
#[must_use]
fn readint(infile: &mut Cursor<&[u8]>) -> i32 {
    let mut b: [u8; 4] = [0; 4];

    /*
    for i in 0..4 {
        qread(infile, &mut b[i..i], 1);
    }
    */
    qread(infile, &mut b, 4);

    let mut a: i32 = i32::from(b[0]);

    for i in 1..4 {
        a = (a << 8) + i32::from(b[i]);
    }
    a
}

/*  ############################################################################  */
/// Read integer A one byte at a time from infile.
///
/// This is portable from Vax to Sun since it eliminates the need for byte-swapping.
///
/// This routine is only called to read the first 3 values
/// in the compressed file, so it doesn't have to be
/// super-efficient
#[must_use]
fn readlonglong(infile: &mut Cursor<&[u8]>) -> i64 {
    let mut b: [u8; 8] = [0; 8];

    for i in 0..8 {
        qread(infile, &mut b[i..], 1);
    }

    let mut a: i64 = i64::from(b[0]);

    for i in 1..8 {
        a = (a << 8) + i64::from(b[i]);
    }
    a
}

/*  ############################################################################  */
/// read n bytes from file into buffer
fn qread(file: &mut Cursor<&[u8]>, buffer: &mut [u8], n: usize) {
    file.copy_to_slice(&mut buffer[0..n]);
}

/*  ############################################################################  */
/// Initialize bit input
#[must_use]
fn start_inputing_bits() -> Buffer2 {
    // Buffer starts out with no bits in it
    Buffer2 {
        buffer2: 0,    // Buffer is empty to start
        bits_to_go: 0, // with
    }
}

/*  ############################################################################  */
/// Input a bit
fn input_bit(infile: &mut Cursor<&[u8]>, b2: &mut Buffer2) -> i32 {
    if b2.bits_to_go == 0 {
        // Read the next byte if no

        b2.buffer2 = infile.get_u8() as i32;
        //b2.nextchar += 1;

        b2.bits_to_go = 8;
    }

    // Return the next bit
    b2.bits_to_go -= 1;
    (b2.buffer2 >> b2.bits_to_go) & 1
}

/*  ############################################################################  */
/// Input n bits (N must be <= 8)
fn input_nbits(infile: &mut Cursor<&[u8]>, n: usize, b2: &mut Buffer2) -> i32 {
    // AND mask for retreiving the right-most n bits
    let mask: [i32; 9] = [0, 1, 3, 7, 15, 31, 63, 127, 255];

    if b2.bits_to_go < n as i32 {
        // need another byte's worth of bits

        b2.buffer2 = (b2.buffer2 << 8) | infile.get_u8() as i32;
        b2.bits_to_go += 8;
    }

    // now pick off the first n bits
    b2.bits_to_go -= n as i32;

    // there was a slight gain in speed by replacing the following line
    //	return( (buffer2>>bits_to_go) & ((1<<n)-1) );
    (b2.buffer2 >> b2.bits_to_go) & (mask[n])
}

/*  ############################################################################  */
/// Input 4 bits
fn input_nybble(infile: &mut Cursor<&[u8]>, b2: &mut Buffer2) -> i32 {
    if b2.bits_to_go < 4 {
        // need another byte's worth of bits

        b2.buffer2 = (b2.buffer2 << 8) | infile.get_u8() as i32;
        b2.bits_to_go += 8;
    }

    // now pick off the first 4 bits
    b2.bits_to_go -= 4;

    (b2.buffer2 >> b2.bits_to_go) & 15
}

/*  ############################################################################  */
/// copy n 4-bit nybbles from infile to the lower 4 bits of array
fn input_nnybble(infile: &mut Cursor<&[u8]>, n: usize, array: &mut [u8], b2: &mut Buffer2) -> i32 {
    let mut _ii: usize;

    //  forcing byte alignment doesn;t help, and even makes it go slightly slower
    // if (bits_to_go != 8) input_nbits(infile, bits_to_go);

    if n == 1 {
        array[0] = input_nybble(infile, b2) as u8;
        return 0;
    }

    if b2.bits_to_go == 8 {
        // already have 2 full nybbles in buffer2, so
        // backspace the infile array to reuse last char

        // TODO
        infile.set_position(infile.position() - 1);
        b2.bits_to_go = 0;
    }

    // bits_to_go now has a value in the range 0 - 7.  After adding
    // another byte, bits_to_go effectively will be in range 8 - 15

    let shift1: i32 = b2.bits_to_go + 4; // shift1 will be in range 4 - 11
    let shift2: i32 = b2.bits_to_go; // shift2 will be in range 0 -  7
    let mut kk: usize = 0;

    // special case
    if b2.bits_to_go == 0 {
        for _ii in 0..(n / 2) {
            // refill the buffer with next byte

            b2.buffer2 = (b2.buffer2 << 8) | infile.get_u8() as i32;
            array[kk] = ((b2.buffer2 >> 4) & 15) as u8;
            array[kk + 1] = ((b2.buffer2) & 15) as u8; // no shift required
            kk += 2;
        }
    } else {
        for _ii in 0..(n / 2) {
            // refill the buffer with next byte
            b2.buffer2 = (b2.buffer2 << 8) | infile.get_u8() as i32;
            array[kk] = ((b2.buffer2 >> shift1) & 15) as u8;
            array[kk + 1] = ((b2.buffer2 >> shift2) & 15) as u8;
            kk += 2;
        }
    }

    let ii = n / 2;
    if ii * 2 != n {
        // have to read last odd byte
        array[n - 1] = input_nybble(infile, b2) as u8;
    }

    (b2.buffer2 >> b2.bits_to_go) & 15
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
        let mut decoder = HCDecoder::new();
        let res = decoder.read(&input, 0, &mut output).unwrap();

        assert_eq!(output.len(), 16);
        assert_eq!(res.0, 4);
        assert_eq!(res.1, 4);

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
            221, 153, 0, 0, 0, 10, 0, 0, 0, 1, 0, 0, 0, 0, 255, 255, 255, 255, 255, 255, 109, 64,
            16, 17, 0, 2, 136, 255, 191, 224, 40, 143, 251, 254, 246, 207, 253, 238, 168, 251, 53,
            238, 168, 255, 223, 247, 253, 255, 127, 223, 247, 253, 255, 127, 223, 247, 253, 255,
            127, 223, 247, 253, 255, 127, 223, 247, 253, 255, 127, 223, 247, 178, 255, 239, 251,
            217, 127, 247, 178, 253, 151, 255, 127, 222, 234, 143, 186, 163, 255, 123, 170, 62,
            234, 143, 186, 163, 238, 168, 255, 192, 100,
        ];
        let mut output: Vec<i32> = vec![0; 10];
        let mut decoder = HCDecoder::new();
        let res = decoder.read(&input, 0, &mut output).unwrap();

        assert_eq!(res.0, 10);
        assert_eq!(res.1, 1);
        assert_eq!(res.2, 0);

        assert_eq!(output.len(), 10);

        assert_eq!(output, [-1, -1, -112, -1, 9983, -28528, -112, -1, -1, -1]);
        assert_eq!(
            input,
            [
                221, 153, 0, 0, 0, 10, 0, 0, 0, 1, 0, 0, 0, 0, 255, 255, 255, 255, 255, 255, 109,
                64, 16, 17, 0, 2, 136, 255, 191, 224, 40, 143, 251, 254, 246, 207, 253, 238, 168,
                251, 53, 238, 168, 255, 223, 247, 253, 255, 127, 223, 247, 253, 255, 127, 223, 247,
                253, 255, 127, 223, 247, 253, 255, 127, 223, 247, 253, 255, 127, 223, 247, 178,
                255, 239, 251, 217, 127, 247, 178, 253, 151, 255, 127, 222, 234, 143, 186, 163,
                255, 123, 170, 62, 234, 143, 186, 163, 238, 168, 255, 192, 100
            ]
        );
    }

    #[test]
    fn test_fits_decompress_strange2() {
        let input: [u8; 106] = [
            221, 153, 0, 0, 0, 12, 0, 0, 0, 1, 0, 0, 0, 0, 255, 255, 255, 255, 255, 251, 193, 64,
            17, 16, 0, 246, 79, 151, 94, 233, 144, 42, 143, 46, 128, 170, 0, 130, 128, 42, 0, 32,
            128, 130, 0, 130, 0, 162, 191, 239, 251, 254, 255, 191, 239, 251, 254, 255, 191, 239,
            251, 254, 255, 191, 239, 251, 254, 255, 191, 239, 251, 254, 255, 191, 224, 40, 47, 46,
            189, 150, 5, 4, 1, 20, 5, 5, 229, 215, 178, 252, 182, 4, 21, 236, 191, 45, 249, 111,
            203, 96, 81, 95, 240, 0, 175, 0,
        ];
        let mut output: Vec<i32> = vec![0; 12];
        let mut decoder = HCDecoder::new();
        let res = decoder.read(&input, 0, &mut output).unwrap();

        assert_eq!(res.0, 12);
        assert_eq!(res.1, 1);
        assert_eq!(res.2, 0);

        assert_eq!(output.len(), 12);

        assert_eq!(
            output,
            [-1, -1, -9584, -28561, -112, -24321, -1, -1, -1, -9584, -28561, -112]
        );
        assert_eq!(
            input,
            [
                221, 153, 0, 0, 0, 12, 0, 0, 0, 1, 0, 0, 0, 0, 255, 255, 255, 255, 255, 251, 193,
                64, 17, 16, 0, 246, 79, 151, 94, 233, 144, 42, 143, 46, 128, 170, 0, 130, 128, 42,
                0, 32, 128, 130, 0, 130, 0, 162, 191, 239, 251, 254, 255, 191, 239, 251, 254, 255,
                191, 239, 251, 254, 255, 191, 239, 251, 254, 255, 191, 239, 251, 254, 255, 191,
                224, 40, 47, 46, 189, 150, 5, 4, 1, 20, 5, 5, 229, 215, 178, 252, 182, 4, 21, 236,
                191, 45, 249, 111, 203, 96, 81, 95, 240, 0, 175, 0
            ]
        );
    }

    #[test]
    fn test_fits_decompress_strange3() {
        let input: [u8; 140] = [
            221, 153, 0, 0, 0, 16, 0, 0, 0, 1, 0, 0, 0, 0, 255, 255, 255, 255, 255, 254, 197, 32,
            19, 17, 0, 246, 207, 253, 245, 68, 211, 250, 117, 84, 126, 153, 116, 5, 21, 31, 78,
            153, 63, 84, 106, 171, 234, 139, 170, 2, 138, 175, 170, 102, 143, 186, 39, 233, 146,
            126, 155, 100, 8, 160, 190, 153, 170, 242, 207, 253, 255, 127, 223, 247, 253, 255, 127,
            223, 247, 253, 255, 127, 223, 247, 253, 255, 127, 223, 247, 253, 255, 127, 222, 75,
            250, 167, 85, 95, 77, 183, 245, 78, 137, 2, 42, 175, 170, 109, 191, 166, 89, 8, 40,
            143, 170, 116, 104, 8, 136, 190, 137, 103, 234, 153, 39, 233, 178, 66, 136, 163, 234,
            155, 53, 244, 77, 31, 248, 0, 158, 248,
        ];
        let mut output: Vec<i32> = vec![0; 16];
        let mut decoder = HCDecoder::new();
        let res = decoder.read(&input, 0, &mut output).unwrap();

        assert_eq!(res.0, 16);
        assert_eq!(res.1, 1);
        assert_eq!(res.2, 0);

        assert_eq!(output.len(), 16);

        assert_eq!(
            output,
            [
                2570, 28560, -5778, -28528, 28816, 28816, 2671, -246, -1, -28417, -5778, -28528,
                28304, -28439, -28528, -5791
            ]
        );
        assert_eq!(
            input,
            [
                221, 153, 0, 0, 0, 16, 0, 0, 0, 1, 0, 0, 0, 0, 255, 255, 255, 255, 255, 254, 197,
                32, 19, 17, 0, 246, 207, 253, 245, 68, 211, 250, 117, 84, 126, 153, 116, 5, 21, 31,
                78, 153, 63, 84, 106, 171, 234, 139, 170, 2, 138, 175, 170, 102, 143, 186, 39, 233,
                146, 126, 155, 100, 8, 160, 190, 153, 170, 242, 207, 253, 255, 127, 223, 247, 253,
                255, 127, 223, 247, 253, 255, 127, 223, 247, 253, 255, 127, 223, 247, 253, 255,
                127, 222, 75, 250, 167, 85, 95, 77, 183, 245, 78, 137, 2, 42, 175, 170, 109, 191,
                166, 89, 8, 40, 143, 170, 116, 104, 8, 136, 190, 137, 103, 234, 153, 39, 233, 178,
                66, 136, 163, 234, 155, 53, 244, 77, 31, 248, 0, 158, 248
            ]
        );
    }

    #[test]
    fn test_fits_decompress_strange4() {
        let input: [u8; 154] = [
            221, 153, 0, 0, 0, 82, 0, 0, 0, 1, 0, 0, 0, 0, 255, 255, 255, 255, 255, 253, 245, 0,
            18, 14, 0, 246, 219, 103, 254, 246, 217, 103, 219, 101, 159, 109, 150, 125, 182, 89,
            246, 221, 50, 79, 182, 233, 209, 100, 10, 32, 136, 130, 0, 160, 160, 128, 0, 136, 35,
            221, 55, 69, 146, 91, 247, 77, 209, 100, 150, 253, 211, 117, 76, 146, 91, 126, 233,
            186, 166, 104, 150, 232, 251, 167, 84, 203, 44, 151, 79, 217, 116, 232, 183, 236, 179,
            76, 255, 223, 247, 253, 255, 127, 223, 247, 253, 255, 127, 223, 247, 253, 255, 127,
            223, 247, 253, 255, 121, 101, 183, 255, 127, 223, 247, 211, 166, 106, 137, 162, 219,
            166, 89, 108, 159, 251, 254, 255, 191, 239, 251, 254, 255, 189, 211, 102, 139, 166, 89,
            255, 128, 141, 59, 83, 9, 0,
        ];
        let mut output: Vec<i32> = vec![0; 82];
        let mut decoder = HCDecoder::new();
        let res = decoder.read(&input, 0, &mut output).unwrap();

        assert_eq!(res.0, 82);
        assert_eq!(res.1, 1);
        assert_eq!(res.2, 0);

        assert_eq!(output.len(), 82);

        assert_eq!(
            output,
            [
                -1, -1, 1, -256, -1, 0, 0, 0, 0, 0, -256, 1, -256, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 256, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, -1, -1,
                -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 255, 0, 0, 0, 0, -256, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, -256, -4097, -1,
            ]
        );
        assert_eq!(
            input,
            [
                221, 153, 0, 0, 0, 82, 0, 0, 0, 1, 0, 0, 0, 0, 255, 255, 255, 255, 255, 253, 245,
                0, 18, 14, 0, 246, 219, 103, 254, 246, 217, 103, 219, 101, 159, 109, 150, 125, 182,
                89, 246, 221, 50, 79, 182, 233, 209, 100, 10, 32, 136, 130, 0, 160, 160, 128, 0,
                136, 35, 221, 55, 69, 146, 91, 247, 77, 209, 100, 150, 253, 211, 117, 76, 146, 91,
                126, 233, 186, 166, 104, 150, 232, 251, 167, 84, 203, 44, 151, 79, 217, 116, 232,
                183, 236, 179, 76, 255, 223, 247, 253, 255, 127, 223, 247, 253, 255, 127, 223, 247,
                253, 255, 127, 223, 247, 253, 255, 121, 101, 183, 255, 127, 223, 247, 211, 166,
                106, 137, 162, 219, 166, 89, 108, 159, 251, 254, 255, 191, 239, 251, 254, 255, 189,
                211, 102, 139, 166, 89, 255, 128, 141, 59, 83, 9, 0
            ]
        );
    }

    #[test]
    fn test_fits_decompress_strange6() {
        let input: [u8; 84] = [
            221, 153, 0, 0, 0, 1, 0, 0, 0, 10, 0, 0, 0, 0, 255, 255, 255, 255, 255, 250, 185, 160,
            20, 9, 0, 246, 215, 253, 237, 95, 253, 238, 237, 123, 87, 238, 237, 127, 223, 247, 182,
            189, 221, 175, 119, 107, 221, 218, 246, 175, 221, 218, 246, 175, 254, 255, 191, 239,
            251, 219, 127, 247, 253, 255, 127, 222, 219, 246, 223, 253, 255, 127, 223, 247, 253,
            255, 127, 223, 247, 253, 255, 127, 192, 128,
        ];
        let mut output: Vec<i32> = vec![0; 10];
        let mut decoder = HCDecoder::new();
        let res = decoder.read(&input, 0, &mut output).unwrap();

        assert_eq!(res.0, 1);
        assert_eq!(res.1, 10);
        assert_eq!(res.2, 0);

        assert_eq!(output.len(), 10);

        assert_eq!(
            output,
            [-28662, -28528, 18761, 18761, 18761, 18761, 18761, 18761, -28528, -28528]
        );
        assert_eq!(
            input,
            [
                221, 153, 0, 0, 0, 1, 0, 0, 0, 10, 0, 0, 0, 0, 255, 255, 255, 255, 255, 250, 185,
                160, 20, 9, 0, 246, 215, 253, 237, 95, 253, 238, 237, 123, 87, 238, 237, 127, 223,
                247, 182, 189, 221, 175, 119, 107, 221, 218, 246, 175, 221, 218, 246, 175, 254,
                255, 191, 239, 251, 219, 127, 247, 253, 255, 127, 222, 219, 246, 223, 253, 255,
                127, 223, 247, 253, 255, 127, 223, 247, 253, 255, 127, 192, 128
            ]
        );
    }

    // Original code fails and has a bad format code
    #[test]
    fn test_fits_decompress_strange7() {
        let input: [u8; 146] = [
            221, 153, 0, 0, 0, 10, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 64, 0, 0, 0, 31, 28, 0, 246,
            255, 239, 251, 254, 255, 191, 239, 251, 254, 255, 191, 239, 251, 254, 255, 191, 239,
            251, 254, 255, 191, 239, 251, 254, 255, 191, 239, 251, 254, 255, 191, 239, 251, 254,
            255, 191, 239, 251, 254, 255, 191, 239, 251, 254, 255, 191, 239, 251, 254, 255, 191,
            239, 251, 254, 255, 191, 239, 251, 254, 255, 191, 239, 251, 254, 255, 191, 239, 251,
            254, 255, 191, 239, 251, 203, 126, 91, 242, 223, 150, 252, 183, 229, 191, 45, 249, 111,
            203, 126, 91, 242, 223, 150, 252, 183, 229, 191, 45, 249, 111, 203, 126, 91, 242, 223,
            150, 252, 183, 229, 191, 45, 249, 111, 203, 126, 91, 242, 223, 252, 0, 0,
        ];
        let mut output: Vec<i32> = vec![0; 10];
        let mut decoder = HCDecoder::new();
        let res = decoder.read(&input, 0, &mut output);

        assert!(res.is_err());
        assert!(res.unwrap_err() == DecodeError::BadFormatCode);
    }

    #[test]
    fn test_decomp_64bit_input1() {
        let input: [u8; 32] = [
            221,153,0,0,0,2,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,4,0,2,0,255,191,239,127,240,0,0
        ];
        let mut output: Vec<i64> = vec![0; 2];
        let mut decoder = HCDecoder::new();
        let res = decoder.read64(&input, 0, &mut output).unwrap();

        assert_eq!(res.0, 2);
        assert_eq!(res.1, 1);
        assert_eq!(res.2, 0);

        assert_eq!(output.len(), 2);

        assert_eq!(output, [0, 1]);
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

        hinv(&mut a[0..], nx, ny, smooth, scale).unwrap();

        assert_eq!(a, [2, 2, 1, 2, 3, 2, 7, 7, 4, 2, 2, 1, 2, 4, 25, 2]);
    }

    #[test]
    fn test_qtree_bitins() {
        let mut a = [0, 2, 0, 0];
        let nx = 6;
        let ny = 1;
        let mut b = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
        let n = 1;
        let bit = 16;
        qtree_bitins(&mut a, nx, ny, &mut b, n, bit);

        let expected = [0, 0, 0, 65536, 0, 0, 0, 0, 0, 0, 0, 0];

        assert_eq!(b, expected);
    }

    #[test]
    fn test_input_nnybble() {
        let n = 4;
        let mut array = [0, 10, 8, 0];
        let mut b2 = Buffer2 {
            bits_to_go: 5,
            buffer2: 2123985925,
        };

        let input: [u8; 140] = [
            221, 153, 0, 0, 0, 16, 0, 0, 0, 1, 0, 0, 0, 0, 255, 255, 255, 255, 255, 254, 197, 32,
            19, 17, 0, 246, 207, 253, 245, 68, 211, 250, 117, 84, 126, 153, 116, 5, 21, 31, 78,
            153, 63, 84, 106, 171, 234, 139, 170, 2, 138, 175, 170, 102, 143, 186, 39, 233, 146,
            126, 155, 100, 8, 160, 190, 153, 170, 242, 207, 253, 255, 127, 223, 247, 253, 255, 127,
            223, 247, 253, 255, 127, 223, 247, 253, 255, 127, 223, 247, 253, 255, 127, 222, 75,
            250, 167, 85, 95, 77, 183, 245, 78, 137, 2, 42, 175, 170, 109, 191, 166, 89, 8, 40,
            143, 170, 116, 104, 8, 136, 190, 137, 103, 234, 153, 39, 233, 178, 66, 136, 163, 234,
            155, 53, 244, 77, 31, 248, 0, 158, 248,
        ];
        let mut infile = Cursor::new(&input[..]);
        infile.set_position(38);

        let res = input_nnybble(&mut infile, n, &mut array, &mut b2);
        assert_eq!(res, 8);
        assert_eq!(b2.bits_to_go, 5);
        assert_eq!(b2.buffer2, 1946490143);
        assert_eq!(array, [2, 8, 10, 8]);

        //n=4  array=[0,10,8,0] bits_to_go=5, buffer2=2123985925, ---->  buffer2=1946490143  bits_to_go=5  array=[2,8,10,8]
    }
}
