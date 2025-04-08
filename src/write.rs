use crate::CODE_MAGIC;

fn ffpmsg(_m: &str) {}

/// Huffman code values and
const CODE: [i32; 16] = [
    0x3e, 0x00, 0x01, 0x08, 0x02, 0x09, 0x1a, 0x1b, 0x03, 0x1c, 0x0a, 0x1d, 0x0b, 0x1e, 0x3f, 0x0c,
];

/// Number of bits in each code
const NCODE: [i32; 16] = [6, 3, 3, 4, 3, 4, 5, 5, 3, 5, 4, 5, 4, 5, 6, 4];

#[derive(Debug)]
pub enum EncodeError {
    DataCompressionError,
}

pub struct HCEncoder {}

impl Default for HCEncoder {
    fn default() -> Self {
        Self::new()
    }
}

impl HCEncoder {
    pub fn new() -> Self {
        HCEncoder {}
    }

    /* ---------------------------------------------------------------------- */
    /// Compress the input image using the H-compress algorithm
    ///
    /// # Arguments
    ///
    /// * `a` - input image array
    /// * `nx` - size of X axis of image
    /// * `ny` - size of Y axis of image
    /// * `scale` - quantization scale factor. Larger values results in more (lossy) compression. 0 does lossless compression
    /// * `output` - pre-allocated array to hold the output compressed stream of bytes
    ///
    /// NOTE: the nx and ny dimensions as defined within this code are reversed from
    /// the usual FITS notation.  ny is the fastest varying dimension, which is
    /// usually considered the X axis in the FITS image display
    ///
    /// NOTE: The input data i32 is a i16 number padded to a i32 number.
    pub fn write(
        &self,
        a: &mut [i32],
        ny: usize,
        nx: usize,
        scale: i32,
        output: &mut Vec<u8>,
    ) -> Result<(), EncodeError> {
        // H-transform
        htrans(a, nx, ny)?;

        // digitize
        digitize(a, nx, ny, scale);

        // encode and write to output array
        // input value is the allocated size of the array
        encode(output, a, nx, ny, scale)
    }

    /* ---------------------------------------------------------------------- */
    /// Compress the input image using the H-compress algorithm
    ///
    /// # Arguments
    ///
    /// * `a` - input image array
    /// * `nx` - size of X axis of image
    /// * `ny` - size of Y axis of image
    /// * `scale` - quantization scale factor. Larger values results in more (lossy) compression. 0 does lossless compression
    /// * `output` - pre-allocated array to hold the output compressed stream of bytes
    ///
    /// NOTE: the nx and ny dimensions as defined within this code are reversed from
    /// the usual FITS notation.  ny is the fastest varying dimension, which is
    /// usually considered the X axis in the FITS image display
    pub fn write64(
        &self,
        a: &mut [i64],
        ny: usize,
        nx: usize,
        scale: i32,
        output: &mut Vec<u8>,
    ) -> Result<(), EncodeError> {
        // H-transform
        htrans64(a, nx, ny)?;

        // digitize
        digitize64(a, nx, ny, scale);

        // encode and write to output array
        // input value is the allocated size of the array
        encode64(output, a, nx, ny, scale)
    }
}

/* ######################################################################### */
/// H-transform of NX x NY integer imag
///
/// # Arguments
///
/// * `a` - input image array
/// * `nx` - size of X axis of image
/// * `ny` - size of Y axis of image
///
fn htrans(a: &mut [i32], nx: usize, ny: usize) -> Result<(), EncodeError> {
    let mut h0: i32;
    let mut hx: i32;
    let mut hy: i32;
    let mut hc: i32;

    let mut oddx: usize;
    let mut oddy: usize;

    let mut s10: usize;
    let mut s00: usize;

    let nmax: usize = if nx > ny { nx } else { ny };

    // log2n is log2 of max(nx,ny) rounded up to next power of 2
    let mut log2n: i32 = ((nmax as f64).ln() / 2.0_f64.ln() + 0.5) as i32;

    if nmax > (1 << log2n) {
        log2n += 1;
    }

    // get temporary storage for shuffling elements
    let mut tmp: Vec<i32> = Vec::new();
    if tmp.try_reserve_exact(nmax.div_ceil(2)).is_err() {
        ffpmsg("htrans: insufficient memory");
        return Err(EncodeError::DataCompressionError);
    } else {
        tmp.resize(nmax.div_ceil(2), 0);
    }

    // set up rounding and shifting masks
    let mut shift: i32 = 0;
    let mut mask: i32 = -2;
    let mut mask2: i32 = mask << 1;
    let mut prnd: i32 = 1;
    let mut prnd2: i32 = prnd << 1;
    let mut nrnd2: i32 = prnd2 - 1;

    /*
     do log2n reductions

     We're indexing a as a 2-D array with dimensions (nx,ny).

    */

    let mut nxtop: usize = nx;
    let mut nytop: usize = ny;

    for _k in 0..log2n {
        oddx = nxtop % 2;
        oddy = nytop % 2;

        for i in (0..(nxtop - oddx)).step_by(2) {
            s00 = i * ny; // s00 is index of a[i,j]
            s10 = s00 + ny; // s10 is index of a[i+1,j]

            for _j in (0..(nytop - oddy)).step_by(2) {
                // Divide h0,hx,hy,hc by 2 (1 the first time through).
                h0 = (a[s10 + 1] + a[s10] + a[s00 + 1] + a[s00]) >> shift;
                hx = (a[s10 + 1] + a[s10] - a[s00 + 1] - a[s00]) >> shift;
                hy = (a[s10 + 1] - a[s10] + a[s00 + 1] - a[s00]) >> shift;
                hc = (a[s10 + 1] - a[s10] - a[s00 + 1] + a[s00]) >> shift;

                /*
                Throw away the 2 bottom bits of h0, bottom bit of hx,hy.
                To get rounding to be same for positive and negative
                numbers, nrnd2 = prnd2 - 1.
                 */
                a[s10 + 1] = hc;
                a[s10] = (if hx >= 0 { hx + prnd } else { hx }) & mask;
                a[s00 + 1] = (if hy >= 0 { hy + prnd } else { hy }) & mask;
                a[s00] = (if h0 >= 0 { h0 + prnd2 } else { h0 + nrnd2 }) & mask2;
                s00 += 2;
                s10 += 2;
            }

            if oddy > 0 {
                /*
                do last element in row if row length is odd
                s00+1, s10+1 are off edge
                 */
                h0 = (a[s10] + a[s00]) << (1 - shift);
                hx = (a[s10] - a[s00]) << (1 - shift);
                a[s10] = (if hx >= 0 { hx + prnd } else { hx }) & mask;
                a[s00] = (if h0 >= 0 { h0 + prnd2 } else { h0 + nrnd2 }) & mask2;
                s00 += 1;
                s10 += 1;
            }
        }

        if oddx > 0 {
            /*
            do last row if column length is odd
            s10, s10+1 are off edge
             */
            s00 = (nxtop - 1) * ny; // i*ny;
            for _j in (0..(nytop - oddy)).step_by(2) {
                h0 = (a[s00 + 1] + a[s00]) << (1 - shift);
                hy = (a[s00 + 1] - a[s00]) << (1 - shift);
                a[s00 + 1] = (if hy >= 0 { hy + prnd } else { hy }) & mask;
                a[s00] = (if h0 >= 0 { h0 + prnd2 } else { h0 + nrnd2 }) & mask2;
                s00 += 2;
            }
            if oddy > 0 {
                /*
                do corner element if both row and column lengths are odd
                s00+1, s10, s10+1 are off edge
                 */
                h0 = a[s00] << (2 - shift);
                a[s00] = (if h0 >= 0 { h0 + prnd2 } else { h0 + nrnd2 }) & mask2;
            }
        }

        // Now shuffle in each dimension to group coefficients by order
        for i in 0..nxtop {
            shuffle(&mut a[ny * i..], nytop, 1, &mut tmp);
        }

        for j in 0..nytop {
            shuffle(&mut a[j..], nxtop, ny, &mut tmp);
        }

        // Image size reduced by 2 (round up if odd)
        nxtop = (nxtop + 1) >> 1;
        nytop = (nytop + 1) >> 1;

        // Divisor doubles after first reduction
        shift = 1;

        // Masks, rounding values double after each iteration
        mask = mask2;
        prnd = prnd2;
        mask2 <<= 1;
        prnd2 <<= 1;
        nrnd2 = prnd2 - 1;
    }

    Ok(())
}

/* ######################################################################### */
/// H-transform of NX x NY integer imag
///
/// # Arguments
///
/// * `a` - input image array
/// * `nx` - size of X axis of image
/// * `ny` - size of Y axis of image
///
fn htrans64(a: &mut [i64], nx: usize, ny: usize) -> Result<(), EncodeError> {
    let mut h0: i64;
    let mut hx: i64;
    let mut hy: i64;
    let mut hc: i64;

    let mut oddx: usize;
    let mut oddy: usize;

    let mut s10: usize;
    let mut s00: usize;

    let nmax: usize = if nx > ny { nx } else { ny };

    // log2n is log2 of max(nx,ny) rounded up to next power of 2
    let mut log2n: i32 = ((nmax as f64).ln() / 2.0_f64.ln() + 0.5) as i32;

    if nmax > (1 << log2n) {
        log2n += 1;
    }

    // get temporary storage for shuffling elements
    let mut tmp: Vec<i64> = Vec::new();
    if tmp.try_reserve_exact(nmax.div_ceil(2)).is_err() {
        ffpmsg("htrans: insufficient memory");
        return Err(EncodeError::DataCompressionError);
    } else {
        tmp.resize(nmax.div_ceil(2), 0);
    }

    // set up rounding and shifting masks
    let mut shift: i32 = 0;
    let mut mask: i64 = -2;
    let mut mask2: i64 = mask << 1;
    let mut prnd: i64 = 1;
    let mut prnd2: i64 = prnd << 1;
    let mut nrnd2: i64 = prnd2 - 1;

    /*
     do log2n reductions

     We're indexing a as a 2-D array with dimensions (nx,ny).

    */

    let mut nxtop: usize = nx;
    let mut nytop: usize = ny;

    for _k in 0..log2n {
        oddx = nxtop % 2;
        oddy = nytop % 2;

        for i in (0..(nxtop - oddx)).step_by(2) {
            s00 = i * ny; // s00 is index of a[i,j]
            s10 = s00 + ny; // s10 is index of a[i+1,j]

            for _j in (0..(nytop - oddy)).step_by(2) {
                // Divide h0,hx,hy,hc by 2 (1 the first time through).
                h0 = (a[s10 + 1] + a[s10] + a[s00 + 1] + a[s00]) >> shift;
                hx = (a[s10 + 1] + a[s10] - a[s00 + 1] - a[s00]) >> shift;
                hy = (a[s10 + 1] - a[s10] + a[s00 + 1] - a[s00]) >> shift;
                hc = (a[s10 + 1] - a[s10] - a[s00 + 1] + a[s00]) >> shift;

                /*
                Throw away the 2 bottom bits of h0, bottom bit of hx,hy.
                To get rounding to be same for positive and negative
                numbers, nrnd2 = prnd2 - 1.
                 */
                a[s10 + 1] = hc;
                a[s10] = (if hx >= 0 { hx + prnd } else { hx }) & mask;
                a[s00 + 1] = (if hy >= 0 { hy + prnd } else { hy }) & mask;
                a[s00] = (if h0 >= 0 { h0 + prnd2 } else { h0 + nrnd2 }) & mask2;
                s00 += 2;
                s10 += 2;
            }

            if oddy > 0 {
                /*
                do last element in row if row length is odd
                s00+1, s10+1 are off edge
                 */
                h0 = (a[s10] + a[s00]) << (1 - shift);
                hx = (a[s10] - a[s00]) << (1 - shift);
                a[s10] = (if hx >= 0 { hx + prnd } else { hx }) & mask;
                a[s00] = (if h0 >= 0 { h0 + prnd2 } else { h0 + nrnd2 }) & mask2;
                s00 += 1;
                s10 += 1;
            }
        }

        if oddx > 0 {
            /*
            do last row if column length is odd
            s10, s10+1 are off edge
             */
            s00 = (nxtop - 1) * ny; // i*ny;
            for _j in (0..(nytop - oddy)).step_by(2) {
                h0 = (a[s00 + 1] + a[s00]) << (1 - shift);
                hy = (a[s00 + 1] - a[s00]) << (1 - shift);
                a[s00 + 1] = (if hy >= 0 { hy + prnd } else { hy }) & mask;
                a[s00] = (if h0 >= 0 { h0 + prnd2 } else { h0 + nrnd2 }) & mask2;
                s00 += 2;
            }
            if oddy > 0 {
                /*
                do corner element if both row and column lengths are odd
                s00+1, s10, s10+1 are off edge
                 */
                h0 = a[s00] << (2 - shift);
                a[s00] = (if h0 >= 0 { h0 + prnd2 } else { h0 + nrnd2 }) & mask2;
            }
        }

        // Now shuffle in each dimension to group coefficients by order
        for i in 0..nxtop {
            shuffle64(&mut a[ny * i..], nytop, 1, &mut tmp);
        }

        for j in 0..nytop {
            shuffle64(&mut a[j..], nxtop, ny, &mut tmp);
        }

        // Image size reduced by 2 (round up if odd)
        nxtop = (nxtop + 1) >> 1;
        nytop = (nytop + 1) >> 1;

        // Divisor doubles after first reduction
        shift = 1;

        // Masks, rounding values double after each iteration
        mask = mask2;
        prnd = prnd2;
        mask2 <<= 1;
        prnd2 <<= 1;
        nrnd2 = prnd2 - 1;
    }

    Ok(())
}

/* ######################################################################### */
/// Digitize H-transform, round to multiple of scale
///
/// # Arguments
///
/// * `a` - Array to round
/// * `nx` - Number of elements in column
/// * `ny` - Number of elements in row
/// * `scale` - Scale factor
///
fn digitize(a: &mut [i32], nx: usize, ny: usize, scale: i32) {
    if scale <= 1 {
        return;
    };

    let d: i32 = (scale + 1) / 2 - 1;

    if d == 0 {
        for p in 0..(nx * ny) {
            a[p] /= scale;
        }
    } else {
        for p in 0..(nx * ny) {
            if a[p] > 0 {
                a[p] = (a[p] + d) / scale;
            } else {
                a[p] = (a[p] - d) / scale;
            }
        }
    }
}

/* ######################################################################### */
/// Digitize H-transform, round to multiple of scale
///
/// # Arguments
///
/// * `a` - Array to round
/// * `nx` - Number of elements in column
/// * `ny` - Number of elements in row
/// * `scale` - Scale factor
///
fn digitize64(a: &mut [i64], nx: usize, ny: usize, scale: i32) {
    if scale <= 1 {
        return;
    };

    let scale = scale as i64;

    let d: i64 = (scale + 1) / 2 - 1;

    if d == 0 {
        for p in 0..(nx * ny) {
            a[p] /= scale;
        }
    } else {
        for p in 0..(nx * ny) {
            if a[p] > 0 {
                a[p] = (a[p] + d) / scale;
            } else {
                a[p] = (a[p] - d) / scale;
            }
        }
    }
}

/* ######################################################################### */
/// Shuffles array by copying odd elements to tmp
///
/// # Arguments
///
/// * `a` - Array to shuffle
/// * `n` - Number of elements to shuffle
/// * `n2` - Second dimension
/// * `tmp` - Scratch storage
///
fn shuffle(a: &mut [i32], n: usize, n2: usize, tmp: &mut [i32]) {
    let mut pt: usize = 0;
    let mut p1: usize = n2;

    for _i in (1..n).step_by(2) {
        tmp[pt] = a[p1];
        pt += 1;
        p1 += n2 + n2;
    }

    // Compress even elements into first half of A
    p1 = n2;
    let mut p2: usize = n2 + n2;
    for _i in (2..n).step_by(2) {
        a[p1] = a[p2];
        p1 += n2;
        p2 += n2 + n2;
    }

    // Put odd elements into 2nd half
    pt = 0;
    for _i in (1..n).step_by(2) {
        a[p1] = tmp[pt];
        p1 += n2;
        pt += 1;
    }
}

/* ######################################################################### */
/// Shuffles array by copying odd elements to tmp
///
/// # Arguments
///
/// * `a` - Array to shuffle
/// * `n` - Number of elements to shuffle
/// * `n2` - Second dimension
/// * `tmp` - Scratch storage
///
fn shuffle64(a: &mut [i64], n: usize, n2: usize, tmp: &mut [i64]) {
    let mut pt: usize = 0;
    let mut p1: usize = n2;

    for _i in (1..n).step_by(2) {
        tmp[pt] = a[p1];
        pt += 1;
        p1 += n2 + n2;
    }

    // Compress even elements into first half of A
    p1 = n2;
    let mut p2: usize = n2 + n2;
    for _i in (2..n).step_by(2) {
        a[p1] = a[p2];
        p1 += n2;
        p2 += n2 + n2;
    }

    // Put odd elements into 2nd half
    pt = 0;
    for _i in (1..n).step_by(2) {
        a[p1] = tmp[pt];
        p1 += n2;
        pt += 1;
    }
}

/* ######################################################################### */
/// Performs H-Compress encoding
///
/// # Arguments
///
/// * `outfile` - Output buffer
/// * `a` - input H-transform array (nx,ny)
/// * `nx` - size of H-transform array
/// * `ny` - size of H-transform array
/// * `scale` - scale factor for digitization
///
fn encode(
    outfile: &mut Vec<u8>,
    a: &mut [i32],
    nx: usize,
    ny: usize,
    scale: i32,
) -> Result<(), EncodeError> {
    let mut nbitplanes: [u8; 3] = [0; 3];
    let mut vmax: [i32; 3] = [0; 3];

    let nel = nx * ny;

    // write magic value
    qwrite(outfile, &CODE_MAGIC, 2);
    writeint(outfile, nx as i32); /* size of image */
    writeint(outfile, ny as i32);
    writeint(outfile, scale); /* scale factor for digitization */

    /*
    write first value of A (sum of all pixels -- the only value
    which does not compress well)
     */
    writelonglong(outfile, i64::from(a[0]));

    a[0] = 0;

    // allocate array for sign bits and save values, 8 per byte (initialize to all zeros)
    let mut signbits: Vec<u8> = Vec::new();
    if signbits.try_reserve_exact(nel.div_ceil(8)).is_err() {
        ffpmsg("encode: insufficient memory");
        return Err(EncodeError::DataCompressionError);
    } else {
        signbits.resize(nel.div_ceil(8), 0);
    }

    let mut nsign = 0;
    let mut bits_to_go = 8;

    for i in 0..nel {
        if a[i] > 0 {
            // positive element, put zero at end of buffer

            signbits[nsign] <<= 1;
            bits_to_go -= 1;
        } else if a[i] < 0 {
            // negative element, shift in a one

            signbits[nsign] <<= 1;
            signbits[nsign] |= 1;
            bits_to_go -= 1;

            //replace a by absolute value
            a[i] = -a[i];
        }

        if bits_to_go == 0 {
            // filled up this byte, go to the next one

            bits_to_go = 8;
            nsign += 1;
        }
    }

    if bits_to_go != 8 {
        /*
        some bits in last element
        move bits in last byte to bottom and increment nsign
         */
        signbits[nsign] <<= bits_to_go;
        nsign += 1;
    }

    /*
    calculate number of bit planes for 3 quadrants

    quadrant 0=bottom left, 1=bottom right or top left, 2=top right,
     */
    vmax.iter_mut().for_each(|x| *x = 0);

    // get maximum absolute value in each quadrant
    let nx2 = nx.div_ceil(2);
    let ny2 = ny.div_ceil(2);
    let mut j = 0; /* column counter	*/
    let mut k = 0; /* row counter		*/
    for i in 0..nel {
        let q = usize::from(j >= ny2) + usize::from(k >= nx2);
        if vmax[q] < a[i] {
            vmax[q] = a[i];
        }

        j += 1;
        if j >= ny {
            j = 0;
            k += 1;
        }
    }

    // now calculate number of bits for each quadrant

    for q in 0..3 {
        nbitplanes[q] = 0;
        while vmax[q] > 0 {
            vmax[q] >>= 1;
            nbitplanes[q] += 1;
        }
    }

    // write nbitplanes
    if qwrite(outfile, &nbitplanes, nbitplanes.len()) == 0 {
        ffpmsg("encode: output buffer too small");
        return Err(EncodeError::DataCompressionError);
    }

    // write coded array
    let stat = do_encode(outfile, a, nx, ny, nbitplanes);

    // write sign bits
    if nsign > 0 && qwrite(outfile, &signbits, nsign) == 0 {
        ffpmsg("encode: output buffer too small");
        return Err(EncodeError::DataCompressionError);
    }

    stat
}

/* ######################################################################### */
/// Performs H-Compress encoding
///
/// # Arguments
///
/// * `outfile` - Output buffer
/// * `a` - input H-transform array (nx,ny)
/// * `nx` - size of H-transform array
/// * `ny` - size of H-transform array
/// * `scale` - scale factor for digitization
///
fn encode64(
    outfile: &mut Vec<u8>,
    a: &mut [i64],
    nx: usize,
    ny: usize,
    scale: i32,
) -> Result<(), EncodeError> {
    let mut nbitplanes: [u8; 3] = [0; 3];
    let mut vmax: [i64; 3] = [0; 3];

    let nel = nx * ny;

    // write magic value
    qwrite(outfile, &CODE_MAGIC, 2);
    writeint(outfile, nx as i32); /* size of image */
    writeint(outfile, ny as i32);
    writeint(outfile, scale); /* scale factor for digitization */

    /*
    write first value of A (sum of all pixels -- the only value
    which does not compress well)
     */
    writelonglong(outfile, a[0]);

    a[0] = 0;

    // allocate array for sign bits and save values, 8 per byte (initialize to all zeros)
    let mut signbits: Vec<u8> = Vec::new();
    if signbits.try_reserve_exact(nel.div_ceil(8)).is_err() {
        ffpmsg("encode64: insufficient memory");
        return Err(EncodeError::DataCompressionError);
    } else {
        signbits.resize(nel.div_ceil(8), 0);
    }

    let mut nsign = 0;
    let mut bits_to_go = 8;

    for i in 0..nel {
        if a[i] > 0 {
            // positive element, put zero at end of buffer

            signbits[nsign] <<= 1;
            bits_to_go -= 1;
        } else if a[i] < 0 {
            // negative element, shift in a one

            signbits[nsign] <<= 1;
            signbits[nsign] |= 1;
            bits_to_go -= 1;

            //replace a by absolute value
            a[i] = -a[i];
        }

        if bits_to_go == 0 {
            // filled up this byte, go to the next one

            bits_to_go = 8;
            nsign += 1;
        }
    }

    if bits_to_go != 8 {
        /*
        some bits in last element
        move bits in last byte to bottom and increment nsign
         */
        signbits[nsign] <<= bits_to_go;
        nsign += 1;
    }

    /*
    calculate number of bit planes for 3 quadrants

    quadrant 0=bottom left, 1=bottom right or top left, 2=top right,
     */
    vmax.iter_mut().for_each(|x| *x = 0);

    // get maximum absolute value in each quadrant
    let nx2 = nx.div_ceil(2);
    let ny2 = ny.div_ceil(2);
    let mut j = 0; /* column counter	*/
    let mut k = 0; /* row counter		*/
    for i in 0..nel {
        let q = usize::from(j >= ny2) + usize::from(k >= nx2);
        if vmax[q] < a[i] {
            vmax[q] = a[i];
        }

        j += 1;
        if j >= ny {
            j = 0;
            k += 1;
        }
    }

    // now calculate number of bits for each quadrant

    for q in 0..3 {
        nbitplanes[q] = 0;
        while vmax[q] > 0 {
            vmax[q] >>= 1;
            nbitplanes[q] += 1;
        }
    }

    // write nbitplanes
    if qwrite(outfile, &nbitplanes, nbitplanes.len()) == 0 {
        ffpmsg("encode64: output buffer too small");
        return Err(EncodeError::DataCompressionError);
    }

    // write coded array
    let stat = do_encode64(outfile, a, nx, ny, nbitplanes);

    // write sign bits
    if nsign > 0 && qwrite(outfile, &signbits, nsign) == 0 {
        ffpmsg("encode64: output buffer too small");
        return Err(EncodeError::DataCompressionError);
    }

    stat
}

/* ######################################################################### */
/// Write n bytes from buffer into file.
///
/// # Arguments
///
/// * `file` - Output data stream
/// * `buffer` - Buffer to write
/// * `n` - Number of bytes to write
///
/// # Returns
///
/// Number of bytes written (=n) if successful, <=0 if not
///
fn qwrite(file: &mut Vec<u8>, buffer: &[u8], n: usize) -> usize {
    file.extend_from_slice(&buffer[0..n]);

    n
}

/// Write integer A one byte at a time to outfile.
///
/// This is portable from Vax to Sun since it eliminates the need for byte-swapping.
///
/// # Arguments
///
/// * `outfile` - Output data stream
/// * `a` - Integer to write
///
fn writeint(outfile: &mut Vec<u8>, a: i32) {
    let mut b: [u8; 4] = [0; 4];
    let mut a = a;

    for i in (0..4).rev() {
        b[i] = a as u8;
        a >>= 8;
    }

    qwrite(outfile, &b, 4);
}

/* ######################################################################### */
/// Write integer A one byte at a time to outfile.
///
/// This is portable from Vax to Sun since it eliminates the need for byte-swapping.
///
/// # Arguments
///
/// * `outfile` - Output data stream
/// * `a` - Integer to write
///
fn writelonglong(outfile: &mut Vec<u8>, a: i64) {
    let mut a = a;
    let mut b: [u8; 8] = [0; 8];

    for i in (0..8).rev() {
        b[i] = (a & 0x0000_00ff) as u8;
        a >>= 8;
    }

    for i in 0..8 {
        qwrite(outfile, &b[i..], 1);
    }
}

/// Encode 2-D array and write stream of characters on outfile
///
/// # Arguments
///
/// * `outfile` - Output data stream
/// * `a` - Array of values to encode
/// * `nx` - Array dimensions [nx][ny]
/// * `ny` - Array dimensions [nx][ny]
/// * `nbitplanes` - Number of bit planes in quadrants
///
/// Note: This version assumes that A is positive.
fn do_encode(
    outfile: &mut Vec<u8>,
    a: &[i32],
    nx: usize,
    ny: usize,
    nbitplanes: [u8; 3],
) -> Result<(), EncodeError> {
    let nx2: usize = nx.div_ceil(2);
    let ny2: usize = ny.div_ceil(2);

    // Initialize bit output
    let mut buffer2 = start_outputing_bits();

    //write out the bit planes for each quadrant
    qtree_encode(outfile, a, ny, nx2, ny2, nbitplanes[0] as i32, &mut buffer2)?;

    qtree_encode(
        outfile,
        &a[ny2..],
        ny,
        nx2,
        ny / 2,
        i32::from(nbitplanes[1]),
        &mut buffer2,
    )?;

    qtree_encode(
        outfile,
        &a[ny * nx2..],
        ny,
        nx / 2,
        ny2,
        i32::from(nbitplanes[1]),
        &mut buffer2,
    )?;

    let mut stat = Ok(());

    if ny * nx2 + ny2 < a.len() {
        stat = qtree_encode(
            outfile,
            &a[(ny * nx2 + ny2)..],
            ny,
            nx / 2,
            ny / 2,
            i32::from(nbitplanes[2]),
            &mut buffer2,
        );
    }

    // Add zero as an EOF symbol
    output_nybble(outfile, 0, &mut buffer2);
    done_outputing_bits(outfile, &mut buffer2);

    stat
}

/// Encode 2-D array and write stream of characters on outfile
///
/// # Arguments
///
/// * `outfile` - Output data stream
/// * `a` - Array of values to encode
/// * `nx` - Array dimensions [nx][ny]
/// * `ny` - Array dimensions [nx][ny]
/// * `nbitplanes` - Number of bit planes in quadrants
///
/// Note: This version assumes that A is positive.
fn do_encode64(
    outfile: &mut Vec<u8>,
    a: &[i64],
    nx: usize,
    ny: usize,
    nbitplanes: [u8; 3],
) -> Result<(), EncodeError> {
    let nx2: usize = nx.div_ceil(2);
    let ny2: usize = ny.div_ceil(2);

    // Initialize bit output
    let mut buffer2 = start_outputing_bits();

    //write out the bit planes for each quadrant
    qtree_encode64(outfile, a, ny, nx2, ny2, nbitplanes[0] as i32, &mut buffer2)?;

    qtree_encode64(
        outfile,
        &a[ny2..],
        ny,
        nx2,
        ny / 2,
        i32::from(nbitplanes[1]),
        &mut buffer2,
    )?;

    qtree_encode64(
        outfile,
        &a[ny * nx2..],
        ny,
        nx / 2,
        ny2,
        i32::from(nbitplanes[1]),
        &mut buffer2,
    )?;

    let mut stat = Ok(());

    if ny * nx2 + ny2 < a.len() {
        stat = qtree_encode64(
            outfile,
            &a[(ny * nx2 + ny2)..],
            ny,
            nx / 2,
            ny / 2,
            i32::from(nbitplanes[2]),
            &mut buffer2,
        );
    }

    // Add zero as an EOF symbol
    output_nybble(outfile, 0, &mut buffer2);
    done_outputing_bits(outfile, &mut buffer2);

    stat
}

/* ######################################################################### */
/* INITIALIZE FOR BIT OUTPUT */
pub struct Buffer3 {
    pub bitbuffer: i32,
    pub bits_to_go: i32,
}

pub struct Buffer2 {
    pub buffer2: i32,
    pub bits_to_go2: i32,
    pub bitcount: usize,
}

#[must_use]
fn start_outputing_bits() -> Buffer2 {
    Buffer2 {
        buffer2: 0,     /* Buffer is empty to start	*/
        bits_to_go2: 8, /* with				*/
        bitcount: 0,
    }
}

/// Output N bits (N must be <= 24)
///
/// # Arguments
///
/// * `outfile` - Output data stream
/// * `bits` - Bits to write
/// * `n` - Number of bits to write
/// * `b2` - Buffer to write to
///
fn output_nbits(outfile: &mut Vec<u8>, bits: i32, n: usize, b2: &mut Buffer2) {
    // AND mask for the right-most n bits
    let mask: [i32; 9] = [0, 1, 3, 7, 15, 31, 63, 127, 255];

    // insert bits at end of buffer
    b2.buffer2 <<= n;
    b2.buffer2 |= bits & mask[n];
    b2.bits_to_go2 -= n as i32;
    while b2.bits_to_go2 <= 0 {
        // buffer full, put out top 8 bits

        // ((lbitbuffer>>(-bits_to_go2)) & 0xff)
        outfile.push((b2.buffer2 >> (-b2.bits_to_go2)) as u8);
        b2.bits_to_go2 += 8;
    }
    b2.bitcount += n;
}

/* ######################################################################### */
/// OUTPUT a 4 bit nybble
///
/// # Arguments
///
/// * `outfile` - Output data stream
/// * `bits` - Bits to write
/// * `buffer` - Buffer to write to
///
fn output_nybble(outfile: &mut Vec<u8>, bits: i32, buffer: &mut Buffer2) {
    output_nbits(outfile, bits, 4, buffer)
}

/*  ############################################################################  */
/// OUTPUT array of 4 BITS
///
/// Pack the 4 lower bits in each element of the array into the outfile array
///
/// # Arguments
///
/// * `outfile` - Output data stream
/// * `n` - Number of elements in array
/// * `array` - Array of bits to write
///
fn output_nnybble(outfile: &mut Vec<u8>, n: usize, array: &[u8], b2: &mut Buffer2) {
    let mut kk = 0;

    if n == 1 {
        output_nybble(outfile, i32::from(array[0]), b2);
        return;
    }

    if b2.bits_to_go2 <= 4 {
        // just room for 1 nybble; write it out separately
        output_nybble(outfile, i32::from(array[0]), b2);
        kk += 1; /* index to next array element */

        if n == 2 {
            // only 1 more nybble to write out
            output_nybble(outfile, i32::from(array[1]), b2);
            return;
        }
    }

    // bits_to_go2 is now in the range 5 - 8
    let shift = 8 - b2.bits_to_go2;

    // now write out pairs of nybbles; this does not affect value of bits_to_go2
    let jj = (n - kk) / 2;

    if b2.bits_to_go2 == 8 {
        // special case if nybbles are aligned on byte boundary
        // this actually seems to make very little difference in speed
        b2.buffer2 = 0;
        for _ii in 0..jj {
            outfile.push(((array[kk] & 15) << 4) | (array[kk + 1] & 15));
            kk += 2;
        }
    } else {
        for _ii in 0..jj {
            b2.buffer2 =
                (b2.buffer2 << 8) | (((array[kk] & 15) << 4) | (array[kk + 1] & 15)) as i32;
            kk += 2;

            // buffer2 full, put out top 8 bits

            outfile.push(((b2.buffer2 >> shift) & 0xff) as u8);
        }
    }

    b2.bitcount += 8 * (jj - 1);

    // write out last odd nybble, if present
    if kk != n {
        output_nybble(outfile, i32::from(array[n - 1]), b2);
    }
}

/* ######################################################################### */
/// FLUSH OUT THE LAST BITS
///
/// # Arguments
///
/// * `outfile` - Output data stream
/// * `buffer` - Buffer to write to
fn done_outputing_bits(outfile: &mut Vec<u8>, buffer: &mut Buffer2) {
    if buffer.bits_to_go2 < 8 {
        outfile.push((buffer.buffer2 << buffer.bits_to_go2) as u8);

        // count the garbage bits too
        buffer.bitcount += buffer.bits_to_go2 as usize;
    }
}

/// Encode values in quadrant of 2-D array using binary quadtree coding for each bit plane.  
///
/// Assumes array is positive.
///
/// # Arguments
///
/// * `outfile` - Output data stream
/// * `a` - Array of values to encode
/// * `n` - Number of elements in array
/// * `nqx` - Number of elements in X axis
/// * `nqy` - Number of elements in Y axis
/// * `nbitplanes` - Number of bit planes to encode
/// * `buffer2` - Buffer to write to
///
fn qtree_encode(
    outfile: &mut Vec<u8>,
    a: &[i32],
    n: usize,
    nqx: usize,
    nqy: usize,
    nbitplanes: i32,
    buffer2: &mut Buffer2,
) -> Result<(), EncodeError> {
    let mut b: usize;
    let _k: i32;

    let mut nx: usize;
    let mut ny: usize;

    let mut b3: Buffer3 = Buffer3 {
        bitbuffer: 0,
        bits_to_go: 0,
    };

    // log2n is log2 of max(nqx,nqy) rounded up to next power of 2
    let nqmax: usize = if nqx > nqy { nqx } else { nqy };
    let mut log2n: i32 = ((nqmax as f32).ln() / 2.0_f32.ln() + 0.5) as i32;

    if nqmax > (1 << log2n) {
        log2n += 1;
    }

    // initialize buffer point, max buffer size
    let nqx2: usize = nqx.div_ceil(2);
    let nqy2: usize = nqy.div_ceil(2);
    let bmax: usize = (nqx2 * nqy2 + 2) / 2;

    /*
    We're indexing A as a 2-D array with dimensions (nqx,nqy).
    Scratch is 2-D with dimensions (nqx/2,nqy/2) rounded up.
    Scr1 is used to store first level of quadtree in case direct
    coding is needed.
    Buffer is used to store string of codes for output.
     */
    let mut scratch: Vec<u8> = Vec::new();
    if scratch.try_reserve_exact(2 * bmax).is_err() {
        ffpmsg("qtree_encode: insufficient memory");
        return Err(EncodeError::DataCompressionError);
    } else {
        scratch.resize(2 * bmax, 0);
    }

    let mut buffer: Vec<u8> = Vec::new();
    if buffer.try_reserve_exact(bmax).is_err() {
        ffpmsg("qtree_encode: insufficient memory");
        return Err(EncodeError::DataCompressionError);
    } else {
        buffer.resize(bmax, 0);
    }

    // now encode each bit plane, starting with the top
    for bit in (0..(nbitplanes as usize)).rev() {
        // initialize bit buffer
        b = 0;
        b3.bitbuffer = 0;
        b3.bits_to_go = 0;

        // on first pass copy A to scr1 array
        qtree_onebit(a, n, nqx, nqy, &mut scratch, bit);
        nx = (nqx + 1) >> 1;
        ny = (nqy + 1) >> 1;

        // copy non-zero values to output buffer, which will be written in reverse order
        if bufcopy(&scratch, nx * ny, &mut buffer, &mut b, bmax, &mut b3) {
            //quadtree is expanding data, change warning code and just fill buffer with bit-map
            write_bdirect(outfile, a, n, nqx, nqy, &mut scratch, bit, buffer2);
            continue;
        }

        // do log2n reductions
        let mut is_continue = false;
        for _ in 1..log2n {
            qtree_reduce(&mut scratch, ny, nx, ny);
            nx = (nx + 1) >> 1;
            ny = (ny + 1) >> 1;
            if bufcopy(&scratch, nx * ny, &mut buffer, &mut b, bmax, &mut b3) {
                // nx=1, ny=1, b=1??, scratch=[10,8,0,0], bmax=2, b3=0/1
                write_bdirect(outfile, a, n, nqx, nqy, &mut scratch, bit, buffer2);
                is_continue = true;
                break; // Break this loop and continue with next
            }
        }

        if is_continue {
            continue;
        }

        // OK, we've got the code in buffer
        // Write quadtree warning code, then write buffer in reverse order
        output_nybble(outfile, 0xF, buffer2);

        if b == 0 {
            if b3.bits_to_go > 0 {
                // put out the last few bits
                output_nbits(
                    outfile,
                    b3.bitbuffer & ((1 << b3.bits_to_go) - 1),
                    b3.bits_to_go as usize,
                    buffer2,
                );
            } else {
                // have to write a zero nybble if there are no 1's in array
                output_nbits(outfile, CODE[0], NCODE[0] as usize, buffer2);
            }
        } else {
            if b3.bits_to_go > 0 {
                // put out the last few bits

                output_nbits(
                    outfile,
                    b3.bitbuffer & ((1 << b3.bits_to_go) - 1),
                    b3.bits_to_go as usize,
                    buffer2,
                );
            }

            // write in blocks of 24 bits to speed things up
            for i in (0..b).rev() {
                output_nbits(outfile, i32::from(buffer[i]), 8, buffer2);
            }
        }
    }
    Ok(())
}

/// Encode values in quadrant of 2-D array using binary quadtree coding for each bit plane.  
///
/// Assumes array is positive.
///
/// # Arguments
///
/// * `outfile` - Output data stream
/// * `a` - Array of values to encode
/// * `n` - Number of elements in array
/// * `nqx` - Number of elements in X axis
/// * `nqy` - Number of elements in Y axis
/// * `nbitplanes` - Number of bit planes to encode
/// * `buffer2` - Buffer to write to
///
fn qtree_encode64(
    outfile: &mut Vec<u8>,
    a: &[i64],
    n: usize,
    nqx: usize,
    nqy: usize,
    nbitplanes: i32,
    buffer2: &mut Buffer2,
) -> Result<(), EncodeError> {
    let mut b: usize;

    let mut nx: usize;
    let mut ny: usize;

    let mut b3: Buffer3 = Buffer3 {
        bitbuffer: 0,
        bits_to_go: 0,
    };

    // log2n is log2 of max(nqx,nqy) rounded up to next power of 2
    let nqmax: usize = if nqx > nqy { nqx } else { nqy };
    let mut log2n: i32 = ((nqmax as f32).ln() / 2.0_f32.ln() + 0.5) as i32;

    if nqmax > (1 << log2n) {
        log2n += 1;
    }

    // initialize buffer point, max buffer size
    let nqx2: usize = nqx.div_ceil(2);
    let nqy2: usize = nqy.div_ceil(2);
    let bmax: usize = (nqx2 * nqy2 + 2) / 2;

    /*
    We're indexing A as a 2-D array with dimensions (nqx,nqy).
    Scratch is 2-D with dimensions (nqx/2,nqy/2) rounded up.
    Scr1 is used to store first level of quadtree in case direct
    coding is needed.
    Buffer is used to store string of codes for output.
     */
    let mut scratch: Vec<u8> = Vec::new();
    if scratch.try_reserve_exact(2 * bmax).is_err() {
        ffpmsg("qtree_encode64: insufficient memory");
        return Err(EncodeError::DataCompressionError);
    } else {
        scratch.resize(2 * bmax, 0);
    }

    let mut buffer: Vec<u8> = Vec::new();
    if buffer.try_reserve_exact(bmax).is_err() {
        ffpmsg("qtree_encode64: insufficient memory");
        return Err(EncodeError::DataCompressionError);
    } else {
        buffer.resize(bmax, 0);
    }

    // now encode each bit plane, starting with the top
    for bit in (0..(nbitplanes as usize)).rev() {
        // initialize bit buffer
        b = 0;
        b3.bitbuffer = 0;
        b3.bits_to_go = 0;

        // on first pass copy A to scr1 array
        qtree_onebit64(a, n, nqx, nqy, &mut scratch, bit);
        nx = (nqx + 1) >> 1;
        ny = (nqy + 1) >> 1;

        // copy non-zero values to output buffer, which will be written in reverse order
        if bufcopy(&scratch, nx * ny, &mut buffer, &mut b, bmax, &mut b3) {
            //quadtree is expanding data, change warning code and just fill buffer with bit-map
            write_bdirect64(outfile, a, n, nqx, nqy, &mut scratch, bit, buffer2);
            continue;
        }

        // do log2n reductions
        let mut is_continue = false;
        for _ in 1..log2n {
            qtree_reduce(&mut scratch, ny, nx, ny);
            nx = (nx + 1) >> 1;
            ny = (ny + 1) >> 1;
            if bufcopy(&scratch, nx * ny, &mut buffer, &mut b, bmax, &mut b3) {
                // nx=1, ny=1, b=1??, scratch=[10,8,0,0], bmax=2, b3=0/1
                write_bdirect64(outfile, a, n, nqx, nqy, &mut scratch, bit, buffer2);
                is_continue = true;
                break; // Break this loop and continue with next
            }
        }

        if is_continue {
            continue;
        }

        // OK, we've got the code in buffer
        // Write quadtree warning code, then write buffer in reverse order
        output_nybble(outfile, 0xF, buffer2);

        if b == 0 {
            if b3.bits_to_go > 0 {
                // put out the last few bits
                output_nbits(
                    outfile,
                    b3.bitbuffer & ((1 << b3.bits_to_go) - 1),
                    b3.bits_to_go as usize,
                    buffer2,
                );
            } else {
                // have to write a zero nybble if there are no 1's in array
                output_nbits(outfile, CODE[0], NCODE[0] as usize, buffer2);
            }
        } else {
            if b3.bits_to_go > 0 {
                // put out the last few bits

                output_nbits(
                    outfile,
                    b3.bitbuffer & ((1 << b3.bits_to_go) - 1),
                    b3.bits_to_go as usize,
                    buffer2,
                );
            }

            // write in blocks of 24 bits to speed things up
            for i in (0..b).rev() {
                output_nbits(outfile, i32::from(buffer[i]), 8, buffer2);
            }
        }
    }
    Ok(())
}

/* ######################################################################### */
/// copy non-zero codes from array to buffer
///
/// # Arguments
///
/// * `a` - Array of values to encode
/// * `n` - Number of elements in array
/// * `buffer` - Buffer to write to
/// * `b` - Index of buffer
/// * `bmax` - Maximum size of buffer
/// * `qtb` - Buffer3
///
/// # Returns
///
/// false if successful, true if buffer is full
fn bufcopy(
    a: &[u8],
    n: usize,
    buffer: &mut [u8],
    b: &mut usize,
    bmax: usize,
    qtb: &mut Buffer3,
) -> bool {
    for i in 0..n {
        if a[i] != 0 {
            // add Huffman code for a[i] to buffer
            qtb.bitbuffer |= CODE[a[i] as usize] << qtb.bits_to_go;
            qtb.bits_to_go += NCODE[a[i] as usize];

            if qtb.bits_to_go >= 8 {
                buffer[*b] = qtb.bitbuffer as u8;
                *b += 1;

                // return warning code if we fill buffer
                if *b >= bmax {
                    return true;
                }

                qtb.bitbuffer >>= 8;
                qtb.bits_to_go -= 8;
            }
        }
    }
    false
}

/// Do first quadtree reduction step on bit BIT of array A.
///
/// Results put into B.
///
/// # Arguments
///
/// * `a` - Array of values to encode
/// * `n` - Number of elements in array
/// * `nx` - Number of elements in X axis
/// * `ny` - Number of elements in Y axis
/// * `b` - Buffer to write to
/// * `bit` - Bit to encode
///
fn qtree_onebit(a: &[i32], n: usize, nx: usize, ny: usize, b: &mut [u8], bit: usize) {
    let mut s10: usize;
    let mut s00: usize;

    // use selected bit to get amount to shift
    let b0: i32 = 1 << bit;
    let b1: i32 = b0 << 1;
    let b2: i32 = b0 << 2;
    let b3: i32 = b0 << 3;
    let mut k: usize = 0;
    let mut ii = 0; /* k is index of b[i/2,j/2]	*/

    if nx == 0 {
        return;
    } // Stops underflow and early return

    for i in (0..(nx - 1)).step_by(2) {
        s00 = n * i; /* s00 is index of a[i,j]	*/
        s10 = s00 + n; /* s10 is index of a[i+1,j]	*/
        let mut ji = 0;

        if ny == 0 {
            continue;
        } // Stops underflow

        for _j in (0..(ny - 1)).step_by(2) {
            b[k] = (((a[s10 + 1] & b0)
                | ((a[s10] << 1) & b1)
                | ((a[s00 + 1] << 2) & b2)
                | ((a[s00] << 3) & b3))
                >> bit) as u8;

            k += 1;
            s00 += 2;
            s10 += 2;
            ji += 2;
        }
        if ji < ny {
            /*
             * row size is odd, do last element in row
             * s00+1,s10+1 are off edge
             */
            b[k] = ((((a[s10] << 1) & b1) | ((a[s00] << 3) & b3)) >> bit) as u8;
            k += 1;
        }
        ii += 2;
    }

    if ii < nx {
        /*
         * column size is odd, do last row
         * s10,s10+1 are off edge
         */
        s00 = n * ii;
        let mut ji = 0;

        if ny == 0 {
            return;
        } // Early return, prevent underflow on next line (ny - 1)

        for _j in (0..(ny - 1)).step_by(2) {
            b[k] = ((((a[s00 + 1] << 2) & b2) | ((a[s00] << 3) & b3)) >> bit) as u8;
            k += 1;
            s00 += 2;
            ji += 2;
        }

        if ji < ny {
            /*
             * both row and column size are odd, do corner element
             * s00+1, s10, s10+1 are off edge
             */
            b[k] = (((a[s00] << 3) & b3) >> bit) as u8;
        }
    }
}

/// Do first quadtree reduction step on bit BIT of array A.
///
/// Results put into B.
///
/// # Arguments
///
/// * `a` - Array of values to encode
/// * `n` - Number of elements in array
/// * `nx` - Number of elements in X axis
/// * `ny` - Number of elements in Y axis
/// * `b` - Buffer to write to
/// * `bit` - Bit to encode
///
fn qtree_onebit64(a: &[i64], n: usize, nx: usize, ny: usize, b: &mut [u8], bit: usize) {
    let mut s10: usize;
    let mut s00: usize;

    // use selected bit to get amount to shift
    let b0: i64 = 1 << bit;
    let b1: i64 = b0 << 1;
    let b2: i64 = b0 << 2;
    let b3: i64 = b0 << 3;
    let mut k: usize = 0;
    let mut ii = 0; /* k is index of b[i/2,j/2]	*/

    if nx == 0 {
        return;
    } // Stops underflow and early return

    for i in (0..(nx - 1)).step_by(2) {
        s00 = n * i; /* s00 is index of a[i,j]	*/
        s10 = s00 + n; /* s10 is index of a[i+1,j]	*/
        let mut ji = 0;

        if ny == 0 {
            continue;
        } // Stops underflow

        for _j in (0..(ny - 1)).step_by(2) {
            b[k] = (((a[s10 + 1] & b0)
                | ((a[s10] << 1) & b1)
                | ((a[s00 + 1] << 2) & b2)
                | ((a[s00] << 3) & b3))
                >> bit) as u8;

            k += 1;
            s00 += 2;
            s10 += 2;
            ji += 2;
        }
        if ji < ny {
            /*
             * row size is odd, do last element in row
             * s00+1,s10+1 are off edge
             */
            b[k] = ((((a[s10] << 1) & b1) | ((a[s00] << 3) & b3)) >> bit) as u8;
            k += 1;
        }
        ii += 2;
    }

    if ii < nx {
        /*
         * column size is odd, do last row
         * s10,s10+1 are off edge
         */
        s00 = n * ii;
        let mut ji = 0;

        if ny == 0 {
            return;
        } // Early return, prevent underflow on next line (ny - 1)

        for _j in (0..(ny - 1)).step_by(2) {
            b[k] = ((((a[s00 + 1] << 2) & b2) | ((a[s00] << 3) & b3)) >> bit) as u8;
            k += 1;
            s00 += 2;
            ji += 2;
        }

        if ji < ny {
            /*
             * both row and column size are odd, do corner element
             * s00+1, s10, s10+1 are off edge
             */
            b[k] = (((a[s00] << 3) & b3) >> bit) as u8;
        }
    }
}

/* ######################################################################### */
/// do one quadtree reduction step on array a
///
/// results put into b (which may be the same as a)
///
/// # Arguments
///
/// * `a` - Array of values to encode
/// * `n` - Number of elements in array
/// * `nx` - Number of elements in X axis
/// * `ny` - Number of elements in Y axis
///
fn qtree_reduce(a: &mut [u8], n: usize, nx: usize, ny: usize) {
    let mut s10;
    let mut s00;
    let mut k = 0;
    let mut ii = 0; /* k is index of b[i/2,j/2]	*/

    if nx == 0 {
        return;
    } // Early return and prevents underflows

    for i in (0..(nx - 1)).step_by(2) {
        s00 = n * i; /* s00 is index of a[i,j]	*/
        s10 = s00 + n; /* s10 is index of a[i+1,j]	*/

        if ny == 0 {
            ii += 2;
            continue;
        } // Early continue and prevents underflows

        for _j in (0..(ny - 1)).step_by(2) {
            a[k] = u8::from(a[s10 + 1] != 0)
                | (u8::from(a[s10] != 0) << 1)
                | (u8::from(a[s00 + 1] != 0) << 2)
                | (u8::from(a[s00] != 0) << 3);
            k += 1;
            s00 += 2;
            s10 += 2;
        }

        if ny % 2 != 0 {
            /*
             * row size is odd, do last element in row
             * s00+1,s10+1 are off edge
             */
            a[k] = (u8::from(a[s10] != 0) << 1) | (u8::from(a[s00] != 0) << 3);
            k += 1;
        }
        ii += 2;
    }

    if nx % 2 != 0 {
        /*
         * column size is odd, do last row
         * s10,s10+1 are off edge
         */
        s00 = n * ii;
        if ny == 0 {
            return;
        } // Early return and prevents underflows

        for _j in (0..(ny - 1)).step_by(2) {
            a[k] = (u8::from(a[s00 + 1] != 0) << 2) | (u8::from(a[s00] != 0) << 3);
            k += 1;
            s00 += 2;
        }

        if ny % 2 != 0 {
            /*
             * both row and column size are odd, do corner element
             * s00+1, s10, s10+1 are off edge
             */
            a[k] = u8::from(a[s00] != 0) << 3;
        }
    }
}

/* ######################################################################### */
#[allow(clippy::too_many_arguments)]
fn write_bdirect(
    outfile: &mut Vec<u8>,
    a: &[i32],
    n: usize,
    nqx: usize,
    nqy: usize,
    scratch: &mut [u8],
    bit: usize,
    buffer2: &mut Buffer2,
) {
    // Write the direct bitmap warning code
    output_nybble(outfile, 0, buffer2);

    // Copy A to scratch array (again!), packing 4 bits/nybble
    qtree_onebit(a, n, nqx, nqy, scratch, bit);

    // write to outfile
    output_nnybble(outfile, nqx.div_ceil(2) * nqy.div_ceil(2), scratch, buffer2);
}

/* ######################################################################### */
#[allow(clippy::too_many_arguments)]
fn write_bdirect64(
    outfile: &mut Vec<u8>,
    a: &[i64],
    n: usize,
    nqx: usize,
    nqy: usize,
    scratch: &mut [u8],
    bit: usize,
    buffer2: &mut Buffer2,
) {
    // Write the direct bitmap warning code
    output_nybble(outfile, 0, buffer2);

    // Copy A to scratch array (again!), packing 4 bits/nybble
    qtree_onebit64(a, n, nqx, nqy, scratch, bit);

    // write to outfile
    output_nnybble(outfile, nqx.div_ceil(2) * nqy.div_ceil(2), scratch, buffer2);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_output_nbits() {
        let mut charb: Vec<u8> = Vec::with_capacity(4);

        let mut b2 = Buffer2 {
            buffer2: 0,
            bits_to_go2: 8,
            bitcount: 0,
        };

        output_nbits(&mut charb, 23, 2, &mut b2);
        output_nbits(&mut charb, 23, 2, &mut b2);
        output_nbits(&mut charb, 23, 4, &mut b2);

        assert_eq!(b2.bitcount, 8);
        assert_eq!(b2.bits_to_go2, 8);
        assert_eq!(b2.buffer2, 247);
        assert_eq!(charb, [247].to_vec());
    }

    #[test]
    fn test_fits_compress() {
        let mut input: [i32; 16] = [2, 2, 1, 2, 3, 2, 7, 7, 4, 2, 2, 1, 2, 4, 25, 2];
        let mut output: Vec<u8> = Vec::with_capacity(16);

        let encoder = HCEncoder::new();
        let _res = encoder.write(&mut input, 4, 4, 0, &mut output);

        assert_eq!(output.len(), 48);
        assert_eq!(
            input,
            [0, 16, 2, 2, 12, 6, 0, 24, 2, 12, 1, 1, 0, 24, 4, 22]
        );
        assert_eq!(
            output,
            [
                221, 153, 0, 0, 0, 4, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 32, 5, 5, 5,
                245, 231, 227, 199, 253, 227, 199, 253, 247, 255, 120, 249, 245, 239, 254, 241,
                255, 124, 120, 251, 0, 68, 200
            ]
            .to_vec()
        );
    }

    #[test]
    fn test_fits_compress_strange_input() {
        let mut input: [i32; 10] = [-1, -1, -112, -1, 9983, -28528, -112, -1, -1, -1];
        let mut output: Vec<u8> = Vec::with_capacity(200);

        let encoder = HCEncoder::new();
        let _res = encoder.write(&mut input, 1, 10, 0, &mut output);
        println!("{:#?}", output);
        assert_eq!(output.len(), 101);

        assert_eq!(
            output,
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
    fn test_fits_compress_strange_input2() {
        let mut input: [i32; 12] = [
            -1, -1, -9584, -28561, -112, -24321, -1, -1, -1, -9584, -28561, -112,
        ];
        let mut output: Vec<u8> = Vec::with_capacity(200);

        let encoder = HCEncoder::new();
        let _res = encoder.write(&mut input, 1, 12, 0, &mut output);

        assert_eq!(output.len(), 106);

        assert_eq!(
            output,
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
    fn test_fits_compress_strange_input3() {
        let mut input: [i32; 16] = [
            2570, 28560, -5778, -28528, 28816, 28816, 2671, -246, -1, -28417, -5778, -28528, 28304,
            -28439, -28528, -5791,
        ];
        let mut output: Vec<u8> = Vec::with_capacity(200);

        let encoder = HCEncoder::new();
        let _res = encoder.write(&mut input, 1, 16, 0, &mut output);

        assert_eq!(output.len(), 140);

        assert_eq!(
            output,
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
    fn test_fits_compress_strange_input4() {
        let mut input: [i32; 82] = [
            -1, -1, 1, -256, -1, 0, 0, 0, 0, 0, -256, 1, -256, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            256, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, -1, -1, -1, -1,
            -1, -1, -1, -1, -1, -1, -1, -1, 255, 0, 0, 0, 0, -256, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, -256, -4097, -1,
        ];
        let mut output: Vec<u8> = Vec::with_capacity(200);

        let encoder = HCEncoder::new();
        let _res = encoder.write(&mut input, 1, 82, 0, &mut output);

        assert_eq!(output.len(), 154);

        assert_eq!(
            output,
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
    fn test_fits_compress_strange_input5() {
        let mut input: [i32; 10] = [61, 14, 0, 0, 0, -23641, -13558, -28528, -28526, -28528];
        let mut output: Vec<u8> = Vec::with_capacity(200);

        let encoder = HCEncoder::new();
        let _res = encoder.write(&mut input, 10, 1, 0, &mut output);

        assert_eq!(output.len(), 104);

        assert_eq!(
            output,
            [
                221, 153, 0, 0, 0, 1, 0, 0, 0, 10, 0, 0, 0, 0, 255, 255, 255, 255, 255, 247, 8, 32,
                19, 16, 0, 246, 215, 180, 255, 224, 64, 143, 109, 123, 104, 16, 35, 254, 246, 215,
                187, 107, 254, 246, 151, 253, 255, 0, 100, 0, 100, 127, 223, 247, 253, 237, 63,
                105, 123, 87, 237, 95, 180, 253, 165, 255, 127, 222, 213, 251, 180, 253, 171, 247,
                118, 253, 182, 16, 145, 238, 219, 255, 191, 239, 251, 254, 255, 191, 239, 251, 254,
                255, 191, 239, 251, 254, 255, 191, 239, 251, 254, 255, 128, 255
            ]
        );
    }

    #[test]
    fn test_fits_compress_strange_input6() {
        let mut input: [i32; 10] = [
            -28662, -28528, 18761, 18761, 18761, 18761, 18761, 18761, -28528, -28528,
        ];
        let mut output: Vec<u8> = Vec::with_capacity(200);

        let encoder = HCEncoder::new();
        let _res = encoder.write(&mut input, 10, 1, 0, &mut output);

        assert_eq!(output.len(), 84);

        assert_eq!(
            output,
            [
                221, 153, 0, 0, 0, 1, 0, 0, 0, 10, 0, 0, 0, 0, 255, 255, 255, 255, 255, 250, 185,
                160, 20, 9, 0, 246, 215, 253, 237, 95, 253, 238, 237, 123, 87, 238, 237, 127, 223,
                247, 182, 189, 221, 175, 119, 107, 221, 218, 246, 175, 221, 218, 246, 175, 254,
                255, 191, 239, 251, 219, 127, 247, 253, 255, 127, 222, 219, 246, 223, 253, 255,
                127, 223, 247, 253, 255, 127, 223, 247, 253, 255, 127, 192, 128
            ]
        );
    }

    // ORIGINAL LIBRARY fails to encode and then decode this input
    // The compression will work, but there is UB
    // This was found via fuzzing
    #[allow(dead_code)]
    fn test_fits_compress_strange_input7() {
        let mut input: [i32; 10] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 134217727];
        let mut output: Vec<u8> = Vec::with_capacity(200);

        let encoder = HCEncoder::new();
        let _res = encoder.write(&mut input, 1, 10, 0, &mut output);

        assert_eq!(output.len(), 146);

        assert_eq!(
            output,
            [
                221, 153, 0, 0, 0, 10, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 64, 0, 0, 0, 31, 28, 0,
                246, 255, 239, 251, 254, 255, 191, 239, 251, 254, 255, 191, 239, 251, 254, 255,
                191, 239, 251, 254, 255, 191, 239, 251, 254, 255, 191, 239, 251, 254, 255, 191,
                239, 251, 254, 255, 191, 239, 251, 254, 255, 191, 239, 251, 254, 255, 191, 239,
                251, 254, 255, 191, 239, 251, 254, 255, 191, 239, 251, 254, 255, 191, 239, 251,
                254, 255, 191, 239, 251, 254, 255, 191, 239, 251, 203, 126, 91, 242, 223, 150, 252,
                183, 229, 191, 45, 249, 111, 203, 126, 91, 242, 223, 150, 252, 183, 229, 191, 45,
                249, 111, 203, 126, 91, 242, 223, 150, 252, 183, 229, 191, 45, 249, 111, 203, 126,
                91, 242, 223, 252, 0, 0
            ]
        );
    }

    #[test]
    fn test_64bit_input1() {
        let mut input: [i64; 2] = [0, 1];
        let mut output: Vec<u8> = Vec::with_capacity(200);

        let encoder = HCEncoder::new();
        let _res = encoder.write64(&mut input, 1, 2, 0, &mut output);

        assert_eq!(output.len(), 32);

        assert_eq!(
            output,
            [
                221, 153, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 2, 0, 255,
                191, 239, 127, 240, 0, 0
            ]
        );
    }

    #[test]
    fn test_htrans() {
        let mut input: [i32; 16] = [2, 2, 1, 2, 3, 2, 7, 7, 4, 2, 2, 1, 2, 4, 25, 2];
        htrans(&mut input, 4, 4).unwrap();

        assert_eq!(
            input,
            [32, 16, -2, 2, 12, 6, 0, -24, 2, 12, -1, -1, 0, 24, 4, -22]
        );
    }

    #[test]
    fn test_digitize_0_scale() {
        let mut input: [i32; 16] = [32, 16, -2, 2, 12, 6, 0, -24, 2, 12, -1, -1, 0, 24, 4, -22];
        digitize(&mut input, 4, 4, 0);
        assert_eq!(
            input,
            [32, 16, -2, 2, 12, 6, 0, -24, 2, 12, -1, -1, 0, 24, 4, -22]
        );
    }

    #[test]
    fn test_qtree_onebit() {
        let nx = 2;
        let ny = 2;
        let bit = 4;
        let n = 4;
        let a: [i32; 16] = [32, 16, -2, 2, 12, 6, 0, -24, 2, 12, -1, -1, 0, 24, 4, -22];
        let mut scratch: [u8; 2] = [0; 2];

        qtree_onebit(&a, n, nx, ny, &mut scratch, bit);

        assert_eq!(
            a,
            [32, 16, -2, 2, 12, 6, 0, -24, 2, 12, -1, -1, 0, 24, 4, -22]
        );
        assert_eq!(scratch, [4, 0]);
    }

    #[test]
    fn test_qtree_onebit_edge() {
        // ny = 0 causes underflow due to usize
        let nx = 5;
        let ny = 0;
        let bit = 16;
        let n = 1;
        let a: [i32; 9] = [37536, 37088, 224, 36864, 0, 222, 77022, 222, 0];
        let mut scratch: [u8; 4] = [99; 4];

        qtree_onebit(&a, n, nx, ny, &mut scratch, bit);

        assert_eq!(a, [37536, 37088, 224, 36864, 0, 222, 77022, 222, 0]);
        assert_eq!(scratch, [99, 99, 99, 99]);
    }

    #[test]
    fn test_bufcopy() {
        let a: [u8; 2] = [4, 0];
        let mut buffer: [u8; 1] = [0];
        let nx = 1;
        let ny = 1;
        let bmax = 1;
        let mut b = 0;
        let mut b3: Buffer3 = Buffer3 {
            bitbuffer: 0,
            bits_to_go: 0,
        };

        let _res = bufcopy(&a, nx * ny, &mut buffer, &mut b, bmax, &mut b3);

        assert_eq!(b3.bitbuffer, 2);
        assert_eq!(b3.bits_to_go, 3);
        assert_eq!(b, 0);
        assert_eq!(buffer[0], 0); // Given bits_to_go is not >= 8, don't expect anything

        // Go again but force buffer to fill
        b3.bits_to_go = 7;
        let _res = bufcopy(&a, nx * ny, &mut buffer, &mut b, bmax, &mut b3);

        assert_eq!(b3.bitbuffer, 258);
        assert_eq!(b3.bits_to_go, 10);
        assert_eq!(b, 1);
        assert_eq!(buffer[0], 2); // Given bits_to_go is not >= 8, don't expect anything
    }

    #[test]
    pub fn test_done_outputting_bits() {
        let mut b2 = Buffer2 {
            buffer2: -137916496,
            bits_to_go2: 4,
            bitcount: 164,
        };

        let mut outfile: Vec<u8> = Vec::new();
        done_outputing_bits(&mut outfile, &mut b2);

        assert_eq!(outfile, [0]);
        assert_eq!(b2.bits_to_go2, 4);
        assert_eq!(b2.bitcount, 168);
        assert_eq!(b2.buffer2, -137916496);
    }

    #[test]
    pub fn test_output_nybble() {
        let mut outfile: Vec<u8> = Vec::new();
        let mut b2 = Buffer2 {
            buffer2: -137916496,
            bits_to_go2: 4,
            bitcount: 164,
        };
        let bits = 4;

        output_nybble(&mut outfile, bits, &mut b2);

        assert_eq!(outfile, [4]);
        assert_eq!(b2.bits_to_go2, 8);
        assert_eq!(b2.bitcount, 168);
        assert_eq!({ b2.buffer2 }, 2088303364); // We use usize so need to truncate
    }

    #[test]
    pub fn test_output_nnybble() {
        let mut outfile: Vec<u8> = Vec::new();
        let mut b2 = Buffer2 {
            buffer2: 0,
            bits_to_go2: 4,
            bitcount: 4,
        };
        let n = 3;
        let array = [2, 8, 8, 85];

        output_nnybble(&mut outfile, n, &array, &mut b2);

        let _t: i8 = -120;

        assert_eq!(outfile, [2, _t as u8]); // Same as 2,136
        assert_eq!(b2.bits_to_go2, 8);
        assert_eq!(b2.bitcount, 8);
        assert_eq!({ b2.buffer2 }, 0); // We use usize so need to truncate
    }
}
