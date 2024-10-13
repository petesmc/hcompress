#include <assert.h>
#include "fits_hcompress.c"


void qtree_onebit_og(int a[], int n, int nx, int ny, unsigned char b[], int bit) {
    int i, *p, *pend;
    unsigned char *pb, *pb0;
    int mbit, bitm2;

	/*
	 * mask to get selected bit
	 */
	mbit = 1<<bit;
	pb = b;
	bitm2 = bit - 2;
	for (i = 0; i<nx; i += 2) {
		pb0 = pb;
		pend = &a[n*i+ny-1];
		switch (bit) {
			case 0:
				for (p = &a[n*i]; p < pend; p += 2, pb += 1)
					*pb = ((*p & mbit)<<3) | ((*(p+1) & mbit)<<2);
				break;
			case 1:
				for (p = &a[n*i]; p < pend; p += 2, pb += 1)
					*pb = ((*p & mbit)<<2) | ((*(p+1) & mbit)<<1);
				break;
			case 2:
				for (p = &a[n*i]; p < pend; p += 2, pb += 1)
					*pb = ((*p & mbit)<<1) | ((*(p+1) & mbit)   );
				break;
			default:
				for (p = &a[n*i]; p < pend; p += 2, pb += 1)
					*pb = ( ((*p & mbit)<<1) | (*(p+1) & mbit) ) >> bitm2;
		}
		if (p == pend) {
			/*
			 * row size is odd, do last element in row
			 * *(p+1) is off edge
			 */
			*pb = ((*p & mbit)<<3) >> bit;
			pb += 1;
		}
		if (i < nx-1) {
			/*
			 * not on last row, add in next row
			 */
			pb = pb0;
			pend = &a[n*(i+1)+ny-1];
			switch (bit) {
				case 0:
					for (p = &a[n*(i+1)]; p < pend; p += 2, pb += 1)
						*pb = ( ((*p & mbit)<<1) |  (*(p+1) & mbit)     ) | *pb;
					break;
				case 1:
					for (p = &a[n*(i+1)]; p < pend; p += 2, pb += 1)
						*pb = (  (*p & mbit)     | ((*(p+1) & mbit)>>1) ) | *pb;
					break;
				default:
					for (p = &a[n*(i+1)]; p < pend; p += 2, pb += 1)
						*pb = ( ( ((*p & mbit)<<1) | (*(p+1) & mbit) ) >> bit) | *pb ;
					break;
			}
			if (p == pend) {
				/* odd row size */
				*pb = ( ((*p & mbit)<<1) >> bit) | *pb ;
				pb += 1;
			}
		}
	}
}

/***** COMPRESSS *********/
void test_output_nbits() {
    unsigned char charb[20];

    noutmax = 20;
	buffer2 = 0;			/* Buffer is empty to start	*/
	bits_to_go2 = 8;		/* with				*/
	bitcount = 0;

    output_nbits(charb ,23, 2);
	output_nbits(charb ,23, 2);
	output_nbits(charb ,23, 4);

    assert(bitcount == 8);
    assert(bits_to_go2 == 8 );
    assert(buffer2 == 247);
    assert(charb[0] ==247);
}

void test_qtree_onebit() {
    
    int a[16] = {32,16,-2,2,12,6,0,-24,2,12,-1,-1,0,24,4,-22};
    int n = 4;
    int nx = 2;
    int ny = 2;
    int bit = 4;

    unsigned char b[2];
    
    qtree_onebit(a, n, nx, ny, b, bit);

    assert(b[0]==4);
    assert(b[1]==0);

}

void test_bufcopy() {

    unsigned char a[2] = {4, 0};
    int nx = 1;
    int ny = 1;
    int bmax = 1;
    int b = 0;
    unsigned char buffer[1] = {0};

    //statics
    bitbuffer = 0;
	bits_to_go3 = 0;

    bufcopy(a,nx*ny,buffer,&b,bmax);

    assert(bitbuffer==2);
    assert(bits_to_go3==3);
    assert(buffer[0]==0);
    assert(b==0);

    bits_to_go3 = 7;
    bufcopy(a,nx*ny,buffer,&b,bmax);

    assert(bitbuffer==258);
    assert(bits_to_go3==10);
    assert(b==1);
    assert(buffer[0]==2);

}

void test_done_outputing_bits() {
    unsigned char outfile[1] = {99}; //Dummy value

    noutmax = 1;
	buffer2 = -137916496;			/* Buffer is empty to start	*/
	bits_to_go2 = 4;		/* with				*/
	bitcount = 164;

    done_outputing_bits(outfile);
	
    assert(outfile[0] == 0);
    assert(bits_to_go2 == 4 );
    assert(bitcount == 168 );
    assert(buffer2 == -137916496);
    
}

void test_output_nybble() {
    unsigned char outfile[1] = {99}; //Dummy value
    int bits = 4;
    noutmax = 1;
	buffer2 = -137916496;			/* Buffer is empty to start	*/
	bits_to_go2 = 4;		/* with				*/
	bitcount = 164;

    output_nybble(outfile, bits);
	
    assert(outfile[0] == 4);
    assert(bits_to_go2 == 8 );
    assert(bitcount == 168 );
    assert(buffer2 == 2088303364);
    
}



/** INTEGRATION ***/
void test_small_input() {
    int input[2] = { -1, 640090111 };
    unsigned char output[200];
    int status = 0;
    long n_bytes = 200;
    int res = fits_hcompress(&input, 1, 2, 0, &output, &n_bytes, &status);

    for (int i=0; i<200; i++) {
        printf("%d,", output[i]);
    }
}

void test_strange_input() {
    int input[10] = { -1,
                -1,
                -112,
                -1,
                9983,
                -28528,
                -112,
                -1,
                -1,
                -1 };
    unsigned char output[200];
    int status = 0;
    long n_bytes = 200;
    int res = fits_hcompress(&input, 1, 10, 0, &output, &n_bytes, &status);

    printf("Bytes: %d\n", n_bytes);

    for (int i=0; i<101; i++) {
        printf("%d,", output[i]);
    }
}

void test_strange_input2() {
    int input[12] = { -1,
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
        -112 };
    unsigned char output[200];
    int status = 0;
    long n_bytes = 200;
    int res = fits_hcompress(&input, 1, 12, 0, &output, &n_bytes, &status);

    printf("Bytes: %d\n", n_bytes);

    for (int i=0; i<106; i++) {
        printf("%d,", output[i]);
    }
}

void test_strange_input3() {
    int input[16] = { 2570,
        28560,
        -5778,
        -28528,
        28816,
        28816,
        2671,
        -246,
        -1,
        -28417,
        -5778,
        -28528,
        28304,
        -28439,
        -28528,
        -5791 };
    unsigned char output[200];
    int status = 0;
    long n_bytes = 200;
    int res = fits_hcompress(&input, 1, 16, 0, &output, &n_bytes, &status);

    printf("Bytes: %d\n", n_bytes);

    for (int i=0; i<140; i++) {
        printf("%d,", output[i]);
    }
}

void test_strange_input4() {
    int input[82] = {  -1,
                -1,
                1,
                -256,
                -1,
                0,
                0,
                0,
                0,
                0,
                -256,
                1,
                -256,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                256,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                -1,
                -1,
                -1,
                -1,
                -1,
                -1,
                -1,
                -1,
                -1,
                -1,
                -1,
                -1,
                -1,
                -1,
                255,
                0,
                0,
                0,
                0,
                -256,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                -256,
                -4097,
                -1 };
    unsigned char output[800];
    int status = 0;
    long n_bytes = 800;
    int res = fits_hcompress(&input, 1, 82, 0, &output, &n_bytes, &status);

    printf("Bytes: %d\n", n_bytes);

    for (int i=0; i<154; i++) {
        printf("%d,", output[i]);
    }
}

void test_strange_input5() {
    int input[10] = { 61,
        14,
        0,
        0,
        0,
        -23641,
        -13558,
        -28528,
        -28526,
        -28528 };
    unsigned char output[200];
    int status = 0;
    long n_bytes = 200;
    int res = fits_hcompress(&input, 10, 1, 0, &output, &n_bytes, &status);

    printf("Bytes: %d\n", n_bytes);

    for (int i=0; i<104; i++) {
        printf("%d,", output[i]);
    }
}

void test_strange_input6() {
    int input[10] = { -28662,
        -28528,
        18761,
        18761,
        18761,
        18761,
        18761,
        18761,
        -28528,
        -28528};
    unsigned char output[200];
    int status = 0;
    long n_bytes = 200;
    int res = fits_hcompress(&input, 10, 1, 0, &output, &n_bytes, &status);

    printf("Bytes: %d\n", n_bytes);

    for (int i=0; i<84; i++) {
        printf("%d,", output[i]);
    }
}

void test_strange_input7() {
    int input[10] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 13421772};
    unsigned char output[200];
    int status = 0;
    long n_bytes = 200;
    int res = fits_hcompress(&input, 1, 2, 0, &output, &n_bytes, &status);

    printf("Bytes: %ld\n", n_bytes);

    for (int i=0; i<n_bytes; i++) {
        printf("%d,", output[i]);
    }
}



void test_qtree_onebit_edge_new() {
/*
 * a = (9) &[37536, 37088, 224, 36864, 0, 222, 77022, 222, 0]
 * n = 1
 * nx = 5
 * ny = 0
 * b = [0,0,1,8]
 * bit = 16
 */

    int a[9] = {37536, 37088, 224, 36864, 0, 222, 77022, 222, 0};
    int n = 1;
    int nx = 5;
    int ny = 0;
    int bit = 16;

    unsigned char b[4] = {99,99,99,99};
    
    qtree_onebit(a, n, nx, ny, b, bit);

    assert(b[0]==99);
    assert(b[1]==99);
    assert(b[2]==99);
    assert(b[3]==99);
}

void test_qtree_onebit_edge() {


    int a[9] = {37536, 37088, 224, 36864, 0, 222, 77022, 222, 0};
    int n = 1;
    int nx = 5;
    int ny = 0;
    int bit = 16;

    unsigned char b[4] = {99,99,99,99};
    
    qtree_onebit_og(a, n, nx, ny, b, bit);

    assert(b[0]==99);
    assert(b[1]==99);
    assert(b[2]==99);
        assert(b[3]==99);

}


int main() {

    //test_output_nbits();
    //test_qtree_onebit();
    //test_bufcopy();
    //test_done_outputing_bits();
    //test_output_nybble();
    //test_small_input();
    //test_strange_input();
    //test_strange_input2();
    //test_strange_input3();
    //test_strange_input4();
    //test_strange_input5();
    //test_strange_input6();
    test_strange_input7();
   //test_qtree_onebit_edge_new();
/*
    int idata[20] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};
    int iidata[20] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};
    int *tmp;
    tmp = malloc(20*sizeof(int));
    int *tmp2;
    tmp2 = malloc(20*sizeof(int));
	unsigned char *tmp3;
	tmp3 = calloc(20, sizeof(unsigned char));
	
	unsigned char charb[20];
	unsigned char a = 42;
	unsigned char b = ( (a != 0) << 2);


    int input[16] = {2,2,1,2,3,2,7,7,4,2,2,1,2,4,25,2};
    //unsigned char *output;
    //output = calloc(100, sizeof(unsigned char));
    char output[100];
    int status = 0;
    long n_bytes = 100;
    int res = fits_hcompress(&input, 4, 4, 0, &output, &n_bytes, &status);
	
    //int res = htrans(idata, 10, 2);
    //htrans_og(iidata, 10, 2);
    for (int i=0; i<48; i++) {
        printf("%d,", output[i]);
    }

    */
    return 0;
}
