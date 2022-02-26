#include <assert.h>
#include "fits_hdecompress.c"


/****** DECOMPRESS *******/
void test_fits_hdecompress() {

    char input[48] = {
            -35,-103,0,0,0,4,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,32,5,5,5,-11,-25,-29,-57,-3,-29,-57,-3,-9,-1,120,-7,-11,-17,-2,-15,-1,124,120,-5,0,68,-56 };

    //int output[16];
    int *output;
    output = calloc(16, sizeof(int));

    int expected[16] = {2,2,1,2,3,2,7,7,4,2,2,1,2,4,25,2};
    int nx;
    int ny;
    int scale;
    int status = 0;
    int res = fits_hdecompress(input, 0, output, &ny, &nx, &scale, &status);

        //assert_eq!(output.len(), 48);
    assert( nx== 4 );
    assert(ny == 4);
    assert(scale == 0);
    assert(!memcmp( output, expected, 16*sizeof(int) ));


    for (int i=0; i<16; i++) {
        printf("%d,", output[i]);
    }

}

void test_fits_hdecompress_strange_input() {

    unsigned char input[101] ={221,153,0,0,0,10,0,0,0,1,0,0,0,0,255,255,255,255,255,255,109,64,16,17,0,2,136,255,191,224,40,143,251,254,246,207,253,238,168,251,53,238,168,255,223,247,253,255,127,223,247,253,255,127,223,247,253,255,127,223,247,253,255,127,223,247,253,255,127,223,247,178,255,239,251,217,127,247,178,253,151,255,127,222,234,143,186,163,255,123,170,62,234,143,186,163,238,168,255,192,100};
            
            //int output[16];
    int *output;
    output = calloc(200, sizeof(int));

    int expected[10] = {-1, -1, -112, -1, 9983, -28528, -112, -1, -1, -1};
    int nx;
    int ny;
    int scale;
    int status = 0;
    int res = fits_hdecompress(input, 0, output, &ny, &nx, &scale, &status);

        //assert_eq!(output.len(), 48);
    assert(nx== 10 );
    assert(ny == 1);
    assert(scale == 0);
    //assert(!memcmp( output, expected, 10*sizeof(int) ));


    for (int i=0; i<10; i++) {
        printf("%d,", output[i]);
    }
}

void test_fits_hdecompress_strange_input2() {

    unsigned char input[106] ={221,153,0,0,0,12,0,0,0,1,0,0,0,0,255,255,255,255,255,251,193,64,17,16,0,246,79,151,94,233,144,42,143,46,128,170,0,130,128,42,0,32,128,130,0,130,0,162,191,239,251,254,255,191,239,251,254,255,191,239,251,254,255,191,239,251,254,255,191,239,251,254,255,191,224,40,47,46,189,150,5,4,1,20,5,5,229,215,178,252,182,4,21,236,191,45,249,111,203,96,81,95,240,0,175,0};
       
            //int output[16];
    int *output;
    output = calloc(200, sizeof(int));

    int expected[12] = {-1,
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
        -112};
    int nx;
    int ny;
    int scale;
    int status = 0;
    int res = fits_hdecompress(input, 0, output, &ny, &nx, &scale, &status);

        //assert_eq!(output.len(), 48);
    assert(nx== 12 );
    assert(ny == 1);
    assert(scale == 0);
    //assert(!memcmp( output, expected, 10*sizeof(int) ));


    for (int i=0; i<12; i++) {
        printf("%d,", output[i]);
    }
}

int main_decompress() {
    //test_fits_hdecompress();
    test_fits_hdecompress_strange_input2();
}

int main() {
    main_decompress();
    return 0;
}

/*
int random() {

    //test_output_nbits();
    //test_qtree_onebit();
    //test_bufcopy();
    //test_done_outputing_bits();
    //test_output_nybble();


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
    unsigned char output[100];
    int status = 0;
    long n_bytes = 100;
    int res = fits_hcompress(&input, 4, 4, 0, &output, &n_bytes, &status);
	
    //int res = htrans(idata, 10, 2);
    //htrans_og(iidata, 10, 2);
    for (int i=0; i<48; i++) {
        printf("%d,", output[i]);
    }
    return 0;
}

*/