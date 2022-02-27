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

void test_fits_hdecompress_strange_input3() {

    unsigned char input[140] ={221,153,0,0,0,16,0,0,0,1,0,0,0,0,255,255,255,255,255,254,197,32,19,17,0,246,207,253,245,68,211,250,117,84,126,153,116,5,21,31,78,153,63,84,106,171,234,139,170,2,138,175,170,102,143,186,39,233,146,126,155,100,8,160,190,153,170,242,207,253,255,127,223,247,253,255,127,223,247,253,255,127,223,247,253,255,127,223,247,253,255,127,222,75,250,167,85,95,77,183,245,78,137,2,42,175,170,109,191,166,89,8,40,143,170,116,104,8,136,190,137,103,234,153,39,233,178,66,136,163,234,155,53,244,77,31,248,0,158,248};
       
            //int output[16];
    int *output;
    output = calloc(200, sizeof(int));

    int expected[16] = {2570,
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
        -5791};
    int nx;
    int ny;
    int scale;
    int status = 0;
    int res = fits_hdecompress(input, 0, output, &ny, &nx, &scale, &status);

        //assert_eq!(output.len(), 48);
    assert(nx== 16 );
    assert(ny == 1);
    assert(scale == 0);
    //assert(!memcmp( output, expected, 10*sizeof(int) ));


    for (int i=0; i<16; i++) {
        printf("%d,", output[i]);
    }
}

void test_fits_hdecompress_strange_input4() {

    unsigned char input[154] ={221, 153, 0, 0, 0, 82, 0, 0, 0, 1, 0, 0, 0, 0, 255, 255, 255, 255, 255, 253, 245,
                0, 18, 14, 0, 246, 219, 103, 254, 246, 217, 103, 219, 101, 159, 109, 150, 125, 182,
                89, 246, 221, 50, 79, 182, 233, 209, 100, 10, 32, 136, 130, 0, 160, 160, 128, 0,
                136, 35, 221, 55, 69, 146, 91, 247, 77, 209, 100, 150, 253, 211, 117, 76, 146, 91,
                126, 233, 186, 166, 104, 150, 232, 251, 167, 84, 203, 44, 151, 79, 217, 116, 232,
                183, 236, 179, 76, 255, 223, 247, 253, 255, 127, 223, 247, 253, 255, 127, 223, 247,
                253, 255, 127, 223, 247, 253, 255, 121, 101, 183, 255, 127, 223, 247, 211, 166,
                106, 137, 162, 219, 166, 89, 108, 159, 251, 254, 255, 191, 239, 251, 254, 255, 189,
                211, 102, 139, 166, 89, 255, 128, 141, 59, 83, 9, 0};
       
            //int output[16];
    int *output;
    output = calloc(200, sizeof(int));

    int expected[82] = {-1, -1, 1, -256, -1, 0, 0, 0, 0, 0, -256, 1, -256, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            256, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, -1, -1, -1, -1,
            -1, -1, -1, -1, -1, -1, -1, -1, 255, 0, 0, 0, 0, -256, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, -256, -4097, -1};
    int nx;
    int ny;
    int scale;
    int status = 0;
    int res = fits_hdecompress(input, 0, output, &ny, &nx, &scale, &status);

        //assert_eq!(output.len(), 48);
    assert(nx== 82 );
    assert(ny == 1);
    assert(scale == 0);
    assert(!memcmp( output, expected, 82*sizeof(int) ));


    for (int i=0; i<16; i++) {
        printf("%d,", output[i]);
    }
}

void test_input_nnybble() {
    int n = 4;
    unsigned char array[4] = {0,10,8,0};
    bits_to_go = 5; //global
    buffer2 = 2123985925;//global
    nextchar = 38; //global

    unsigned char infile[140] ={221,153,0,0,0,16,0,0,0,1,0,0,0,0,255,255,255,255,255,254,197,32,19,17,0,246,207,253,245,68,211,250,117,84,126,153,116,5,21,31,78,153,63,84,106,171,234,139,170,2,138,175,170,102,143,186,39,233,146,126,155,100,8,160,190,153,170,242,207,253,255,127,223,247,253,255,127,223,247,253,255,127,223,247,253,255,127,223,247,253,255,127,222,75,250,167,85,95,77,183,245,78,137,2,42,175,170,109,191,166,89,8,40,143,170,116,104,8,136,190,137,103,234,153,39,233,178,66,136,163,234,155,53,244,77,31,248,0,158,248};
    
    int res = input_nnybble(&infile, n, array);

    assert(res==8);
    assert(bits_to_go==5);
    assert(buffer2=1946490143);

    unsigned char expected[4] = {2, 8, 10, 8};
    assert(!memcmp( array, expected, 4*sizeof(unsigned char) ));

}

int main_decompress() {
    //test_fits_hdecompress();
    //test_fits_hdecompress_strange_input3();
    //test_input_nnybble();
    test_fits_hdecompress_strange_input4();
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