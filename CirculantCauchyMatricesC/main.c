//
//  main.c
//  Circulant Cauchy
//
//  Created by Christian Schindelhauer on 27.04.12.
//  Copyright (c) 2012 Rechnernetze und Telematik an der Albert-Ludwigs-Universität Freiburg. All rights reserved.
//

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "artin-constants.h"

long *newBlock(int wordsize); 

void randomBlock(long *block, int wordsize);
void zeroBlock(long *block, int wordsize);
void oneBlock(long *block, int wordsize);

void printBlock(long *block, int wordsize);
void printFirstBit(long *block, int wordsize);
void printlnFirstBit(long *block, int wordsize);

void rotateBlock(long *block, int wordsize, int offset);
void copy_and_rotate_Block(long *a, long *b, int wordsize, int offset);
void copyBlock(long *a, long *b, int wordsize);
void addBlock(long *a, long *b, int wordsize);

void testBlock(void);

void copy_and_multiply_with_sum_of_power_of_twos(long *a, long *b, int artinprime, int p1, int p2);
void multiply_with_sum_of_power_of_twos(long *a, int artinprime, int p1, int p2);
void copy_and_divide_by_sum_of_power_of_twos(long *a, long *b, int artinprime, int p1, int p2);  
void divide_by_sum_of_power_of_twos(long *a, int artinprime, int p1, int p2);  

void extendToArtinPrime(long *a, long *b, int wordsize); 
void reduceToWordSize(long *a, long *b, int wordsize);   

int X(int i,int artinprime);
int Y(int i,int artinprime);

long *newVector(int wordsize, int n);

void randomVector(long *data, int wordsize, int n);
void addVector(long *data, long *operand, int wordsize, int n);
void vectorPrint(long *data, int wordsize, int n);

void testVector(void);

void innerCauchy(long *code, long *data, int artin, int n, int m); 
void innerCauchyXY(long *data, long *code, int wordsize, int n, int m, int *Xs, int *Ys);

void testCauchy(void);

void divide_by_A(long *result, int k, int artinprime, int n, int *Xs);
void divide_by_B(long *result, int k, int wordsize, int n, int *Ys);
void multiply_by_E(long *result, int k, int wordsize, int n, int *Xs, int *Ys);
void multiply_by_F(long *result, int k, int wordsize, int n, int *Xs, int *Ys);

void innerInverseCauchy(long *data, long *code, int wordsize, int n, int *Xs, int *Ys);
// int data[n][artin]
// int code[n][artin]
// int Xs[n] = powers of two (NOT!! the index)
// int Ys[n] = powers of two (NOT!! the index)
void innerReconstructData(long *data, long *code, int artinprime, int n, int m, int *Xindex, int *Ylist);
    // reconstructs the data
    // data consists partially of the available data 
    // boolean array Xindex[n] indicates wether data is available X[i]=1 <=> data is there
    // Ylist is the list of the Y-indices
    // n is the full data length n
    // m is the number of given encodings

void extendToArtinPrime(long *a, long *b, int wordsize); 
void reduceToWordSize(long *a, long *b, int wordsize);   
void extendToArtinPrimeVector(long *cir, long *data, int wordsize, int n);
void reduceToWordSizeVector(long *code, long *cir, int wordsize, int n);

int isArtin(int i); 
int ArtinSplitter(int *w,int wordsize, int maxSplitters); // w is int[maxSplitters]
void testArtin(void);

void CauchyMDSArtin(long *code, long *data, int wordsize, int n, int m); // wordsize is the real wordsize
void reconstructCauchyMDSArtin(long *data, long *code, 
                                       int wordsize, 
                                       int n, int m, 
                                       int *Xindex, int *Ylist);
void testMDSARtin(void);

void CauchyMDS(long *code, long *data, int wordsize, int n, int m); 
void reconstructCauchyMDS(long *data, long *code, 
                          int wordsize, 
                          int n, int m, 
                          int *Xindex, int *Ylist);
void testMDS(void);

void createDataFileNames (char * datafilename,  int n, const char * origFileName);
void analyzeFilename(char filename[], char masterfile[], int * fitsScheme, 
                     int * isData, long * index, 
                     int * n, int * blocksize, int * word);
// IMPLEMENTATION


long *newBlock(int wordsize) {
    // allocates a blog of wordsize integers
    // this allows to use the word parallelism of the processor
    // the size of the block is not stored
    long *block;
    block = malloc(wordsize * sizeof(long));
    return block;
}

void printBlock(long *block, int wordsize) {
    // for debugging: print all bits of integers of a block
    int i,j;
    long z;
    for (i=0; i<wordsize; i++) {
        z = 1;
        for (j=0; j<sizeof(long)*8; j++) {
            if ((block[i] & z) != 0) {
                printf("1");
            } else {
                printf("0");
            }
            z <<= 1;
        }
        printf("\n");
    }
    printf("\n");    
}

void printFirstBit(long *block, int wordsize) {
    // This prints only the first bit of each block
    // for debugging of the number theoretic aspect of Circulent Boolean matrices
    for (int i=0; i<wordsize; i++) {
        if ((block[i] & 1) != 0) {
            printf("1");
        } else {
            printf("0");
        }
    }
}
void printlnFirstBit(long *block, int wordsize) {
    // same as printFirstBit (same as printfirstBit but with new line
    printFirstBit(block,  wordsize); 
    printf("\n");
}

void randomBlock(long *block, int wordsize) {
    // generates a random block for debugging
    // note that it is not initialized and delivers always the same series
    for (int i=0; i<wordsize; i++) {
        block[i]=(long) rand()+ ((long) rand()<<33);
    }
}

void zeroBlock(long *block, int wordsize) {
    // sets block to 0
    for (int i=0; i<wordsize; i++) {
        block[i]=0;
    }
}

void oneBlock(long *block, int wordsize) {
    // sets block to Cir(1)
    block[0]= -1;
    for (int i=1; i<wordsize; i++) {
        block[i]=0;
    }    
}


void copyBlock(long *a, long *b, int wordsize) {
    // Block a = Block b
    memcpy(a,b,wordsize*sizeof(long));
}

void addBlock(long *a, long *b, int wordsize) {
    // Cir(a) = Cir(a) + Cir(b)
    // addition is Xor (^)
    for (int i=0; i<wordsize; i++) {
        a[i] ^= b[i];
    }
}

void copy_and_rotate_Block(long *a, long *b, int wordsize, int offset) {
    // required: 0 <= offset  < wordsize
    // wordsize ist prim mit offset
    for (int i=0;i<wordsize; i++) {
        a[(i+offset) % wordsize] = b[i];
    }
}

void rotateBlock(long *block, int artinprime, int offset) {
    // rotates block by an offset 
    // this is the multiplication with the circulant matrix Cir(2^offset)
    // IMPORTANT REQUISITE: offset >= 0 
    // otherwise the horrendous stupid definition of the C Modulo returns an error
    // 0 <= offset 
    // wordsize is prime w.r.t. offset
    
    long *h;
    h = newBlock(artinprime);
    copy_and_rotate_Block(h, block, artinprime, offset);
    copyBlock(block,h,artinprime);
    free(h);
}

void copy_and_multiply_with_sum_of_power_of_twos(long *a, long *b, int wordsize, int p1, int p2) {
    // Cir(a) = Cir(b) * Cir(2^p1+2^p2)
    copy_and_rotate_Block(a, b, wordsize, p1);
    for (int i=0;i<wordsize; i++) {
        a[(i+p2) % wordsize] ^= b[i];
    }    
}

void copy_and_divide_by_sum_of_power_of_twos(long *b, long *a, int wordsize, int p1, int p2) {
    // Cir(b) = Cir(a) / Cir(2^p2+2^p2);
    b[0] = 0;
    int delta,k;
    if (p1<p2) {
        delta = (p2-p1) % wordsize;
        k = p1 % wordsize;
    } else {
        delta = (p1-p2) % wordsize;
        k = p2 % wordsize;
    }
    for(int i=0;i<wordsize-1;i++) {
        b[(2*(i+1)*delta) % wordsize] =     b[(     2*i*delta)     % wordsize] ^
                                            a[((2*(i+1)*delta) +k) % wordsize] ^
                                            a[(((2*i+1)*delta) +k) % wordsize];
    }
}


void multiply_with_sum_of_power_of_twos(long *a, int wordsize, int p1, int p2) {
    long *h;
    // Cir(a) = Cir(a) * Cir(2^p1+2^p2)
    h = newBlock(wordsize);
    copy_and_multiply_with_sum_of_power_of_twos(h, a, wordsize, p1, p2);
    copyBlock(a, h, wordsize);
    free(h);
}

void divide_by_sum_of_power_of_twos(long *a, int wordsize, int p1, int p2) {
    // Cir(a) = Cir(a) / Cir(2^p2+2^p2);
    long *h;
    h = newBlock(wordsize);
    copy_and_divide_by_sum_of_power_of_twos(h, a, wordsize, p1, p2);
    copyBlock(a, h, wordsize);
    free(h);
}

void extendToArtinPrime(long *a, long *b, int wordsize) { 
    // extends a block of wordsize to wordsize+1 (= artin prime) for the circulent Boolean matrices
    copyBlock(a,b,wordsize);
    a[wordsize]=0;
}

void reduceToWordSize(long *a, long *b, int wordsize) {   
    // reduces a block of wordsize+1 (= artin prime) to wordsize for output
    copyBlock(a,b,wordsize);
    for(int i=0;i<wordsize;i++) {
        a[i] ^= b[wordsize];
    }
}


void testBlock(void) {
    // test of block operations
    printf("Testing Blocks ... \n");
    long *c, *b, *a;
    printf("Generating 10000 blocks of size 4096 and freeing memory ... \n");
    for(int i=0;i<10000;i++) {
        b = newBlock(4096);
        randomBlock(b, 4096);
        free(b);
    }
    printf("done ! \n");
    
    printf("zero Block\n");
    a = newBlock(3);
    randomBlock(a,3);
    zeroBlock(a,3);
    printBlock(a,3);
    printf("one Block\n");
    oneBlock(a,3);
    printBlock(a,3);
    free(a);
    
    printf("Testing Rotation of a block of size 3\n");
    a = newBlock(3);
    randomBlock(a,3);
    printBlock(a,3);
    for(int i=0;i<3;i++) {
        printf("Rotation by %d steps\n",i);
        rotateBlock(a,3,i);
        printBlock(a,3);
    }
    free(a);
    
    printf("Testing Xor two blocks of size 3\n");
    a = newBlock(3);
    b = newBlock(3);
    randomBlock(a,3);
    randomBlock(b,3);

    printf("a\n");
    printBlock(a,3);
    printf("b\n");
    printBlock(b,3);
    addBlock(a,b,3);
    printf("a xor b\n");
    printBlock(a,3);
    free(a);
    free(b);
    
    
    printf("Testing multiplying with sum of powers of two \n");
    a = newBlock(11);
    b = newBlock(11);
    printf("a\n");
    randomBlock(a,11);
    copyBlock(b,a,11);
    printBlock(a,11);
    printf("Cir(a) = Cir(a) * Cir(2^2+2^4)\n");
    multiply_with_sum_of_power_of_twos(a,11,2,4);
    printBlock(a,11);
    
    printf("Testing division modulo artin prime \n");
    divide_by_sum_of_power_of_twos(a,11,2,4);
    printBlock(a,11);
     
    printf("Xoring with a should give rows of 0 and 1\n");
    addBlock(a,b,11);
    printBlock(a,11);

    free(a);
    free(b);

    
    int artin = Artin[3];
    int ws = artin-1;
    printf("Testing multiplication from original wordsize of Artin prime %d and wordsize %d\n",artin,ws);
    printf("Data\n");
    a = newBlock(ws);
    b = newBlock(artin);
    c=newBlock(ws);

    randomBlock(a,ws);
    printBlock(a,ws);
    
    printf("Extended to artin prime %d \n",artin);
    extendToArtinPrime(b, a, ws);
    printBlock(b,artin);
    
    printf("Cir(a) = Cir(a) * Cir(2^6+2^7)\n");
    multiply_with_sum_of_power_of_twos(b,artin,6,7);
    printBlock(b,artin);

    printf("Cir(a) = Cir(a) / Cir(2^6+2^7)\n");
    divide_by_sum_of_power_of_twos(b,artin,6,7);
    printBlock(b,artin);
    
    printf("Reduce to word size %d\n", ws);
    reduceToWordSize(c,b,ws);
    printBlock(c,ws);
    
    printf("Xor with original data\n");
    addBlock(a,c,ws);
    printBlock(a,ws);
    
    free(a);
    free(b);
    free(c);
    
    printf("Test of Block succeeded\n");

}




int isArtin(int i) {
    // Tests whether a number is an Artin number by using the constant table Artin
    // see constants.h
    int upper, lower,test;
    lower=0;
    upper=ArtinSize-1;
    while (upper>lower) {
        test = (upper+lower)/2;
        if (Artin[test]>i) {
            upper=test-1;
            if (upper < 0) upper=0;
        } else  if (Artin[test]<i) {
            lower=test+1;
        } else {
            return Artin[test] == i;
        }
    }
    return Artin[upper] == i;
}

int ArtinSplitter(int *w, int wordsize, int maxSplitters) {
    // splits a word into Artin numbers -1 by using a table
    // WARNING: The splitting is only guarranteed to be optimal 
    //           (i.e. uses Artin numbers of maxmimum size)
    //          only if w is a power of two and w <= 2^ ArtinPowerSplitSize
    //               or w <= ArtinSplitSize 
    // The splitting defines the encoding. Improvements of this
    // functions may defunct existing encodings using strange large wordsizes
    //
    // Output is the number of summands
    // if the output is 0 then there is no splitting within maxSplitter summands
    
    if ((wordsize%2 !=0) || (wordsize ==0) || (maxSplitters==0) || (wordsize<0)) {
        return 0;
    } else {
        int wrest = wordsize;
        int index = 0;
        while ((wrest > 0) && (index < maxSplitters)) {
            if (isArtin(wrest+1)) {
                w[index] = wrest;
                wrest = 0;
                index ++;
            } else if (wrest <= ArtinSplitSize*2) {
                w[index] = ArtinSplit[wrest/2-1][1];
                wrest -= w[index];
                index ++;
            } else if (wrest <= ArtinPowerSplit[ArtinPowerSplitSize-1][0]) {
                int i = ArtinPowerSplitSize-1;
                while ((i>0)&&(ArtinPowerSplit[i][0]>wrest)) {
                    i--;
                }
                w[index] = ArtinPowerSplit[i][1];
                wrest -= w[index];
                index++;
            }
        }
        if (wrest > 0 ) {
            return 0;
        } else {
            for (int i=index; i< maxSplitters; i++) {
                w[i]=0;
            }
            return index;
        }
    }
}

void testArtin(void) {
    printf("Artin number Test\n");
    const int ws=10;
    int www[ws];
    for(int wort=62; wort<413100; wort+=6414){
        if (ArtinSplitter(www,wort,ws)) {
            printf("Splitting %d = ", wort);
            for (int i=0; i<ws; i++) {
                printf("%d + ",www[i]);
            }   
            printf(" nothing \n");
        } else {
            printf(" %d cannot be split into maximum 10 terms \n",wort);
        }
    }
    for(int i=1; i<200; i++) {
        if (isArtin(i)) {
            printf("%d is Artin number\n",i);
        }
    }
}


int X(int i, int artin) { 
    // computes the power of two for the data in the Cauchy matrix
    // i can have values 0,1,..., artin- #code
    return i;
}
int Y(int i, int artin) { 
    // computes the power of two for the code in the Cauchy matrix
    // i can have values 0,1,..., artin- #data
    return artin-i-1;
}


long *newVector(int wordsize, int n) {
    // produces new space for a vector of n blocks of wordsize integers
    // must be freed by free
    long *block;
    block = malloc(n * wordsize * sizeof(long));
    return block;

}

void addVector(long *data, long *operand, int wordsize, int n){
    for (int i=0; i< n*wordsize; i++) {
        data[i] ^= operand[i];
    }
}

void randomVector(long *data, int wordsize, int n) {
    // produces random vector for testing. 
    // Not initialized the series is always the same

    for (int i=0; i<n; i++) {
        randomBlock(&data[i*wordsize],wordsize);
    }
}
void vectorPrint(long *data, int wordsize, int n) {
    // prints the vector of the first bits
    // to verify the correctness of the mathematics
    printf("( ");
    for (int i=0; i<n; i++) {
        printFirstBit(&data[i*wordsize],wordsize);
        printf(" ");
    }
    printf(")\n");
}


void testVector() {
    // tests crucial properties of vectors
    long *v1;
    v1 = newVector(5,3);
    vectorPrint(v1,5,3);
    randomVector(v1,5,3);
    vectorPrint(v1,5,3);
    free(v1);
}

void innerCauchy(long *code, long *data, int artinprime, int n, int m) { 
// computes the Cauchy product 
// int data[n][artin]
// int code[n][artin]
    long *h, *s;
    s = newBlock(artinprime);
    h = newBlock(artinprime);

    for(int j=0;j<m;j++) {
        zeroBlock(&code[j*artinprime],artinprime);
        for(int i=0;i<n;i++) {
            copy_and_divide_by_sum_of_power_of_twos(h, &data[i*artinprime], artinprime, X(i,artinprime),Y(j,artinprime));
            addBlock(&code[j*artinprime],h,artinprime);
        } 
    }
    free(s);
    free(h);
}
void innerCauchyXY(long *code, long *data, int artinprime, int n, int m, int *Xs, int *Ys){
    // computes the Cauchy product for flexible array
    // int data[n][artin]
    // int code[n][artin]
    long *h, *s;
    s = newBlock(artinprime);
    h = newBlock(artinprime);
    
    for(int j=0;j<m;j++) {
        zeroBlock(&code[j*artinprime],artinprime);
        for(int i=0;i<n;i++) {
            copy_and_divide_by_sum_of_power_of_twos(h, &data[i*artinprime], artinprime, Xs[i],Ys[j]);
            addBlock(&code[j*artinprime],h,artinprime);
        } 
    }
    free(s);
    free(h);
}



void divide_by_A(long *result, int k, int artinprime, int n, int *Xs) { 
    // divides result by a_k   =  \prod_{i\neq k} (x_i + x_k)
    for(int i=0;i<n;i++) {
        if (i!=k) {
            divide_by_sum_of_power_of_twos(result, artinprime, Xs[i], Xs[k]);
        }
    } 
}
void divide_by_B(long *result, int k, int wordsize, int n, int *Ys) {
    // divides result by b_k =  \prod_{i\neq k} (y_i + y_k)
    for(int i=0;i<n;i++) {
        if (i!=k) {
            divide_by_sum_of_power_of_twos(result, wordsize, Ys[i],Ys[k]);
        }
    } 
}
void multiply_by_E(long *result, int k, int wordsize, int n, int *Xs, int *Ys) {
    // multiplies result by e_k = \prod_{i=1}^n (x_k+y_i) 
    for(int i=0;i<n;i++) {
        multiply_with_sum_of_power_of_twos(result, wordsize, Xs[k],Ys[i]);
    }
}
void multiply_by_F(long *result, int k, int wordsize, int n, int *Xs, int *Ys) {
    // multiplies result by f_k = \prod_{i=1}^n (x_i+y_k)
    for(int i=0;i<n;i++) {
        multiply_with_sum_of_power_of_twos(result, wordsize, Xs[i],Ys[k]);
    }
}

void innerInverseCauchy(long *data, long *code, int artinprime, int n, int *Xs, int *Ys) {
    // Given code computes the data (or its complement) corresponding to the code 
    // code int[artinprime * n] 
    // data int[artinprime * n]
    // int[n] Xs = powers of twos corresponding to the missing data
    // int[n] Ys = powers of twos corresponding of the given code
    long *h, *Cs;
    h = newBlock(artinprime);
    Cs = newVector(artinprime,n);
    for(int j=0;j<n;j++) {
        copyBlock    (&Cs[j*artinprime], &code[j*artinprime], artinprime);
        multiply_by_F(&Cs[j*artinprime], j, artinprime, n, Xs, Ys);
        divide_by_B  (&Cs[j*artinprime], j, artinprime, n, Ys);
    }
    for(int i=0;i<n;i++) {
        zeroBlock(&data[i*artinprime],artinprime);
    }
    for(int i=0;i<n;i++) {
        for(int j=0;j<n;j++) {
            copy_and_divide_by_sum_of_power_of_twos(h, &Cs[j*artinprime], artinprime, Ys[j], Xs[i]);
            addBlock(&data[i*artinprime],h,artinprime);
        } 
    }
    for(int i=0;i<n;i++) {
        multiply_by_E(&data[i*artinprime], i, artinprime, n, Xs, Ys);
        divide_by_A  (&data[i*artinprime], i, artinprime, n, Xs);
    }
    free(h);
    free(Cs);
}

void extendToArtinPrimeVector(long *cir, long *data, int wordsize, int n){
    // adds to each block a zero integer for Circulant matrices
    for(int i=0;i<n;i++) {
        extendToArtinPrime(&cir[i*(wordsize+1)], &data[i*wordsize], wordsize);
    }
}
void reduceToWordSizeVector(long *code, long *cir, int wordsize, int n) {
    // removes extra block by evaluating parity for each block
    for(int i=0;i<n;i++) {
        reduceToWordSize(&code[i*wordsize], &cir[i*(wordsize+1)], wordsize);
    }    
}

void innerReconstructData(long *data, long *code, int artinprime, int n, int m, int *Xindex, int *Ylist) {
    // reconstructs the data
    // data consists partially of the available data 
    // boolean array Xindex[n] indicates wether data is available X[i]=1 <=> data is there
    // Ylist is the list of the Y-indices
    // n is the full data length 
    // m is the number of given encodings
    
    int *Yindex;
    long *preparedCode;
    long *reconstructData;
    long *availableDataVector;
    long *missingDataVector;
    long *availableDataIndexList;
    int *availableDataX;
    long *missingDataIndexList;
    int *missingDataX;
    int nAvail;
    int nMissing;
    
    availableDataVector = newVector(artinprime,n);
    missingDataVector = newVector(artinprime,n);
    
    availableDataIndexList = newBlock(n);
    availableDataX = calloc(n, sizeof(int));
    
    missingDataIndexList = newBlock(n);
    missingDataX = calloc(n, sizeof(int));
    
    nAvail=0;
    nMissing =0;
    for (int i=0;i<n;i++) {
        if (Xindex[i]) {
            // gather available data and remember index
            copyBlock(&availableDataVector[nAvail*artinprime],&data[i*artinprime],artinprime);
            availableDataIndexList[nAvail]=i;
            // compute the corresponding powers of two
            availableDataX[nAvail] = X(i,artinprime);
            nAvail++;
        } else {
            // remember the missing data indices
            missingDataIndexList[nMissing]=i;
            missingDataX[nMissing] = X(i,artinprime);
            nMissing++;
        }
    }

    Yindex = calloc(nMissing, sizeof(int));
    for (int i=0;i<nMissing;i++) {
        Yindex[i] = Y(Ylist[i],artinprime);
    }
    
    
    preparedCode  = newVector(artinprime, nMissing);
    innerCauchyXY(preparedCode, availableDataVector, artinprime, nAvail, nMissing, availableDataX, Yindex);
    addVector(preparedCode, code, artinprime, nMissing);
    
    reconstructData = newVector(nMissing, artinprime);
    innerInverseCauchy(reconstructData, preparedCode,artinprime,nMissing, missingDataX, Yindex);
    for (int i=0;i<nMissing;i++) {
        copyBlock(&data[missingDataIndexList[i]*artinprime],&reconstructData[i*artinprime],artinprime);
    }
    
    free(availableDataVector);
    free(missingDataVector);
    
    free(availableDataIndexList);
    free(availableDataX);
    
    free(missingDataIndexList);
    free(missingDataX); 
    
    free(Yindex);
    free(preparedCode);
    free(reconstructData);
}

void testCauchy(void) {
    // prints out the Cauchy matrix

    int m=3;
    int n=2;
    int k;
    int artinprime = Artin[1];
    long *v1,*v2;
    v1 = newVector(artinprime,n);
    v2 = newVector(artinprime,m);
    for(int i=0; i<n;i++) {
        for(int j=0; j<artinprime;j++) {
            for( k=0;k<artinprime*n;k++){
                v1[k] = 0;
            }
            v1[i*artinprime] = 1;
            rotateBlock(&v1[i*artinprime], artinprime, j);
 //           vectorPrint(v1,artinprime,n); this will print the unit matrix
            innerCauchy(v2,v1,artinprime,n,m);
            vectorPrint(v2,artinprime,m);
        }
        printf("\n");
    }
    printf("Cauchy Result should be: \n");
    printf("00101 00011 01100\n");
    printf("01101 01110 00110\n");
    printf("01001 00111 00011\n");
    printf("01011 01100 01110\n");
    printf("01010 00110 00111\n\n");
    printf("00110 00111 01011\n");
    printf("00011 01100 01010\n");
    printf("01110 00110 00101\n");
    printf("00111 00011 01101\n");
    printf("01100 01110 01001\n\n");

    free(v1);
    free(v2);
    
    long *b, *a;
    b = newBlock(artinprime);
    a = newBlock(artinprime);
    k=0;
    const int nc=3;
    
    int Xs[nc];
    Xs[0]=X(0,artinprime);
    Xs[1]=X(1,artinprime);
    Xs[2]=X(2,artinprime);
    int Ys[nc];
    Ys[0]=Y(0,artinprime);
    Ys[1]=Y(1,artinprime);
    Ys[2]=Y(2,artinprime);
    
    printf("Test of A \n");
    for(int j=0; j<artinprime;j++) {
        oneBlock(b, artinprime);
        rotateBlock(b, artinprime, j);
//        printFirstBit(b,artinprime);
//        printf("\n");
        divide_by_A(b, k, artinprime, nc, Xs);
        printFirstBit(b,artinprime);
        
        //        (v1,artinprime,n); this will print the unit matrix
        printf("\n");
    }
    printf("Correct result\n");
    printf("01000\n");
    printf("00100\n");
    printf("00010\n");
    printf("00001\n");
    printf("01111\n\n");
    
    
    printf("Test of B \n");
    for(int j=0; j<artinprime;j++) {
        oneBlock(b, artinprime);
        rotateBlock(b, artinprime, j);
        //        printFirstBit(b,artinprime);
        //        printf("\n");
        divide_by_B(b, k, artinprime, nc, Ys);
        printFirstBit(b,artinprime);
        
        //        (v1,artinprime,n); this will print the unit matrix
        printf("\n");
    }
    printf("Correct result\n");
    printf("01000\n");
    printf("00100\n");
    printf("00010\n");
    printf("00001\n");
    printf("01111\n\n");

    printf("Test of E \n");
    for(int j=0; j<artinprime;j++) {
        oneBlock(b, artinprime);
        rotateBlock(b, artinprime, j);
        //        printFirstBit(b,artinprime);
        //        printf("\n");
        multiply_by_E(b, k, artinprime, nc, Xs,Ys);
        printFirstBit(b,artinprime);
        
        //        (v1,artinprime,n); this will print the unit matrix
        printf("\n");
    }
    printf("Correct result\n");
    printf("01010\n");
    printf("00100\n");
    printf("10010\n");
    printf("01001\n");
    printf("10100\n\n");

    printf("Test of F \n");
    for(int j=0; j<artinprime;j++) {
        oneBlock(b, artinprime);
        rotateBlock(b, artinprime, j);
        //        printFirstBit(b,artinprime);
        //        printf("\n");
        multiply_by_F(b, k, artinprime, nc, Xs,Ys);
        printFirstBit(b,artinprime);
        
        //        (v1,artinprime,n); this will print the unit matrix
        printf("\n");
    }
    printf("Correct result\n");
    printf("01001\n");
    printf("10100\n");
    printf("01010\n");
    printf("00101\n");
    printf("10010\n\n");

    free(a);
    free(b);

    artinprime=11;
    n=3;
    printf("Inverse Cauchy matrix\n");
    Xs[0]=X(0,artinprime);
    Xs[1]=X(1,artinprime);
    Xs[2]=X(2,artinprime);
    Ys[0]=Y(0,artinprime);
    Ys[1]=Y(1,artinprime);
    Ys[2]=Y(2,artinprime);

    v1 = newVector(artinprime,n);
    v2 = newVector(artinprime,n);
    for(int i=0; i<n;i++) {
        for(int j=0; j<artinprime;j++) {
            for( k=0;k<artinprime*n;k++){
                v1[k] = 0;
            }
            v1[i*artinprime] = 1;
            rotateBlock(&v1[i*artinprime], artinprime, j);
            //vectorPrint(v1,artinprime,n); //this will print the unit matrix
            innerInverseCauchy(v2,v1,artinprime,n,Xs,Ys);
            vectorPrint(v2,artinprime,m);
        }
        printf("\n");
    }
    printf("Result should be: \n");
    
    printf("00011111000 00101001010 00010010010\n");
    printf("00001111100 00010100101 00001001001\n");
    printf("00000111110 01110101101 01111011011\n");
    printf("00000011111 01000101001 01000010010\n");
    printf("01111110000 01011101011 00100001001\n");
    printf("00111111000 01010001010 01101111011\n");
    printf("00011111100 00101000101 01001000010\n");
    printf("00001111110 01101011101 00100100001\n");
    printf("00000111111 01001010001 01101101111\n");
    printf("01111100000 01011010111 01001001000\n");
    printf("00111110000 01010010100 00100100100\n\n");
  
    printf("01010010100 00100000100 01000000001\n");
    printf("00101001010 00010000010 01011111111\n");
    printf("00010100101 00001000001 01010000000\n");
    printf("01110101101 01111011111 00101000000\n");
    printf("01000101001 01000010000 00010100000\n");
    printf("01011101011 00100001000 00001010000\n");
    printf("01010001010 00010000100 00000101000\n");
    printf("00101000101 00001000010 00000010100\n");
    printf("01101011101 00000100001 00000001010\n");
    printf("01001010001 01111101111 00000000101\n\n");
  
    printf("01011010111 01000001000 01111111101\n");
    printf("01001001000 01111111101 00111111100\n");
    printf("00100100100 01000000001 00011111110\n");
    printf("00010010010 01011111111 00001111111\n");
    printf("00001001001 01010000000 01111000000\n");
    printf("01111011011 00101000000 00111100000\n");
    printf("01000010010 00010100000 00011110000\n");
    printf("00100001001 00001010000 00001111000\n");
    printf("01101111011 00000101000 00000111100\n");
    printf("01001000010 00000010100 00000011110\n");
    printf("00100100001 00000001010 00000001111\n");
    printf("01101101111 00000000101 01111111000\n\n");
    free(v1);
    free(v2);
    
    printf("Reconstructing from the first n Parities\n");
    long *data;
    artinprime = Artin[4];
    n=3;
    int nn=3;
    m=4;
    data = newVector(artinprime-1,n);
    randomVector(data,artinprime-1,n);
    printf("Starting data\n");
    vectorPrint(data,artinprime-1,n);
    
    long *artindata;
    artindata = newVector(artinprime,n);
    extendToArtinPrimeVector(artindata,data,artinprime-1,n);
    printf("Artin extended data\n");
    vectorPrint(artindata,artinprime,n);

    long *artincode;
    artincode = newVector(artinprime,m);
    innerCauchy(artincode,artindata,artinprime,n,m);
    printf("Artin code\n");
    vectorPrint(artincode,artinprime,m);

    int *XX, *YY;
    XX = calloc(sizeof(int),n);
    YY = calloc(sizeof(int),n);
    long *artinrecode;
    artinrecode = newVector(artinprime,nn);
    for (int i=0;i<n;i++){
        XX[i] = X(i,artinprime);
        YY[i] = Y(i,artinprime);
        copyBlock(&artinrecode[i*artinprime],&artincode[(i)*artinprime],artinprime);
    }
    printf("Chosen sub-code\n");
    vectorPrint(artinrecode,artinprime,nn);

    long *artinredata;
    artinredata = newVector(artinprime,nn);
    innerInverseCauchy(artinredata, artinrecode, artinprime, n, XX, YY);
    printf("Artin reconstructed data\n");
    vectorPrint(artinrecode,artinprime,nn);

    long *redata;
    redata = newVector(artinprime+1,nn);
    reduceToWordSizeVector(redata, artinredata, artinprime-1, n);
    printf("reconstructed data\n");

    vectorPrint(redata,artinprime-1,nn);
    
    free(redata);
    free(artinredata);
    free(artinrecode);
    free(artincode);
    free(artindata);
    free(data);
    free(XX);
    free(YY);
    
    
  
    printf("Reconstructing from any Parities\n");
    artinprime = Artin[4];
    n=4;
    m=4;
    data = newVector(artinprime-1,n);
    randomVector(data,artinprime-1,n);
    printf("Starting data\n");
    vectorPrint(data,artinprime-1,n);
        
    artindata = newVector(artinprime,n);
    extendToArtinPrimeVector(artindata,data,artinprime-1,n);
    printf("Artin extended data\n");
    vectorPrint(artindata,artinprime,n);
    
    artincode = newVector(artinprime,m);
    innerCauchy(artincode,artindata,artinprime,n,m);
    printf("Artin code\n");
    vectorPrint(artincode,artinprime,m);
    
    int *Xindex;
    XX = calloc(sizeof(int), n);
    Xindex = calloc(sizeof(int), n);
    YY = calloc(sizeof(int), 3);
    artinrecode = newVector(artinprime,m);
    artinredata = newVector(artinprime,n);
    
    copyBlock(&artinrecode[0 * artinprime], &artincode[0 * artinprime],artinprime);
    copyBlock(&artinrecode[1 * artinprime], &artincode[1 * artinprime],artinprime);
    copyBlock(&artinrecode[2 * artinprime], &artincode[2 * artinprime],artinprime);
//    copyBlock(&artinrecode[3 * artinprime], &artincode[3 * artinprime],artinprime);
    YY[0] = 0;
    YY[1] = 1;
    YY[2] = 2;
//    YY[3] = 3;
    Xindex[0] = 1;
    Xindex[1] = 0;
    Xindex[2] = 0;
    Xindex[3] = 0;
    copyBlock(artinredata,artindata,artinprime*n);
//    zeroBlock(&artinredata[0* artinprime],artinprime);
    zeroBlock(&artinredata[1* artinprime],artinprime);
    zeroBlock(&artinredata[2* artinprime],artinprime);
    zeroBlock(&artinredata[3* artinprime],artinprime);

    printf("Available data\n");
    vectorPrint(artinredata,artinprime,4);
    printf("Chosen sub-code\n");
    vectorPrint(artinrecode,artinprime,4);
    
    innerReconstructData(artinredata, artinrecode, artinprime, n, 4,  Xindex, YY);
    printf("Artin reconstructed data\n");
    vectorPrint(artinredata,artinprime,n);
    
    redata = newVector(artinprime-1,n);
    reduceToWordSizeVector(redata, artinredata, artinprime-1, n);
    printf("reconstructed data\n");
    
    vectorPrint(redata,artinprime-1,n);
    
    free(redata);
    free(artinredata);
    free(artinrecode);
    free(artincode);
    free(artindata);
    free(data);
    free(XX);
    }

void CauchyMDSArtin(long *code, long *data, int wordsize, int n, int m) {
    // given the data with wordsize = Artin prime -1
    // we do not test whether wordsize+1 is an Artin prime
    
    long *dataCir;
    long *codeCir;
    
    dataCir = newVector(wordsize+1,n);
    codeCir = newVector(wordsize+1,m);
    extendToArtinPrimeVector(dataCir, data, wordsize,n );
    innerCauchy(codeCir, dataCir, wordsize+1, n, m);
    reduceToWordSizeVector(code, codeCir, wordsize,m);

    free(dataCir);
    free(codeCir);
}
void reconstructCauchyMDSArtin(long *data, long *code, 
                               int wordsize, 
                               int n, int m, 
                               int *Xindex, int *Ylist) {
    // adds missing data to *data 
    // wordsize is Artin prime-1 (will be trusted)
    // boolean array Xindex[n] indicates wether data is available X[i]=1 <=> data is there
    // Ylist is the list of the Y-indices
    // n is the full data length 
    // m is the number of given encodings

    long *dataCir;
    long *codeCir;
    int artinprime=wordsize+1;
    dataCir = newVector(artinprime,n);
    codeCir = newVector(artinprime,m);
    extendToArtinPrimeVector(dataCir, data, wordsize,n );
    extendToArtinPrimeVector(codeCir, code, wordsize,m );
    
    innerReconstructData(dataCir, codeCir, artinprime, n, m, Xindex, Ylist);
    reduceToWordSizeVector(data, dataCir, wordsize,n);
    free(dataCir);
    free(codeCir);
}

void testMDSARtin(void) {
    // tests the MDS Codes 
    long *code;
    long *data;
    int wordsize = Artin[6]-1;
    int n=4;
    int m=3;
    data = newVector(wordsize,n);
    code = newVector(wordsize, m);
    randomVector(data,wordsize,n);
    printf("Input data\n");
    vectorPrint(data, wordsize, n);
    CauchyMDSArtin(code, data, wordsize, n, m);
    printf("Code\n");
    vectorPrint(code, wordsize, m);
    
    int *XX,*Xindex, *YY;
    long *recode, *redata;
    XX = calloc(sizeof(int),n);
    Xindex = calloc(sizeof(int),n);
    YY = calloc(sizeof(int),3);
    recode = newVector(wordsize,m);
    redata = newVector(wordsize,n);
    
    copyBlock(&recode[0 * wordsize], &code[0 * wordsize],wordsize);
    copyBlock(&recode[1 * wordsize], &code[1 * wordsize],wordsize);
    copyBlock(&recode[2 * wordsize], &code[2 * wordsize],wordsize);
    //copyBlock(&recode[3 * wordsize], &code[3 * wordsize],wordsize);
    YY[0] = 0;
    YY[1] = 1;
    YY[2] = 2;
    //YY[3] = 3;
    Xindex[0] = 1;
    Xindex[1] = 0;
    Xindex[2] = 0;
    Xindex[3] = 0;
    copyBlock(redata,data,wordsize*n);
    //zeroBlock(&redata[0* wordsize],wordsize);
    zeroBlock(&redata[1* wordsize],wordsize);
    zeroBlock(&redata[2* wordsize],wordsize);
    zeroBlock(&redata[3* wordsize],wordsize);
    
    printf("Available data\n");
    vectorPrint(redata,wordsize,n);
    printf("Chosen sub-code\n");
    vectorPrint(recode,wordsize,3);
    
    reconstructCauchyMDSArtin(redata, recode, wordsize, n, 3,  Xindex, YY);
    printf("Reconstructed data\n");
    vectorPrint(redata,wordsize,n);
    
    free(XX); free(Xindex);free(YY);free(recode);free(redata);
    free(code); free(data);

    
    
}



void CauchyMDS(long *code, long *data, int wordsize, int n, int m) { 
    // produces m Boolean Circulant Cauchy blocks of size wordsize
    // code must exist as an int array [0..wordsize*m] and will be overwritten
    
    long *subCodeVector[40];
    long *subDataVector[40];
    int subWordsize[40];
    int subWordNumber;
    int position;
    
    subWordNumber = ArtinSplitter(subWordsize, wordsize, 40);
    
    position =0;
    for(int i=0;i<subWordNumber;i++) {
        subDataVector[i] = newVector(subWordsize[i],n);
        subCodeVector[i] = newVector(subWordsize[i],m);
        
        for(int k=0;k<n;k++) {
            copyBlock(subDataVector[i]+k*subWordsize[i], 
                      &data[k*wordsize+position] , 
                      subWordsize[i]);  
        }
        position += subWordsize[i];
    }
    for(int i=0;i<subWordNumber;i++) {
        CauchyMDSArtin(subCodeVector[i], subDataVector[i], subWordsize[i], n, m);
    }
    position =0;
    for(int i=0;i<subWordNumber;i++) {
        for(int k=0;k<m;k++) {
            copyBlock(&code[k*wordsize+position], 
                      subCodeVector[i]+k*subWordsize[i] , 
                      subWordsize[i]);  
        }
        position += subWordsize[i];
    }
    for(int i=0;i<subWordNumber;i++) {
        free(subDataVector[i]);
        free(subCodeVector[i]);
    }
}

void reconstructCauchyMDS(long *data, long *code, 
                          int wordsize, 
                          int n, int m, 
                          int *Xindex, int *Ylist) {  
    // data is a vector where the missing entries will be completed
    // data i*wordsize ... (i+1)*wordsize-1  = i-th data block for i=0 ... n
    // missing entries are marked with Xindex[i] == 0 
    // existing ones are marked with Xindex[i] == 1
    // code is a vector of available code blocks
    // Ylist[i] denotes the position of the code blocks and m its number
    //      if the 2nd and the 4th code block are availabe
    //      Ylist[0]=1; Ylist[1]=3; m=2
    // 
    // reconstructCauchyMDS does not check whether the data can be reconstructed
    // if it is not possible it produces an interesting pseudo random sequence
    // otherwise it reconstructs the data by overwriting the missing entries
    
    int subWordsize[40];
    int subWordNumber;    
    long * subCodeVector[40];
    long * subDataVector[40];
    int position;

    // First Artin Splitting
    subWordNumber = ArtinSplitter(subWordsize, wordsize, 40);

    // Prepare Matrices for each artin number
    position =0;
    for(int i=0;i<subWordNumber;i++) {
        subDataVector[i] = newVector(subWordsize[i],n);
        subCodeVector[i] = newVector(subWordsize[i],m);
        
       for(int k=0;k<n;k++) {
            copyBlock(subDataVector[i]+k*subWordsize[i], 
                      &data[k*wordsize+position] , 
                      subWordsize[i]);  
        }
        for(int k=0;k<m;k++) {
            copyBlock(subCodeVector[i]+k*subWordsize[i], 
                      &code[k*wordsize+position], 
                      subWordsize[i]);  
        }
        
        position += subWordsize[i];
    }
    // reconstructCircularCauchyMDSArtin
    for(int i=0;i<subWordNumber;i++) {
        reconstructCauchyMDSArtin(subDataVector[i], subCodeVector[i], 
                                  subWordsize[i],
                                  n, m, Xindex, Ylist);
    }
    position =0;
    for(int i=0;i<subWordNumber;i++) {        
        for(int k=0;k<n;k++) {
            copyBlock(&data[k*wordsize+position],
                      subDataVector[i]+k*subWordsize[i], 
                      subWordsize[i]);  
        }
//wtf        for(int k=0;k<m;k++) {
// wtf           copyBlock(subCodeVector[i]+k*subWordsize[i], 
// wtf                     &code[k*wordsize+position], 
// wtf                     subWordsize[i]);  
//  wtf      }
        position += subWordsize[i];
    }

 /*   for(int i=0;i<subWordNumber;i++) {        
        free(subCodeVector[i]);
        free(subDataVector[i]);
    }*/
}
void testMDS(void) {
    long *code;
    long *data;
    int wordsize = 64;
    int www[40];
    if (ArtinSplitter(www,wordsize,40)) {
        printf("Splitting %d = ", wordsize);
        for (int i=0; i<ArtinSplitter(www,wordsize,40); i++) {
            printf("%d + ",www[i]);
        }   
        printf(" nothing \n");
    } else {
        printf(" %d cannot be split into maximum 10 terms \n",wordsize);
    }
    

    int n=4;
    int m=5;
    data = newVector(wordsize,n);
    code = newVector(wordsize, m);
    randomVector(data,wordsize,n);
    printf("Input data\n");
    vectorPrint(data, wordsize, n);
    CauchyMDS(code, data, wordsize, n, m);
    printf("Code\n");
    vectorPrint(code, wordsize, m);
    
    int *Xindex, *YY;
    long *recode, *redata;
    Xindex = calloc(sizeof(int),n);
    YY = calloc(sizeof(int),3);
    recode = newVector(wordsize,m);
    redata = newVector(wordsize,n);
    
    //copyBlock(&recode[0 * wordsize], &code[0 * wordsize],wordsize);
    copyBlock(&recode[0 * wordsize], &code[1 * wordsize],wordsize);
    copyBlock(&recode[1 * wordsize], &code[2 * wordsize],wordsize);
    copyBlock(&recode[2 * wordsize], &code[3 * wordsize],wordsize);
    YY[0] = 1;
    YY[1] = 2;
    YY[2] = 3;
    //YY[3] = 3;
    Xindex[0] = 0;
    Xindex[1] = 1;
    Xindex[2] = 0;
    Xindex[3] = 0;
    copyBlock(redata,data,wordsize*n);
    zeroBlock(&redata[0* wordsize],wordsize);
    //zeroBlock(&redata[1* wordsize],wordsize);
    zeroBlock(&redata[2* wordsize],wordsize);
    zeroBlock(&redata[3* wordsize],wordsize);
    
    printf("Available data\n");
    vectorPrint(redata,wordsize,n);
    printf("Chosen sub-code\n");
    vectorPrint(recode,wordsize,3);
    
    reconstructCauchyMDS(redata, recode, wordsize, n, 3,  Xindex, YY);
    printf("Reconstructed data\n");
    vectorPrint(redata,wordsize,n);
    
    free(Xindex);free(YY);free(recode);free(redata);free(code);free(data);
}

const int MAX_FILE_LENGTH= 200;

void createDataFileNames (char * datafilename,  int n, const char * origFileName) {    
        sprintf(datafilename,"%s%s%d%s",origFileName,".D",n,".mds");
        printf("%s\n",datafilename);
}

void analyzeFilename(char filename[], char masterfile[], int * fitsScheme, int * isData, long * index, int * n, int * blocksize, int * word) {    // test whether to encode or decode - indicated by the suffix mds
    int tfits, tindex, tn, tblock, tword, tisdata;
    tfits = 0; 
    if ((strlen(  filename) >=4 ) &&
        (strncmp( filename + strlen(filename)-4 ,".mds",4)==0)) { 
        // Now we try to read out all the decoding parameters
        // FORMAT    Original-File-Name.D4N8B1024W8.mds
        // means 4th out of 8 data file with block size 1024 on 8x8 wordsize encoding
        // FORMAT    Original-File-Name.P4N8B1024W8.mds
        // means 4th out of 8 parity file with block size 1024 on 8x8 wordsize encoding
        int j = (int) strlen(filename)-5;
        if ((j>0) && (filename[j] >= '0') && (filename[j]<='9')) {
            tword = filename[j]-'0';
            j--;
            if ((j>0) && (filename[j] == 'W')) {
//                printf("64 Bit parallel code found \n");
                j--;
                tblock=0;
                int factor = 1;
                while((j>0) && (filename[j] >= '0') && (filename[j] <= '9')) {
                    tblock = tblock + factor * (int) (filename[j] -'0');
                    factor *=10;
                    j--;
                }
//                printf("Block size is %d\n", tblock);
                if ((j>0) && (filename[j] == 'B')) {
                    j--;
                    tn=0;
                    factor = 1;
                    while((j>0) && (filename[j] >= '0') && (filename[j] <= '9')) {
                        tn = tn + factor * (int) (filename[j] -'0');
                        factor *=10;
                        j--;
                    }
//                    printf("Number of data words is %d\n", tn);
                    if ((j>0) && (filename[j] == 'N')) {
                        j--;
                        tindex=0;
                        factor = 1;
                        while((j>0) && (filename[j] >= '0') && (filename[j] <= '9')) {
                            tindex = tindex + factor * (int) (filename[j] -'0');
                            factor *=10;
                            j--;
                        }
//                        printf("Index is %d\n", tindex);
                        if (j>0) {
                            if (filename[j] == 'D') {
                                tisdata = 1;
                                tfits = 1;
//                                printf("It is a data word\n");
                            } else if (filename[j] == 'P') {
                                tisdata = 0;
                                tfits = 1;
//                                printf("It is a parity word\n");
                            } else if (filename[j] == 'S') {
                                tisdata = 2;
                                tfits =1;
//                                printf("It is a size indicating file\n");                                
                            } else {
//                                printf("I do not know what it is \n");
                                tfits = 0;
                            }
                            j--;
                            if (tfits && (j>0) && (filename[j] == '.')) {
//                                printf("All is fine\n");
                                memcpy(masterfile, filename,j+1);
                                masterfile[j] = '\0';
//                                printf("Masterfile :%s: \n",masterfile);
                            }
                        }
                    }
                }
            }
        }
    }
    *fitsScheme = tfits;
    if (tfits) {
        *isData     = tisdata;
        *index      = tindex;
        *n          = tn;
        *blocksize  = tblock;
        *word       = tword;
    }
}

int main (int argc, const char * argv[])
{
//   testBlock();
//    testArtin(); 
//    testVector();
//   testCauchy();
//    testMDS();
    char * filename;
    char * masterfile;
    char * tmpmasterfile;
    int tn;
    int tword;
    int tblocksize;
    int fitsScheme, isData, n, blocksize, word;
    long index;
    
    
    if (argc<=1) {
/*        testBlock();
        testArtin(); 
        testVector();
        testCauchy();
        testMDS();*/
        printf("Circulant Cauchy Codes Version 0.3 25.05.2012\n");
        printf("Copyright Christian Schindelhauer\n");
        printf("University of Freiburg, Germany\n\n");
        printf("Usage:\n");
        printf("CirCauchy sourcefile [n [m [w]]]\n");
        printf("     -> produces n+m MDS-files of with word length w \n");
        printf("CirCauchy file1.mds file2.mds ... filen.mds\n");
        printf("     -> reconstructs the source file from any combination of n MDS code\n");
        return 0;
    }
    filename = malloc(strlen(argv[1]));
    masterfile  = malloc(strlen(argv[1]));
    strcpy(filename,argv[1]); 
    analyzeFilename(filename, masterfile, &fitsScheme, 
                    &isData, &index, &n, &blocksize, &word);
    if (fitsScheme) {  // Reconstructing
        int *dataIndex;  // boolean denotes which data files are present
        dataIndex = calloc(sizeof(int),n);
        int *dataArgIndex;  // which argument is used as file name
        dataArgIndex = calloc(sizeof(int),n);
        FILE *dataFileIndex[1024]; // which file corresponds to the index
        int nAvailableData = 0; // how many files do we have
        int m=0; // number of available parity files
        int *ParityIndex;  // index of available parity file
        ParityIndex = calloc(sizeof(int),argc);
        FILE *parityFileIndex[1024]; // which file corresponds to the index
        int *parityArgIndex;  // which argument is used as file name
        parityArgIndex = calloc(sizeof(int),argc);       
        int fileSizeIndicatorFound = 0;
        long numberOfBytes;
        
        for (int i=1; i<argc; i++) {
            free(filename); 
            filename = malloc(strlen(argv[i]));
            tmpmasterfile  = malloc(strlen(argv[i]));
            strcpy(filename,argv[i]); 
            analyzeFilename(filename, masterfile, &fitsScheme, 
                            &isData, &index, &tn, &tblocksize, &tword);
            fitsScheme =  fitsScheme && ((n == tn) && (tblocksize == blocksize) && (tword == word));
            free(tmpmasterfile);
            if (fitsScheme) {
                if ((isData ==1) && (index>=0) && (index<n) && !dataIndex[index]) {
                    dataIndex[index] = 1;
                    dataArgIndex[index] = i;
                    nAvailableData ++;
                } else if (isData == 0) { // is parity
                    ParityIndex[m] = (int) index;
                    parityArgIndex[m] = i;
                    m ++;
                } else if (isData == 2) { // file size indicator found
                    fileSizeIndicatorFound = 1;
                    numberOfBytes = index;
                }
            }
        }
       if (nAvailableData+m < n) {
            printf("Only %d data and %d parity files found. %d are necessary.\n",nAvailableData,m,n);
            exit(0);
        }
        m = n-nAvailableData;
        for (int i=0; i<n; i++) {
            if (dataIndex[i]) {
                printf("Opening %s\n",argv[dataArgIndex[i]]);
                dataFileIndex[i] = fopen(argv[dataArgIndex[i]],"rb");
                if (dataFileIndex[i] == NULL) {
                    printf("Cannot open %s\n",argv[dataArgIndex[i]]);
                    return 1;
                }
            }
        }
        for (int i=0; i<m; i++) {
            printf("Opening %s\n",argv[parityArgIndex[i]]);
            parityFileIndex[i] = fopen(argv[parityArgIndex[i]],"rb");
            if (parityFileIndex[i] == NULL) {
                printf("Cannot open %s\n",argv[parityArgIndex[i]]);
                return 1;
            }
        }
        // Open Output File
        FILE *OutputFile;
        OutputFile = fopen(masterfile,"wb");
        if (OutputFile==NULL) {
            printf("Cannot open file %s\n",masterfile);
            exit(1);
        }
        printf("Open file %s\n",masterfile);        
        // Reconstruct
        int EOFfound = 0;
        long position = 0;
        long maxRead;
        long readbytes;
        long *datavector;
        long *parityvector;
        datavector   = newVector(blocksize,n);
        parityvector = newVector(blocksize,m);
        
        while (!EOFfound && (!fileSizeIndicatorFound || (position<numberOfBytes))) {
            maxRead = sizeof(long)*blocksize;
            
            zeroBlock(datavector,   blocksize*n);
            zeroBlock(parityvector, blocksize*m);
            for (int i=0; i<n; i++) {
                if (dataIndex[i]) {
                    readbytes = fread(&datavector[i*blocksize], 1 , sizeof(long)*blocksize, dataFileIndex[i]);
                    position += readbytes; 
                    if (readbytes < maxRead) {
                        maxRead = 0;
                        EOFfound = 1;
                    }
                }
            }
            for (int i=0; i<m; i++) {
                readbytes = fread(&parityvector[i*blocksize], 1 , sizeof(long)*blocksize,parityFileIndex[i]);
                position += readbytes; 
                if (readbytes < maxRead) {
                    maxRead = 0;
                    EOFfound = 1;
                }
            }
            reconstructCauchyMDS(datavector, parityvector, blocksize, n, m,  dataIndex, ParityIndex);
            //vectorPrint(datavector, blocksize,n);
            if (fileSizeIndicatorFound && (position > numberOfBytes)) {
                fwrite(datavector, 1, numberOfBytes % (sizeof(long)* blocksize*n) , OutputFile);
            } else {
                fwrite(datavector, 1, sizeof(long)* blocksize*n , OutputFile);
            }
            printf("\b\bProcessing %ld Bytes %ld Blocks\r",position *  sizeof(long), position/blocksize/n);
        }
        printf("Processed %ld Bytes %ld Blocks\n",position *  sizeof(long), position/blocksize/n);
        
        // Close Files
        fclose(OutputFile);
        printf("Closing file %s\n",masterfile);

        
        
        for (int i=0; i<n; i++) {
            if (dataIndex[i]) {
                fclose(dataFileIndex[i]);       
                printf("Closing %s\n",argv[dataArgIndex[i]]);
            }
        }  

        for (int i=0; i<m; i++) {
            fclose(parityFileIndex[i]);
            printf("Closing %s\n",argv[parityArgIndex[i]]);
            
        }          
        // decoding END
        free(dataIndex); 
        free(dataArgIndex); 
        free(parityArgIndex); free(ParityIndex);
        return 0;
    } else { // encoding
        FILE * datafile[1024];
        FILE * parityfile[1024];
        char * datafilename[1024];
        char * parityfilename[1024];

        
        FILE *InputFile;
        InputFile = fopen(argv[1],"rb");
        if (InputFile==NULL) {
            printf("Cannot open file %s\n",argv[1]);
            exit(1);
        }
        printf("Open file %s\n",argv[1]);
        
        long *inputvector;
        int n,m, wordsize;
        if ((argc > 2) && (atoi(argv[2]) > 0 )) {
            n=atoi(argv[2]);
        } else {
            n = 4;
        }
        if ((argc > 3) && (atoi(argv[3]) > 0 )) {
            m=atoi(argv[3]);
        } else {
            m = 4;
        }
        if ((argc > 4) && (atoi(argv[4]) > 0 )) {
            wordsize=atoi(argv[4]);
            if (wordsize%2!=0) {
                wordsize = 1024;
            }
        } else {
            wordsize = 1024;
        }

        inputvector = newVector(wordsize,n);
        for (int i=0; i<n; i++) {
            datafilename[i]=malloc(strlen(argv[1])+40);
            sprintf(datafilename[i],"%s.D%dN%dB%dW%d.mds",argv[1],i,n,wordsize, (int) sizeof(long));
            printf("Opening %s\n",datafilename[i]);
            datafile[i] = fopen(datafilename[i],"wb");
            if (datafile[i]==NULL) {
                printf("Cannot open %s\n",datafilename[i]);
                return 1;
            }
        }
        for (int i=0; i<m; i++) {
            parityfilename[i]=malloc(strlen(argv[1])+40);
            sprintf(parityfilename[i],"%s.P%dN%dB%dW%d.mds",argv[1],i,n,wordsize, (int) sizeof(long));
            printf("Opening %s\n",parityfilename[i]);
            parityfile[i] = fopen(parityfilename[i],"wb");
            if (parityfile[i]==NULL) {
                printf("Cannot open %s\n",parityfilename[i]);
                return 1;
            }
        }
        
        unsigned long lastReadSize= sizeof(long)*wordsize*n;
        long *data;
        data = newVector(wordsize,n);
        long *parity;
        parity = newVector(wordsize,m);
        
        long position = 0;
        while(lastReadSize == sizeof(long)*wordsize*n) {
            zeroBlock(data,wordsize*n);
            lastReadSize = fread(data, 1, sizeof(long)* wordsize*n , InputFile);
            printf("Processing %ld Bytes %ld Blocks\r",position, position/wordsize/n);
            position += lastReadSize;
            for (int i=0; i<n; i++) {
                fwrite(&data[i*wordsize], sizeof(long), wordsize, datafile[i]);
            }
            CauchyMDS(parity, data, wordsize, n, m);
            for (int i=0; i<m; i++) {
                fwrite(&parity[i*wordsize], sizeof(long), wordsize, parityfile[i]);
            }            
        }
        
        for (int i=0; i<n; i++) {
            fclose(datafile[i]);       
            printf("Closing %s\n",datafilename[i]);
            free(datafilename[i]);
            
        }  
        
        for (int i=0; i<m; i++) {
            fclose(parityfile[i]);
            printf("Closing %s\n",parityfilename[i]);
            free(parityfilename[i]);
            
        }  
        fclose(InputFile);
        
        FILE * sizeFile;
        char sizeFileName[1024]; 
        sprintf(sizeFileName,"%s.S%ldN%dB%dW%ld.mds",argv[1],position,n,wordsize, sizeof(long));
        printf("Opening %s\n",sizeFileName);
        sizeFile = fopen(sizeFileName,"wb");
        fclose(sizeFile);

        return 0;
    }
    return 0;
    printf("Enjoy the result!\n");
}

