/* 
 * trans.c - Matrix transpose B = A^T
 * 
 * Each transpose function must have a prototype of the form:
 * void trans(int M, int N, int A[N][M], int B[M][N]);
 *
 * A transpose function is evaluated by counting the number of misses
 * on a 1KB direct mapped cache with a block size of 32 bytes.
 *
 * Nathan Walzer, Lucas Varella
 * nwalzer, lnvarella
 * Group: lnvarella-nwalzer
 */ 
#include <stdio.h>
#include "cachelab.h"

int is_transpose(int M, int N, int A[N][M], int B[M][N]);

/* 
 * transpose_submit - This is the solution transpose function that you
 *     will be graded on for Part B of the assignment. Do not change
 *     the description string "Transpose submission", as the driver
 *     searches for that string to identify the transpose function to
 *     be graded. 
 */
char transpose_submit_desc[] = "Transpose submission";
void transpose_submit(int M, int N, int A[N][M], int B[M][N])
{
    int row, col, rowBlock, colBlock, temp, diagIdx, stride;
    if(N == 32 && M == 32){
	//treat 32x32 in 8x8 blocks
	stride = 8;
	for(rowBlock = 0; rowBlock < N; rowBlock += stride){
	    for(colBlock = 0; colBlock < M; colBlock += stride){
		for(row = rowBlock; row < rowBlock+stride; row++){
		    for(col = colBlock; col < colBlock+stride; col++){
			if(row != col){
			    B[col][row] = A[row][col];
			} else {
			    //if we have a diagonal we do not need to immediately transpose it
			    temp = A[row][row];
			    diagIdx = row;
			}
		    }
		    if(rowBlock == colBlock){
		        B[diagIdx][diagIdx] = temp;
		    }
		}
	    }
	}
    } else if (N == 64 && M == 64){
	int xOffset8x8, yOffset8x8, tempXOff, tempYOff;
	//treating this as 8x8 blocks, but splitting those into 4x4
	for(rowBlock = 0; rowBlock < N; rowBlock += 8){
	    for(colBlock = 0; colBlock < M; colBlock += 8){
		xOffset8x8 = colBlock-rowBlock;
		yOffset8x8 = rowBlock-colBlock;
		
		//Put first and second rows of A into third and fourth rows of B (avoid 1->1 and 2->2 evictions)
		tempXOff = xOffset8x8+2;
		for(row=rowBlock; row<(rowBlock+2); row++){
                    for(col=colBlock; col<(colBlock+8); col++){
                    	B[row+tempXOff][col+yOffset8x8]=A[row][col];
                    }
            	}

		//Put second and fourth rows of A into first and second rows of B (avoid 3->3 and 4->4 evictions)
		tempXOff = xOffset8x8 - 2;		
		for(row=rowBlock+2; row<(rowBlock+4); row++){
                    for(col=colBlock; col<(colBlock+8); col++){
                    	B[row+tempXOff][col+yOffset8x8]=A[row][col];
                    }
            	}

		//1-to-1 translation of A onto B (for given 8x8) (top left quad)
		tempXOff = xOffset8x8 + 2;		
		for(row=rowBlock; row<(rowBlock+2); row++){
                    for(col=colBlock; col<(colBlock+8); col++){
			temp = B[row+xOffset8x8][col+yOffset8x8];
                    	B[row+xOffset8x8][col+yOffset8x8]=B[row+tempXOff][col+yOffset8x8];
			B[row+tempXOff][col+yOffset8x8] = temp;
                    }
            	}
		
		//transpose upper left quad (into final position)
		for(row=rowBlock+xOffset8x8; row<(rowBlock+xOffset8x8+4); row++){
                    for(col=colBlock+yOffset8x8; col<(colBlock+yOffset8x8+4); col++){
                    	if((col-(colBlock+yOffset8x8))<(row-(rowBlock+xOffset8x8))){
                            temp=B[row][col];
                            B[row][col]=B[col+xOffset8x8][row+yOffset8x8];
                            B[col+xOffset8x8][row+yOffset8x8]=temp;
                        }
                    }
            	}

		//transpose upper right quad (right of prev quad final position)
		tempXOff = xOffset8x8-4;
		tempYOff = yOffset8x8+4;
		for(row=rowBlock+xOffset8x8; row<(rowBlock+xOffset8x8+4); row++){
                    for(col=colBlock+yOffset8x8+4; col<(colBlock+yOffset8x8+8); col++){
                    	if((col-(colBlock+yOffset8x8+4))<(row-(rowBlock+xOffset8x8))){
                            temp=B[row][col];
                            B[row][col]=B[col+tempXOff][row+tempYOff];
                            B[col+tempXOff][row+tempYOff]=temp;
                        }
                    }
            	}

		//Swap upper left's rows 1&2 with its own 3&4 respectively
		tempXOff = xOffset8x8+2;
		for(row=rowBlock; row<(rowBlock+2); row++){
                    for(col=colBlock+4; col<(colBlock+8); col++){
                    	temp=B[row+xOffset8x8][col+yOffset8x8];
                    	B[row+xOffset8x8][col+yOffset8x8]=B[row+tempXOff][col+yOffset8x8];
                    	B[row+tempXOff][col+yOffset8x8]=temp;
                    }
            	}
		
		//put fifth and sixth rows of A into 7 and 8 in B (avoids evictions
		tempXOff = xOffset8x8+2;
            	for(row=rowBlock+4; row<(rowBlock+6); row++){
                    for(col=colBlock; col<(colBlock+8); col++){
                    	B[row+tempXOff][col+yOffset8x8]=A[row][col];
                    }
		}

		//put A7 and 8 into B 5 and 6 (avoids evictions)
		tempXOff = xOffset8x8-2;
            	for(row=rowBlock+6; row<(rowBlock+8); row++){
                    for(col=colBlock; col<(colBlock+8); col++){
                    	B[row+tempXOff][col+yOffset8x8]=A[row][col];
                    }
            	}

            	//swap B 5 and 6 with 7 and 8 respectively
		tempXOff = xOffset8x8+2;
            	for(row=rowBlock+4; row<(rowBlock+6); row++){
                    for(col=colBlock; col<(colBlock+8); col++){
                    	temp=B[row+xOffset8x8][col+yOffset8x8];
                    	B[row+xOffset8x8][col+yOffset8x8]=B[row+tempXOff][col+yOffset8x8];
                    	B[row+tempXOff][col+yOffset8x8]=temp;
                    }
            	}

            	//transpose left quadrant of lower half
		tempXOff = xOffset8x8+4;
		tempYOff = yOffset8x8-4;
            	for(row=rowBlock+4+xOffset8x8; row<(rowBlock+8+xOffset8x8); row++){
                    for(col=colBlock+yOffset8x8; col<(colBlock+4+yOffset8x8); col++){
                    	if((col-(colBlock+yOffset8x8))<(row-(rowBlock+4+xOffset8x8))){
                            temp=B[row][col];
                            B[row][col]=B[col+tempXOff][row+tempYOff];
                            B[col+tempXOff][row+tempYOff]=temp;
                    	}
                    }
            	}

            	//transpose right quadrant of lower half
            	for(row=rowBlock+4+xOffset8x8; row<(rowBlock+8+xOffset8x8); row++){
                    for(col=colBlock+4+yOffset8x8; col<(colBlock+8+yOffset8x8); col++){
                    	if((col-(colBlock+4+yOffset8x8))<(row-(rowBlock+4+xOffset8x8))){
                            temp=B[row][col];
                            B[row][col]=B[col+xOffset8x8][row+yOffset8x8];
                            B[col+xOffset8x8][row+yOffset8x8]=temp;
                        }
                    }
            	}

            	//switch lower left quadrant, B5&6 with 7&8
		tempXOff = xOffset8x8+2;
            	for(row=rowBlock+4; row<(rowBlock+6); row++){
                    for(col=colBlock; col<(colBlock+4); col++){
                    	temp=B[row+xOffset8x8][col+yOffset8x8];
                        B[row+xOffset8x8][col+yOffset8x8]=B[row+tempXOff][col+yOffset8x8];
                        B[row+tempXOff][col+yOffset8x8]=temp;
                    }
            	}

		//swap upper part of lower left with bottom part of upper right
		tempXOff = xOffset8x8-2;
		tempYOff = yOffset8x8+4;
            	for(row=rowBlock+4; row<(rowBlock+6); row++){
                    for(col=colBlock; col<(colBlock+4); col++){
                    	temp=B[row+xOffset8x8][col+yOffset8x8];
                    	B[row+xOffset8x8][col+yOffset8x8]=B[row+tempXOff][col+tempYOff];
                    	B[row+tempXOff][col+tempYOff]=temp;
                    }
		}

            	//swap lower part of lower left with upper part of upper right
		tempXOff = xOffset8x8-6;
		tempYOff = yOffset8x8+4;
            	for(row=rowBlock+6; row<(rowBlock+8); row++){
                    for(col=colBlock+0; col<(colBlock+4); col++){
                    	temp=B[row+xOffset8x8][col+yOffset8x8];
                    	B[row+xOffset8x8][col+yOffset8x8]=B[row+tempXOff][col+tempYOff];
                    	B[row+tempXOff][col+tempYOff]=temp;
                    }
            	}
	    }
	}
    } else {
	int hasDiag = 0;
	//treat it as 16x4 blocks (up to array bounds)
	for(rowBlock = 0; rowBlock < N; rowBlock += 16){
	    for(colBlock = 0; colBlock < M; colBlock += 4){
		for(row = rowBlock; row < rowBlock+16 && row < N; row++){
		    for(col = colBlock; col < colBlock+4 && col < M; col++){
			if(row != col){
			    B[col][row] = A[row][col];
			} else {
			    //if on a diagonal, signal it hasDiag
			    temp = A[row][row];
			    diagIdx = row;
			    hasDiag = 1;
			}
		    }
		    if(hasDiag){
		        B[diagIdx][diagIdx] = temp;
			hasDiag = 0;
		    }
		}
	    }
	}
    }
    
}

char transpose_desc32[] = "FOR 32x32";
void transpose_32(int M, int N, int A[N][M], int B[M][N])
{
    int row, col, rowBlock, colBlock, temp, diagIdx, stride;
	//treat 32x32 in 8x8 blocks
	stride = 8;
	for(rowBlock = 0; rowBlock < N; rowBlock += stride){
	    for(colBlock = 0; colBlock < M; colBlock += stride){
		for(row = rowBlock; row < rowBlock+stride; row++){
		    for(col = colBlock; col < colBlock+stride; col++){
			if(row != col){
			    B[col][row] = A[row][col];
			} else {
			    //if we have a diagonal we do not need to immediately transpose it
			    temp = A[row][row];
			    diagIdx = row;
			}
		    }
		    if(rowBlock == colBlock){
		        B[diagIdx][diagIdx] = temp;
		    }
		}
	    }
	}
}

char transpose_desc64[] = "FOR 64x64";
void transpose_64(int M, int N, int A[N][M], int B[M][N])
{
    int row, col, rowBlock, colBlock, stride, temp;
    int xOffset8x8, yOffset8x8, tempXOff, tempYOff;
    stride = 8;
	//treating this as 8x8 blocks, but splitting those into 4x4
	for(rowBlock = 0; rowBlock < N; rowBlock += stride){
	    for(colBlock = 0; colBlock < M; colBlock += stride){
		xOffset8x8 = colBlock-rowBlock;
		yOffset8x8 = rowBlock-colBlock;
		
		//Put first and second rows of A into third and fourth rows of B (avoid 1->1 and 2->2 evictions)
		tempXOff = xOffset8x8+2;
		for(row=rowBlock; row<(rowBlock+2); row++){
                    for(col=colBlock; col<(colBlock+8); col++){
                    	B[row+tempXOff][col+yOffset8x8]=A[row][col];
                    }
            	}

		//Put second and fourth rows of A into first and second rows of B (avoid 3->3 and 4->4 evictions)
		tempXOff = xOffset8x8 - 2;		
		for(row=rowBlock+2; row<(rowBlock+4); row++){
                    for(col=colBlock; col<(colBlock+8); col++){
                    	B[row+tempXOff][col+yOffset8x8]=A[row][col];
                    }
            	}

		//1-to-1 translation of A onto B (for given 8x8) (top left quad)
		tempXOff = xOffset8x8 + 2;		
		for(row=rowBlock; row<(rowBlock+2); row++){
                    for(col=colBlock; col<(colBlock+8); col++){
			temp = B[row+xOffset8x8][col+yOffset8x8];
                    	B[row+xOffset8x8][col+yOffset8x8]=B[row+tempXOff][col+yOffset8x8];
			B[row+tempXOff][col+yOffset8x8] = temp;
                    }
            	}
		
		//transpose upper left quad (into final position)
		for(row=rowBlock+xOffset8x8; row<(rowBlock+xOffset8x8+4); row++){
                    for(col=colBlock+yOffset8x8; col<(colBlock+yOffset8x8+4); col++){
                    	if((col-(colBlock+yOffset8x8))<(row-(rowBlock+xOffset8x8))){
                            temp=B[row][col];
                            B[row][col]=B[col+xOffset8x8][row+yOffset8x8];
                            B[col+xOffset8x8][row+yOffset8x8]=temp;
                        }
                    }
            	}

		//transpose upper right quad (right of prev quad final position)
		tempXOff = xOffset8x8-4;
		tempYOff = yOffset8x8+4;
		for(row=rowBlock+xOffset8x8; row<(rowBlock+xOffset8x8+4); row++){
                    for(col=colBlock+yOffset8x8+4; col<(colBlock+yOffset8x8+8); col++){
                    	if((col-(colBlock+yOffset8x8+4))<(row-(rowBlock+xOffset8x8))){
                            temp=B[row][col];
                            B[row][col]=B[col+tempXOff][row+tempYOff];
                            B[col+tempXOff][row+tempYOff]=temp;
                        }
                    }
            	}

		//Swap upper left's rows 1&2 with its own 3&4 respectively
		tempXOff = xOffset8x8+2;
		for(row=rowBlock; row<(rowBlock+2); row++){
                    for(col=colBlock+4; col<(colBlock+8); col++){
                    	temp=B[row+xOffset8x8][col+yOffset8x8];
                    	B[row+xOffset8x8][col+yOffset8x8]=B[row+tempXOff][col+yOffset8x8];
                    	B[row+tempXOff][col+yOffset8x8]=temp;
                    }
            	}
		
		//put fifth and sixth rows of A into 7 and 8 in B (avoids evictions
		tempXOff = xOffset8x8+2;
            	for(row=rowBlock+4; row<(rowBlock+6); row++){
                    for(col=colBlock; col<(colBlock+8); col++){
                    	B[row+tempXOff][col+yOffset8x8]=A[row][col];
                    }
		}

		//put A7 and 8 into B 5 and 6 (avoids evictions)
		tempXOff = xOffset8x8-2;
            	for(row=rowBlock+6; row<(rowBlock+8); row++){
                    for(col=colBlock; col<(colBlock+8); col++){
                    	B[row+tempXOff][col+yOffset8x8]=A[row][col];
                    }
            	}

            	//swap B 5 and 6 with 7 and 8 respectively
		tempXOff = xOffset8x8+2;
            	for(row=rowBlock+4; row<(rowBlock+6); row++){
                    for(col=colBlock; col<(colBlock+8); col++){
                    	temp=B[row+xOffset8x8][col+yOffset8x8];
                    	B[row+xOffset8x8][col+yOffset8x8]=B[row+tempXOff][col+yOffset8x8];
                    	B[row+tempXOff][col+yOffset8x8]=temp;
                    }
            	}

            	//transpose left quadrant of lower half
		tempXOff = xOffset8x8+4;
		tempYOff = yOffset8x8-4;
            	for(row=rowBlock+4+xOffset8x8; row<(rowBlock+8+xOffset8x8); row++){
                    for(col=colBlock+yOffset8x8; col<(colBlock+4+yOffset8x8); col++){
                    	if((col-(colBlock+yOffset8x8))<(row-(rowBlock+4+xOffset8x8))){
                            temp=B[row][col];
                            B[row][col]=B[col+tempXOff][row+tempYOff];
                            B[col+tempXOff][row+tempYOff]=temp;
                    	}
                    }
            	}

            	//transpose right quadrant of lower half
            	for(row=rowBlock+4+xOffset8x8; row<(rowBlock+8+xOffset8x8); row++){
                    for(col=colBlock+4+yOffset8x8; col<(colBlock+8+yOffset8x8); col++){
                    	if((col-(colBlock+4+yOffset8x8))<(row-(rowBlock+4+xOffset8x8))){
                            temp=B[row][col];
                            B[row][col]=B[col+xOffset8x8][row+yOffset8x8];
                            B[col+xOffset8x8][row+yOffset8x8]=temp;
                        }
                    }
            	}

            	//switch lower left quadrant, B5&6 with 7&8
		tempXOff = xOffset8x8+2;
            	for(row=rowBlock+4; row<(rowBlock+6); row++){
                    for(col=colBlock; col<(colBlock+4); col++){
                    	temp=B[row+xOffset8x8][col+yOffset8x8];
                        B[row+xOffset8x8][col+yOffset8x8]=B[row+tempXOff][col+yOffset8x8];
                        B[row+tempXOff][col+yOffset8x8]=temp;
                    }
            	}

		//swap upper part of lower left with bottom part of upper right
		tempXOff = xOffset8x8-2;
		tempYOff = yOffset8x8+4;
            	for(row=rowBlock+4; row<(rowBlock+6); row++){
                    for(col=colBlock; col<(colBlock+4); col++){
                    	temp=B[row+xOffset8x8][col+yOffset8x8];
                    	B[row+xOffset8x8][col+yOffset8x8]=B[row+tempXOff][col+tempYOff];
                    	B[row+tempXOff][col+tempYOff]=temp;
                    }
		}

            	//swap lower part of lower left with upper part of upper right
		tempXOff = xOffset8x8-6;
		tempYOff = yOffset8x8+4;
            	for(row=rowBlock+6; row<(rowBlock+8); row++){
                    for(col=colBlock+0; col<(colBlock+4); col++){
                    	temp=B[row+xOffset8x8][col+yOffset8x8];
                    	B[row+xOffset8x8][col+yOffset8x8]=B[row+tempXOff][col+tempYOff];
                    	B[row+tempXOff][col+tempYOff]=temp;
                    }
            	}
	    }
	}   
}
    
char transpose_descPrime[] = "FOR PRIME NUMBERS";
void transpose_prime(int M, int N, int A[N][M], int B[M][N]){
    int row, col, rowBlock, colBlock, temp, diagIdx;
    int hasDiag = 0;
	//treat it as 16x4 blocks (up to array bounds)
	for(rowBlock = 0; rowBlock < N; rowBlock += 16){
	    for(colBlock = 0; colBlock < M; colBlock += 4){
		for(row = rowBlock; row < rowBlock+16 && row < N; row++){
		    for(col = colBlock; col < colBlock+4 && col < M; col++){
			if(row != col){
			    B[col][row] = A[row][col];
			} else {
			    //if on a diagonal, signal it hasDiag
			    temp = A[row][row];
			    diagIdx = row;
			    hasDiag = 1;
			}
		    }
		    if(hasDiag){
		        B[diagIdx][diagIdx] = temp;
			hasDiag = 0;
		    }
		}
	    }
	}
}

/*char transpose_desc2[] = "USING 2";
void transpose_2(int M, int N, int A[N][M], int B[M][N])
{

}

char transpose_desc[] = "USING 32";
void transpose_6(int M, int N, int A[N][M], int B[M][N])
{
    int i, j, a, b, aOff, bOff;
    for(a = 0; a < N/32; a++){
	aOff = a*32;
	for(b = 0; b < M/32; b++){
	    bOff = b*32;
    	    for(i = 0; i < 32; i++){
		for(j = 0; j < 32; j++){
	 	    B[j+bOff][i+aOff] = A[i+aOff][j+bOff];
		}
   	    }
	}
    }
    
}

char transpose_desc16[] = "USING 16";
void transpose_16(int M, int N, int A[N][M], int B[M][N])
{
    int i, j, a, b, aOff, bOff;
    for(a = 0; a < N/16; a++){
	aOff = a*16;
	for(b = 0; b < M/16; b++){
	    bOff = b*16;
    	    for(i = 0; i < 16; i++){
		for(j = 0; j < 16; j++){
	 	    B[j+bOff][i+aOff] = A[i+aOff][j+bOff];
		}
   	    }
	}
    }
    
}

char transpose_desc32x8[] = "USING 32x8";
void transpose_32x8(int M, int N, int A[N][M], int B[M][N])
{
    int i, j, a, b, aOff, bOff;
    for(a = 0; a < N/32; a++){
	aOff = a*32;
	for(b = 0; b < M/8; b++){
	    bOff = b*8;
    	    for(i = 0; i < 32; i++){
		for(j = 0; j < 8; j++){
	 	    B[j+bOff][i+aOff] = A[i+aOff][j+bOff];
		}
   	    }
	}
    }
    
}

char transpose_desc32x16[] = "USING 32x16";
void transpose_32x16(int M, int N, int A[N][M], int B[M][N])
{
    int i, j, a, b, aOff, bOff;
    for(a = 0; a < N/32; a++){
	aOff = a*32;
	for(b = 0; b < M/16; b++){
	    bOff = b*16;
    	    for(i = 0; i < 32; i++){
		for(j = 0; j < 16; j++){
	 	    B[j+bOff][i+aOff] = A[i+aOff][j+bOff];
		}
   	    }
	}
    }
    
}

char transpose_desc32x4[] = "USING 32x4";
void transpose_32x4(int M, int N, int A[N][M], int B[M][N])
{
    int i, j, a, b, aOff, bOff;
    for(a = 0; a < N/32; a++){
	aOff = a*32;
	for(b = 0; b < M/4; b++){
	    bOff = b*4;
    	    for(i = 0; i < 32 && i+aOff < N; i++){
		for(j = 0; j < 4 && j+bOff < M; j++){
	 	    B[j+bOff][i+aOff] = A[i+aOff][j+bOff];
		}
   	    }
	}
    }
    
}

char transpose_desc8x32[] = "USING 8x32";
void transpose_8x32(int M, int N, int A[N][M], int B[M][N])
{
    int i, j, a, b, aOff, bOff;
    for(a = 0; a < N/8; a++){
	aOff = a*8;
	for(b = 0; b < M/32; b++){
	    bOff = b*32;
    	    for(i = 0; i < 8; i++){
		for(j = 0; j < 32; j++){
	 	    B[j+bOff][i+aOff] = A[i+aOff][j+bOff];
		}
   	    }
	}
    }
    
}

char transpose_desc8x4[] = "USING 8x4";
void transpose_8x4(int M, int N, int A[N][M], int B[M][N])
{
    int i, j, a, b, aOff, bOff;
    for(a = 0; a < N/8; a++){
	aOff = a*8;
	for(b = 0; b < M/4; b++){
	    bOff = b*4;
    	    for(i = 0; i < 8; i++){
		for(j = 0; j < 4; j++){
	 	    B[j+bOff][i+aOff] = A[i+aOff][j+bOff];
		}
   	    }
	}
    }
    
}

char transpose_desc8x8[] = "USING 8x8";
void transpose_8x8(int M, int N, int A[N][M], int B[M][N])
{
    int i, j, a, b, aOff, bOff;
    for(a = 0; a < N/8; a++){
	aOff = a*8;
	for(b = 0; b < M/8; b++){
	    bOff = b*8;
    	    for(i = 0; i < 8; i++){
		for(j = 0; j < 8; j++){
	 	    B[j+bOff][i+aOff] = A[i+aOff][j+bOff];
		}
   	    }
	}
    }
    
}

char transpose_desc8x16[] = "USING 8x16";
void transpose_8x16(int M, int N, int A[N][M], int B[M][N])
{
    int i, j, a, b, aOff, bOff;
    for(a = 0; a < N/8; a++){
	aOff = a*8;
	for(b = 0; b < M/16; b++){
	    bOff = b*16;
    	    for(i = 0; i < 8; i++){
		for(j = 0; j < 16; j++){
	 	    B[j+bOff][i+aOff] = A[i+aOff][j+bOff];
		}
   	    }
	}
    }
    
}*/
/* 
 * You can define additional transpose functions below. We've defined
 * a simple one below to help you get started. 
 */ 

/* 
 * trans - A simple baseline transpose function, not optimized for the cache.
 */
/*
char trans_desc[] = "Simple row-wise scan transpose";
void trans(int M, int N, int A[N][M], int B[M][N])
{
    int i, j, tmp;

    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            tmp = A[i][j];
            B[j][i] = tmp;
        }
    }    

}
*/
/*
 * registerFunctions - This function registers your transpose
 *     functions with the driver.  At runtime, the driver will
 *     evaluate each of the registered functions and summarize their
 *     performance. This is a handy way to experiment with different
 *     transpose strategies.
 */
void registerFunctions()
{
    /* Register your solution function */
    registerTransFunction(transpose_submit, transpose_submit_desc); 

    /* Register any additional transpose functions */
    //registerTransFunction(trans, trans_desc); 
    registerTransFunction(transpose_prime, transpose_descPrime); 
    //registerTransFunction(transpose_16, transpose_desc16);
    registerTransFunction(transpose_32, transpose_desc32);  
    registerTransFunction(transpose_64, transpose_desc64); 
    /*registerTransFunction(transpose_8x4, transpose_desc8x4);
    registerTransFunction(transpose_8x8, transpose_desc8x8);
    registerTransFunction(transpose_8x16, transpose_desc8x16);
    registerTransFunction(transpose_8x32, transpose_desc8x32);
    registerTransFunction(transpose_32x4, transpose_desc32x4); 
    registerTransFunction(transpose_32x8, transpose_desc32x8);
    registerTransFunction(transpose_32x16, transpose_desc32x16); 
    registerTransFunction(transpose_2, transpose_desc2); */

}

/* 
 * is_transpose - This helper function checks if B is the transpose of
 *     A. You can check the correctness of your transpose by calling
 *     it before returning from the transpose function.
 */
int is_transpose(int M, int N, int A[N][M], int B[M][N])
{
    int i, j;

    for (i = 0; i < N; i++) {
        for (j = 0; j < M; ++j) {
            if (A[i][j] != B[j][i]) {
                return 0;
            }
        }
    }
    return 1;
}

