
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <time.h>
#include <conio.h>
#define MAX 6 // 48 * 32 integers - 1536 bit max
#define THREAD_NUMBER 5
#define INT_SIZE 32

unsigned int	    
		
		ZERO_MEMORY[MAX],
		
		KKK[33] = { 3, 5, 17, 23, 27, 29, 39, 47, 57, 65, 71, 75, 77, 83, 93, 99, 105, 107, 113, 117, 129, 135, 143, 149, 153, 159, 167, 173, 185, 189, 195, 199, 203 },
		HHH[33] = { 285, 297, 299, 303, 309, 315, 323, 327, 329, 339, 353, 359, 363, 365, 369, 383, 387, 395, 413, 419, 429, 437, 453, 465, 467, 479, 483, 485, 489, 497, 507, 509, 513 },

		Base[THREAD_NUMBER],
		Base_[THREAD_NUMBER],

		N_Me,					
		Me,						
		InverseByModuleM_Me,	
		InverseByModuleMMMe,	
		input[THREAD_NUMBER],	
		input_[THREAD_NUMBER],
		N[THREAD_NUMBER] ,
		N_[THREAD_NUMBER] ,
		Mi[THREAD_NUMBER][MAX],     // Large
		Mi_[THREAD_NUMBER][MAX],    // Large
		M[ MAX ],					//
		M_[ MAX ],					//
		$N[MAX],
		MiInv[THREAD_NUMBER] ,
		Mi_Inv[THREAD_NUMBER] ,
		SQR_M_MOD_N_RNS[THREAD_NUMBER],
		SQR_M_MOD_N_RNS_[THREAD_NUMBER] ,
		InverseByModuleMM_[THREAD_NUMBER] ,
		InverseByModuleNM[THREAD_NUMBER] ,
		Mi_InMe[THREAD_NUMBER] ,
		MiInMe [THREAD_NUMBER],
		M_InBase[THREAD_NUMBER],
		MiInBase_[THREAD_NUMBER][THREAD_NUMBER],
		Mi_InBase[THREAD_NUMBER][THREAD_NUMBER],
		ABmodN[2][THREAD_NUMBER];

	bool power[ MAX * INT_SIZE ];
	bool e[ MAX * INT_SIZE ];
	bool d[ MAX * INT_SIZE ];
	
	
	int iterationsCounter,
		lastValueOfStateInformation,
		currentOperationIndex,
		eBitsCount, 
		numberOfIterationsForE,
		dBitsCount, 
		numberOfIterationsForD;

int profiler_Inversions = 0, 
	profiler_extended_euclid_Iterations[1000000], 
	profiler_eeI_counter = -1;


void clearLongVariable ( unsigned int *a ) {
	
	memcpy( a, ZERO_MEMORY, MAX * sizeof(int) );

}

void copyVariable ( unsigned int *source, unsigned int *dest ) {
	
	memcpy( dest, source, MAX * sizeof(int) );

}

// зсунути 32-бітне число вліво на 1 біт і повернути carry flag
int shl (unsigned int &input, int rightBit){
	
	int CF;
	CF = ( input & ( 1 << INT_SIZE - 1 ) ) == ( 1 << INT_SIZE - 1 );
	input <<= 1;             
	input &= 0xFFFFFFFE;	
	input |= rightBit;		
	return CF;
}

// зсунути 32-бітне число вправо на 1 біт і повернути carry flag
int shr (unsigned int &input, int leftBit){
	
	int CF;
	CF =  input & 0x1 ;  
	input >>= 1;
		
	if ( leftBit == 1 ) { 
		input |= 0x80000000;	
	} else { 
		input &= 0x7FFFFFFF;	
	}
	
	return CF;
}

// зсунути 64-бітне число вправо на 1 біт і повернути carry flag
int shr_long (unsigned long long &input, int leftBit){
	
	int CF;
	CF =  input & 0x1 ;  
	input >>= 1;
		
	if ( leftBit == 1 ) { 
		input |= 0x8000000000000000;	
	} else { 
		input &= 0x7FFFFFFFFFFFFFFF;	
	}
	
	return CF;
}

//зсунути "довге число" вліво на 1 біт і повернути carry flag
int shiftToLeftVariable ( unsigned int *input ) {
	int CF = 0;
	for (int i = MAX - 1; i >= 0; i--) {
		CF = shl( input[ i ], CF );
	}
	return CF;
}

//зсунути "довге число" вправо на 1 біт і повернути carry flag
int shiftToRightVariable ( unsigned int *input ) {
	int CF = 0;
	for (int i = 0; i < MAX; i++) {
		CF = shr( input[ i ], CF );
	}

	return CF;
}

void convertToNormalForm (char input[], unsigned int *a){
	
	clearLongVariable (a);

	int temp = 0;
    for( int i = MAX - 1, 
		 int j = strlen(input) - 1, 
		 int p = 0; 
					j >= 0; j-- ) {
		
			if ( input[j] == 32 ) { continue; } // дозволено пробіл для зручності
			else if ( toupper(input[j]) == 65 ) { temp |= 0xA << p; }
			else if ( toupper(input[j]) == 66 ) { temp |= 0xB << p; }
			else if ( toupper(input[j]) == 67 ) { temp |= 0xC << p; }
			else if ( toupper(input[j]) == 68 ) { temp |= 0xD << p; }
			else if ( toupper(input[j]) == 69 ) { temp |= 0xE << p; }
			else if ( toupper(input[j]) == 70 ) { temp |= 0xF << p; }
			else { temp |= ( toupper(input[j]) - 48 ) << p; }
			 
	    p+=4;
		if( p >= 32 || j == 0){ p = 0; a[i] = temp; temp = 0; i--; }
	}

}

void convertFromNormalForm (char output[], unsigned int *a) {
	
	int temp;
	int j = ( MAX * 8 ) - 1;
	output [ j + 1 ] = '\0';
    for( int i = MAX - 1; i >= 0; i-- ) {
		
		int mask = 0xF;

		for ( int p = 0; p < 8; p++) {
			temp = a[ i ] ;
			temp >>= 4*p;
			temp &=  mask ;
			output [ j ] = ( temp > 9 ) ? ( temp + 55 ) : ( temp + 48 );
			j--;
		}
	}
}

void convertPowerToArrayOfBits( unsigned int* $power, bool *power, int *bin_digits_in_power, int *totalIterationsNeeded ){
	
	int c = 0;

	for(int i = MAX * INT_SIZE - 1, bit; i > 0 ; i--){

		bit = shiftToRightVariable($power);
		if (bit == 1){
			*bin_digits_in_power = MAX * INT_SIZE - i;
			c++;
		}
		power[i] = (bit == 1) ? true : false;
	}
	c *= 2;
	c += 2 * (*bin_digits_in_power);


	*totalIterationsNeeded = c;

}




void add ( unsigned int *a,  
		   unsigned int *b,  
		   unsigned int *c) 
{
	int CF = 0;
	
	for (int i = MAX - 1; i >= 0; i--) {

		c[ i ] = a[ i ] + b[ i ] + CF;
		
		if ( ( c[ i ] < a [ i ] ) || ( c[ i ] < b [ i ] ) ) { 
			
			CF = 1; 

		} else {

			CF = 0;
		}
	}
}


void sub ( unsigned int *a,  
		   unsigned int *b,  
		   unsigned int *c) 
{

		
	unsigned int tempResult [MAX]; 
	clearLongVariable ( tempResult );
	
	int CF = 0;

	for (int i = MAX - 1; i >= 0; i--) {
		tempResult[i] = a[i] - b[i] - CF;

		if(b[i] == 0xFFFFFFFF){ 
			CF = 1; 
			continue;
		}
		if(a[i] < (b[i] + CF) ){
			CF = 1;
		} else {
			CF=0;
		}

	}
	copyVariable ( tempResult, c );

}



//
//  1 a > b
//  0 a = b
// -1 a < b
// знакове порівняння чисел
int cmp (unsigned int *a, unsigned int *b) 
{
	
	bool a_positive = !(a[0] & 0x80000000);
	bool b_positive = !(b[0] & 0x80000000);

	// 1. + +
	// 1. + -
	// 1. - +
	// 1. - -

	if ( a_positive && b_positive ){

		for (int i = 0; i < MAX; i++) {
			if ( a[ i ] > b[ i ] ) { 
				return 1; 
			}

			if ( a[ i ] < b[ i ] ) { 
				return -1; 
			}

		}
		return 0;

	} else if(a_positive && !b_positive) {

		return 1;

	} else if(!a_positive && b_positive) {

		return -1;

	} else { // two numbers are negative:
	
		for (int i = 0; i < MAX; i++) {

			if ( a[ i ] > b[ i ] ) { 
				return 1; 
			}
			
			if ( a[ i ] < b[ i ] ) { 
				return -1; 
			}

		}
		return 0;
	}
}




void mod (  unsigned int *a,  
			unsigned int *b, 
			unsigned int *c) 
{
	
		unsigned int temp_A[ MAX ];
		unsigned int R[ MAX ];	// remainder

		copyVariable ( a, temp_A );


		if ( cmp ( temp_A, b) <= 0) {
			copyVariable ( temp_A, c );
			return;
		}

		
		clearLongVariable ( R );   // initialize remainder to zero

		// Integer division (unsigned) with remainder
		for (int i = MAX * INT_SIZE - 1; i >= 0 ; i-- ) {
		
			shiftToLeftVariable ( R );					   // left-shift R by 1 bit  
			R[ MAX - 1] |= shiftToLeftVariable ( temp_A ); // set the least-significant bit of R equal to bit i of the numerator

			if ( cmp( R, b ) >= 0 ) {
				sub ( R, b, R );	
			}
		}
		copyVariable ( R, c );   
}




void mod_ ( unsigned int *a,  
			unsigned int b,  
			unsigned int* c) 
{

		unsigned int temp_A[ MAX ];
		unsigned long long R;	// остача

		copyVariable ( a, temp_A );
		
		R = 0;   

		for (int i = MAX * INT_SIZE - 1; i >= 0 ; i-- ) {
		
			R <<= 1;					    
			R |= shiftToLeftVariable ( temp_A ); 

			if ( R >= b ) {
				R -= b;	
			}
				
		}

		*c = (unsigned int) R; 

}


// division
void div (  unsigned int *a,  
			unsigned int *b,  
			unsigned int *c) 
{

		unsigned int temp_A[ MAX ];
		unsigned int temp_B[ MAX ];
		unsigned int Q[ MAX ];	// остача 
		unsigned int R[ MAX ];	// залишок 
		unsigned int zero [MAX];
		unsigned int mask = 0x80000000;
		bool aIsLessThenZero = false, 
			 bIsLessThenZero = false;

		for(int i=0; i< MAX;i++){
			zero[i] = 0;
		}

		if (cmp (a,zero) < 0) {
			sub (zero, a, temp_A);
			aIsLessThenZero = true;
		} else {
			copyVariable ( a, temp_A );
		}

		if (cmp (b,zero) < 0) {
			sub (zero, b, temp_B);
			bIsLessThenZero = true;
		} else {
			copyVariable ( b, temp_B );
		}


		clearLongVariable ( Q );   
		clearLongVariable ( R );   
		
		//Integer division (unsigned) with remainder
		//http://en.wikipedia.org/wiki/Division_algorithm
		int j;
		for (int i = MAX * INT_SIZE - 1; i >= 0 ; i-- ) {
		
			j =  i / INT_SIZE ;
			shiftToLeftVariable ( R );					    
			R[ MAX - 1] |= shiftToLeftVariable ( temp_A ); 

			if ( cmp( R, temp_B ) >= 0 ) {
			
				sub ( R, temp_B, R );	
				Q[MAX - 1 - j] |= mask;
			}

			mask >>= 1;
			if(mask == 0) mask = 0x80000000;
				
		}

		// - / - = +
		// + / - = -
		// - / + = -
		// + / + = +
		if (aIsLessThenZero ^ bIsLessThenZero) {
			sub (zero, Q, Q);
		}

		copyVariable ( Q, c );   

}


void div_ ( unsigned int *a,  
			unsigned int b,  
			unsigned int *c) 
{

		unsigned int temp_A[ MAX ];
		unsigned long long temp_B;
		unsigned int Q[ MAX ];	// остача 
		unsigned long long R;	// залишок 
		unsigned int zero [MAX];
		bool aIsLessThenZero = false, 
			 bIsLessThenZero = false;

		for(int i=0; i< MAX;i++){
			zero[i] = 0;
		}

		if (cmp (a,zero) < 0) {
			sub (zero, a, temp_A);
			aIsLessThenZero = true;
		} else {
			copyVariable ( a, temp_A );
		}

		if ( b < 0 ) {
			temp_B = b * (-1);
			bIsLessThenZero = true;
		} else {
			temp_B = b;
		}


		clearLongVariable ( Q );  
		R = 0 ;   
		
		unsigned int mask = 0x80000000;
		
		//Integer division (unsigned) with remainder
		//http://en.wikipedia.org/wiki/Division_algorithm
		int j;
		for (int i = MAX * INT_SIZE - 1; i >= 0 ; i-- ) {
		
			j =  i / INT_SIZE ;
			R <<= 1;					    
			R |= shiftToLeftVariable ( temp_A ); 

			if (  R >= temp_B ) {
				R -= temp_B;	
				Q[MAX - 1 - j] |= mask;
			}

			mask >>= 1;
			if(mask == 0) mask = 0x80000000;
				
		}

		// - / - = +
		// + / - = -
		// - / + = -
		// + / + = +
		if (aIsLessThenZero ^ bIsLessThenZero) {
			sub (zero, Q, Q);
		}

		copyVariable ( Q, c );   

}




void mul (unsigned int *firstFactor,  
		  unsigned int *secondFactor,  
		  unsigned int *product) 
{

	unsigned int tempFirstFactor[ MAX ]; 
	unsigned int tempSecondFactor[ MAX ];
	unsigned int tempProduct[ MAX ];
	unsigned int zero [MAX];

	bool firstFactorIsLessThenZero = false, 
		 secondFactorIsLessThenZero = false;

	for(int i=0; i < MAX;i++){
		zero[i] = 0;
	}

	if (cmp (firstFactor,zero) < 0) {
		sub (zero, firstFactor, tempFirstFactor);
		firstFactorIsLessThenZero = true;
	} else {
		copyVariable ( firstFactor, tempFirstFactor );
	}

	if (cmp (secondFactor,zero) < 0) {
		sub (zero, secondFactor, tempSecondFactor);
		secondFactorIsLessThenZero = true;
	} else {
		copyVariable ( secondFactor, tempSecondFactor );
	}

	int CF;
	clearLongVariable (product);

	for (int i = 0; i < MAX * INT_SIZE; i++ ) {
		
		CF = shiftToRightVariable ( tempSecondFactor );	

		if ( CF == 1 ) {
			add ( tempFirstFactor, product, product);
		}

		shiftToLeftVariable ( tempFirstFactor );	

	}


	// - / - = +
	// + / - = -
	// - / + = -
	// + / + = +

	copyVariable(product, tempProduct);

	if (firstFactorIsLessThenZero ^ secondFactorIsLessThenZero) {
		sub (zero, product, product);
	}


}



// a * x = -1 mod b
// a * x + b * y = gcd(a,b)
void extended_euclid( unsigned int *a,  
					  unsigned int *b,  
					  unsigned int *x,  
					  unsigned int *y,  
					  unsigned int *d )
{
	unsigned int x1[ MAX ];
	unsigned int x2[ MAX ];
	unsigned int y1[ MAX ];
	unsigned int y2[ MAX ];
	unsigned int q[ MAX ];
	unsigned int r[ MAX ];
	unsigned int zero[ MAX ];
	unsigned int temp[ MAX ];
	unsigned int temp_a[ MAX ];
	unsigned int temp_b[ MAX ];
	
	copyVariable( a, temp_a );
	copyVariable( b, temp_b );
	clearLongVariable ( x1 );
	clearLongVariable ( x2 );
	clearLongVariable ( y1 );
	clearLongVariable ( y2 );
	clearLongVariable ( q );
	clearLongVariable ( r );
	clearLongVariable ( zero );
	clearLongVariable ( temp );

	x2[ MAX-1 ] = 1; // x2 = 1
	x1[ MAX-1 ] = 0; // x1 = 0
	y2[ MAX-1 ] = 0; // y2 = 0
	y1[ MAX-1 ] = 1; // y1 = 1

	profiler_eeI_counter++;
	profiler_extended_euclid_Iterations[profiler_eeI_counter] = 0;

	while ( cmp( temp_b, zero ) > 0) {

		profiler_extended_euclid_Iterations[profiler_eeI_counter]++;

		div ( temp_a, temp_b, q);		//q = a / b, 
		mul ( q, temp_b, temp );
		sub ( temp_a, temp, r);			//r = a - q * b;
		clearLongVariable ( temp );
		mul ( q, x1, temp );
		sub ( x2, temp, x );			//*x = x2 - q * x1, 
		clearLongVariable ( temp );
		mul( q, y1, temp );
		sub( y2, temp, y );				//*y = y2 - q * y1;
		copyVariable( temp_b, temp_a );	//a = b, 
		copyVariable( r, temp_b );		//b = r;
		copyVariable( x1, x2 );			//x2 = x1,  
		copyVariable( x, x1 );			//x1 = *x,  // t
		copyVariable( y1, y2 );			//y2 = y1,  //
		copyVariable( y, y1 );			//y1 = *y;
	 
	}

	copyVariable( temp_a, d );	//*d = a,
	copyVariable( x2, x );		//*x = x2,
	copyVariable( y2, y );		//*y = y2;

}


//input * result = 1 mod module
void InverseByModule( unsigned int *input,  
					  unsigned int *module,  
					  unsigned int *result )
{
	unsigned int x[MAX];
	unsigned int y[MAX];
	unsigned int d[MAX];
	unsigned int zero[MAX];
	
	profiler_Inversions++;

	clearLongVariable(x);
	clearLongVariable(y);
	clearLongVariable(d);
	clearLongVariable(zero);

	extended_euclid(input, module, x, y, d);

	if (cmp ( x, zero ) < 0 ) { 
		add ( x, module, zero );
		copyVariable( zero, result );
	} else {
		copyVariable( x, result );
	}
	//return x < 0 ? x + module : x;
}


void InverseByModule_( unsigned int *input,  
					   unsigned int module,  
					   unsigned int *result )
{
	unsigned int x[MAX];
	unsigned int y[MAX];
	unsigned int d[MAX];
	unsigned int temp[MAX];
	unsigned int zero[MAX];
	
	profiler_Inversions++;

	clearLongVariable(x);
	clearLongVariable(y);
	clearLongVariable(d);
	clearLongVariable(temp);
	clearLongVariable(zero);

	temp[ MAX - 1 ] = module;

	extended_euclid(input, temp, x, y, d);

	if (cmp ( x, zero ) < 0 ) { 
		add ( x, temp, zero );
		 *result = zero[MAX-1];
	} else {
		 *result = x[MAX-1];
	}
}



void getNumberInRNSByModMe( unsigned int* input, unsigned long long* result){

	unsigned int X[MAX];
	unsigned int R[MAX];
	unsigned int tempProduct[MAX];
	unsigned int tempM[MAX];
	unsigned int tempFirstFactor[MAX];
	unsigned long long temp;
	

	for (int i = 0; i < MAX; i++){
		tempM[ i ] = M[ i ];
		X[ i ] = 0;
		R[ i ] = 0;
	}

	// X[i] = INPUT[i] * MiInv[i] * Mi[i]
	for (int i = 0; i < THREAD_NUMBER; i++){
		
		for (int j = 0; j < MAX; j++){
			tempProduct[j] = 0;
			tempFirstFactor [j] = Mi[i][j];
		}

		temp = input[i];
		temp *= MiInv[i];

		while( temp !=0 ){

			if ( shr_long (temp, 0) == 1 ) {
				add ( tempFirstFactor, tempProduct, tempProduct);
			}

			shiftToLeftVariable ( tempFirstFactor );	
		}

		add ( X, tempProduct, X);

	}

	//X[i] % M
	for (int i = MAX * INT_SIZE - 1; i >= 0 ; i-- ) {
		
		shiftToLeftVariable ( R );					    
		R[ MAX - 1] |= shiftToLeftVariable ( X );

		if ( cmp( R, tempM ) >= 0 ) {
			sub ( R, tempM, R );	
		}
	}

	//X[i] % M % Me (Me = 2^6)
	*result = 0;
	*result += ( shiftToRightVariable ( R )) == 1 ? 1 : 0;	
	*result += ( shiftToRightVariable ( R )) == 1 ? 2 : 0;	
	*result += ( shiftToRightVariable ( R )) == 1 ? 4 : 0;	
	*result += ( shiftToRightVariable ( R )) == 1 ? 8 : 0;	
	*result += ( shiftToRightVariable ( R )) == 1 ? 16 : 0;
	*result += ( shiftToRightVariable ( R )) == 1 ? 32 : 0;

}

void convertFromRNS( unsigned int* input, unsigned int* result){

	unsigned int X[MAX];
	unsigned int R[MAX];
	unsigned int tempInput[THREAD_NUMBER];
	unsigned int tempProduct[MAX];
	unsigned int tempM[MAX];
	unsigned int tempFirstFactor[MAX];
	unsigned long long temp;
	

	for (int i = 0; i < MAX; i++){
		tempM[ i ] = M[ i ];
		X[ i ] = 0;
		R[ i ] = 0;
	}

	
	// X[i] = INPUT[i] * MiInv[i] * Mi[i]
	for (int i = 0; i < THREAD_NUMBER; i++){
		
		for (int j = 0; j < MAX; j++){
			tempProduct[j] = 0;
			tempFirstFactor [j] = Mi[i][j];
		}

		tempInput[i] = input[i];

		temp = tempInput[i];
		temp *= MiInv[i];

		while( temp !=0 ){

			if ( shr_long (temp, 0) == 1 ) {
				add ( tempFirstFactor, tempProduct, tempProduct);
			}

			shiftToLeftVariable ( tempFirstFactor );	
		}

		add ( X, tempProduct, X);

	}

	//X[i] % M
	for (int i = MAX * INT_SIZE - 1; i >= 0 ; i-- ) {
		
		shiftToLeftVariable ( R );					   
		R[ MAX - 1] |= shiftToLeftVariable ( X ); 

		if ( cmp( R, tempM ) >= 0 ) {
			sub ( R, tempM, R );	
		}
	}
	

	mod(R, $N, result);


}

void convertFromRNS_( unsigned int* input, unsigned int* result){

	unsigned int X[MAX];
	unsigned int R[MAX];
	unsigned int tempProduct[MAX];
	unsigned int tempInput[THREAD_NUMBER];
	unsigned int tempM[MAX];
	unsigned int tempFirstFactor[MAX];
	unsigned long long temp;
	

	for (int i = 0; i < MAX; i++){
		tempM[ i ] = M_[ i ];
		X[ i ] = 0;
		R[ i ] = 0;
	}

	

	// X[i] = INPUT[i] * MiInv[i] * Mi[i]
	for (int i = 0; i < THREAD_NUMBER; i++){
		
		for (int j = 0; j < MAX; j++){
			tempProduct[j] = 0;
			tempFirstFactor [j] = Mi_[i][j];
		}

		tempInput[i] = input[i];
		
		temp = tempInput[i];
		temp *= Mi_Inv[i];

		while( temp !=0 ){

			if ( shr_long (temp, 0) == 1 ) {
				add ( tempFirstFactor, tempProduct, tempProduct);
			}

			shiftToLeftVariable ( tempFirstFactor );	
		}

		add ( X, tempProduct, X);

	}

	//X[i] % M
	for (int i = MAX * INT_SIZE - 1; i >= 0 ; i-- ) {
		
		shiftToLeftVariable ( R );					  
		R[ MAX - 1] |= shiftToLeftVariable ( X ); 

		if ( cmp( R, tempM ) >= 0 ) {
			sub ( R, tempM, R );	
		}
	}
	

	mod(R, $N, result);


}

void printCurrentStateInformation(){

	
	//system("cls");
	if( currentOperationIndex == 0 ){
		if( lastValueOfStateInformation != iterationsCounter * 100 /  numberOfIterationsForE ){
			lastValueOfStateInformation = iterationsCounter * 100 /  numberOfIterationsForE;
			printf("Encryption... %i %% done.\n", iterationsCounter * 100 /  numberOfIterationsForE );	
		}

	} else if( currentOperationIndex == 1 ){
		if( lastValueOfStateInformation != iterationsCounter * 100 / numberOfIterationsForD ){
			lastValueOfStateInformation = iterationsCounter * 100 / numberOfIterationsForD;
			printf("Decryption... %i %% done.\n", iterationsCounter * 100 /  numberOfIterationsForD );	
		}

	}

	iterationsCounter++ ;

}




void MM( unsigned int* A, 
		 unsigned int* A_, 
		 unsigned int* B, 
		 unsigned int* B_,
		 unsigned int* R, 
		 unsigned int* R_ )
{

	unsigned long long Q[THREAD_NUMBER];
	unsigned long long Xi[THREAD_NUMBER];
	unsigned long long Q_Me, R_Me, A_Me, B_Me, A_Me2, B_Me2 ;
	unsigned long long Q_[THREAD_NUMBER]; 
	unsigned long long s, s1, s2, currSum, r; 
	unsigned long long Sig[THREAD_NUMBER];
	unsigned long long Beta, temp;


	unsigned long long tempA[THREAD_NUMBER];
	unsigned long long tempB[THREAD_NUMBER];

	int rightBit;

	clock_t begin;
	begin = clock();


	//computation of Q
	for (int i=0; i<THREAD_NUMBER; i++){

		Q[i] = ( Base[ i ]  - A[ i ] ) % Base[i];
		Q[i] *= (B[ i ] % Base[i]);
		Q[i] %= Base[i];
		Q[i] *= (InverseByModuleNM[i] % Base[i]);
		Q[i] %= Base[i] ;
	}

	

	// First Base extension:
	//****************************************************************************

	for (int i=0; i<THREAD_NUMBER; i++){
		Sig [i] = Q[i] % Base [i];
		Sig [i] *= MiInv[i] % Base [i];
		Sig [i] %= Base [i];
	}

	

	for (int i=0; i<THREAD_NUMBER; i++){
		s=0;
		for (int j=0; j<THREAD_NUMBER; j++){
			s += ((MiInBase_[i][j] % Base_ [i]) * (Sig[j] % Base_ [i])) % Base_ [i];
		}


		Q_[i] = s % Base_ [i];
	}

	

	//Extra modulus computation;
	// Only first thread:
	s=0;

	for (int i=0; i<THREAD_NUMBER; i++){
		s += ((MiInMe[i] % Me) * ((Q[i] * MiInv[i]) %  Base [i]) % Me ) % Me;

	}

	Q_Me = s % Me;


	for (int i=0; i<THREAD_NUMBER; i++){
		tempA[i] = A[i];
		tempB[i] = B[i];
	}





	unsigned long long jjj = 1;
	unsigned long long j;
	unsigned int ttempA, ttempB;

	for (int i=0; i<THREAD_NUMBER; i++){

		InverseByModule_(Mi[i], 64, &ttempA);

		

		

		jjj *= ttempA;
	}
		
	
	A_Me2 = jjj % Me;

	

	getNumberInRNSByModMe(A, &A_Me);
	getNumberInRNSByModMe(B, &B_Me);

	//printf("\n[%4.2f] Q_Me, A_Me, B_Me\n",(double)(clock() - begin) / CLOCKS_PER_SEC);

	// Computing R in Base_:

	for (int i=0; i<THREAD_NUMBER; i++){ 
		temp = ( (Q_[i] * N_[i]) % Base_[i] + ((unsigned long long)A_[i] * (unsigned long long)B_[i]) % Base_[i] ) % Base_[i] ;
		temp *= InverseByModuleMM_[i] % Base_[i];
		temp %=  Base_[i];
		R_[i] = (unsigned int)temp;

	}

	R_Me = (((A_Me * B_Me)  + (Q_Me * N_Me) )  * InverseByModuleMMMe ) % Me ;

	//printf("\n[%4.2f] R_, R_Me\n",(double)(clock() - begin) / CLOCKS_PER_SEC);

	//Second base extension:
	//*********************************************************************************************

	for (int i=0; i<THREAD_NUMBER; i++){

		Xi [i] = R_[i]  % Base_ [i];
		Xi [i] *= Mi_Inv[i]  % Base_ [i];
		Xi [i] %= Base_ [i];
	}


	//Computing beta:
	s=0;
	for (int j=0; j<THREAD_NUMBER; j++){
		s +=   Xi [j] * Mi_InMe[j];
	}

	Beta = ( ((s - R_Me) % Me) * InverseByModuleM_Me) % Me;

	//printf("\n[%4.2f] Beta\n",(double)(clock() - begin) / CLOCKS_PER_SEC);

	for (int i=0; i<THREAD_NUMBER; i++){
		s1=0; 
		s2=0; 
		currSum=0; //[s1, s2] 128 bit max;

		for (int j=0; j<THREAD_NUMBER; j++){

			s2 += Mi_InBase[i][j] * Xi[j];

			//overflow
			if(s2 < currSum ) {
				s1++;
			} 

			currSum = s2;

		}


		if ( s1==0 ){

			R[i] = (  s2 - (( Beta % Base[i] ) * M_InBase[i] ) % Base[i]    ) % Base[i];

		} else {

			r = 0;

			s2 -= (( Beta % Base[i] ) * M_InBase[i] ) % Base[i] ;

			for (int j = 4 * INT_SIZE - 1; j > 0 ; j-- ) {

				rightBit = ( s2 & ( 0x8000000000000000 ) ) == ( 0x8000000000000000 ) ? 1 : 0; 
				r <<= 1;
				s2 <<=1;
				s1 <<= 1;            
				s1 &= 0xFFFFFFFFFFFFFFFE;	
				s1 |= rightBit;		
				r |= ( s1 & ( 0x8000000000000000 ) ) == ( 0x8000000000000000 ) ? 1 : 0; 

				if ( r > Base[i] ) {
					r -= Base[i];	
				}
			}

			R[i] = r;

		}

	}

	//printf("\n[%4.2f] R\n",(double)(clock() - begin) / CLOCKS_PER_SEC);

	//printCurrentStateInformation();

	//(A * B * M^-1 mod N + Beta * N) mod M
}



void showProfilerInfo(){
printf("\nInversions:			%i", profiler_Inversions);
}



void generateStaticData(){

	clock_t begin;
	unsigned int temp[ MAX ];
	unsigned int $SQR_M_MOD_N[ MAX ];
	unsigned int $InverseByModuleMM_[ MAX ];
	unsigned int $InverseByModule$NM[ MAX ];

	begin = clock();
	// M = Base[0] * Base[1] * ... * Base[n]
	// M_ = Base_[0] * Base_[1] * ... * Base_[n]
	convertToNormalForm( "1" , M); 
	convertToNormalForm( "1" , M_); 

	for(int i=0; i < THREAD_NUMBER; i++){
		clearLongVariable( temp );
		temp[ MAX - 1 ] = Base[i];
		mul ( M, temp, M );

		temp[ MAX - 1 ] = Base_[i];
		mul ( M_, temp, M_ );

	}

	printf("\n[%4.2f] M and M_ done.\n",(double)(clock() - begin) / CLOCKS_PER_SEC);

	clearLongVariable( $SQR_M_MOD_N );
	clearLongVariable( $InverseByModuleMM_ );
	clearLongVariable( $InverseByModule$NM );
	clearLongVariable( temp );
	
	// $SQR_M_MOD_N = M*M % $N;
	//mul ( M, M, temp );
	//mod (temp,  $N,  $SQR_M_MOD_N); 
	 
	
	//sqrMmodN(temp, $N, $SQR_M_MOD_N);

	mod (M,  $N,  temp);
	mul (temp,  temp,  $SQR_M_MOD_N); 
	copyVariable($SQR_M_MOD_N,temp);
	mod (temp,$N,$SQR_M_MOD_N);

	

	printf("\n[%4.2f] $SQR_M_MOD_N done.\n",(double)(clock() - begin) / CLOCKS_PER_SEC);

	InverseByModule(M, M_, $InverseByModuleMM_);
	printf("\n[%4.2f] $InverseByModuleMM_ done [%i iterations].\n",(double)(clock() - begin) / CLOCKS_PER_SEC, profiler_extended_euclid_Iterations[0]);
	
	InverseByModule ($N, M , $InverseByModule$NM);
	printf("\n[%4.2f] $InverseByModule$NM done [%i iterations].\n",(double)(clock() - begin) / CLOCKS_PER_SEC, profiler_extended_euclid_Iterations[1]);

	for( int i = 0; i < THREAD_NUMBER; i++){
			 
		ABmodN[0][i] = 1;
		ABmodN[1][i] = 1;
		
		mod_($N, Base[i], &N[i]);
		mod_($N, Base_[i], &N_[i]);
		
		div_(M, Base[i], Mi[i]);
		div_(M_, Base_[i], Mi_[i]);
		
		InverseByModule_(Mi[i], Base[i], &MiInv[i]);
		InverseByModule_(Mi_[i], Base_[i], &Mi_Inv[i]);
		
		mod_($SQR_M_MOD_N, Base[i], &SQR_M_MOD_N_RNS[i]);
		mod_($SQR_M_MOD_N, Base_[i], &SQR_M_MOD_N_RNS_[i]);

		mod_($InverseByModuleMM_, Base_[i], &InverseByModuleMM_[i]);
		mod_($InverseByModule$NM, Base[i], &InverseByModuleNM[i]);

		mod_(Mi_[i], Me, &Mi_InMe[i]);
		mod_(Mi[i], Me, &MiInMe[i]);

		mod_(M_, Base[i], &M_InBase[i]);

		printf("\n[%4.2f] #%i thread computation done.\n",(double)(clock() - begin) / CLOCKS_PER_SEC, i);
	}


	for (int i = 0; i < THREAD_NUMBER; i++){
		for (int j = 0; j < THREAD_NUMBER; j++){

			mod_( Mi [j], Base_ [i], &MiInBase_[i][j] );
			mod_( Mi_ [j], Base [i], &Mi_InBase[i][j] );

		}
	}

	printf("\n[%4.2f] &MiInBase_ and &Mi_InBase done.\n",(double)(clock() - begin) / CLOCKS_PER_SEC);

	mod_($N, Me, &N_Me);

	clearLongVariable( temp );
	
	//(M_^-1 mod Me) % Me;
	unsigned int t;
	InverseByModule_(M_, Me, &t);
	temp[ MAX - 1] = t;
	mod_(temp, Me, &InverseByModuleM_Me);
	
	printf("\n[%4.2f] &InverseByModuleM_Me done.\n",(double)(clock() - begin) / CLOCKS_PER_SEC);


	//(( M_ * Me ) ^ -1 mod Me) % Me;
	unsigned int tempM_Me[ MAX ];
	clearLongVariable( tempM_Me );
	clearLongVariable( temp );
	temp[MAX - 1] = Me;
	mul(M_, temp, tempM_Me);
	clearLongVariable( temp );
	InverseByModule(M, tempM_Me, temp);
	mod_(temp, Me, &InverseByModuleMMMe);
	
	printf("\n[%4.2f] &InverseByModuleMMMe done.\n",(double)(clock() - begin) / CLOCKS_PER_SEC);

	printf("\nTOTAL TIME SPENT: %f s\n",(double)(clock() - begin) / CLOCKS_PER_SEC);

}


int main(int argc, char **argv)
{
	
	//argv[1] - input message
	//argv[2] - action message
	//argv[3] - method message

	int mode; //encr/decr
	int unit; //CPU/GPU

	unsigned int $d[MAX];
	unsigned int $e[MAX];
	unsigned int q[MAX];
	unsigned int p[MAX];
	//unsigned int n[MAX];
	unsigned int p_minus_1[MAX];
	unsigned int q_minus_1[MAX];
	unsigned int phi[MAX];
	unsigned int temp[ MAX ];
	unsigned int $input[ MAX ];
	unsigned int R[THREAD_NUMBER], R_[THREAD_NUMBER];
	clock_t begin;
	cudaError_t cudaStatus;
	char output [ MAX * 10 ];

	if(
	   argv[2][0] == *("e") && 
	   argv[2][1] == *("n") && 
	   argv[2][2] == *("c") && 
	   argv[2][3] == *("r") && 
	   argv[2][4] == *("y") && 
	   argv[2][5] == *("p") && 
	   argv[2][6] == *("t") ) mode = 1;
   else if( 
	   argv[2][0] == *("d") && 
	   argv[2][1] == *("e") && 
	   argv[2][2] == *("c") && 
	   argv[2][3] == *("r") && 
	   argv[2][4] == *("y") && 
	   argv[2][5] == *("p") && 
	   argv[2][6] == *("t") ) mode = 2;
   else {
	   printf("Wrong mode: must be \"encrypt\" or \"decrypt\" (case sensitive) ");
	   goto Error;
   }
	
   if( argv[3][0] == *("C") && 
	   argv[3][1] == *("P") && 
	   argv[3][2] == *("U")) unit = 1;
   else if( 
	   argv[3][0] == *("G") && 
	   argv[3][1] == *("P") && 
	   argv[3][2] == *("U")) unit = 2;
   else {
	   printf("Wrong unit: must be \"CPU\" or \"GPU\" (case sensitive) ");
	   goto Error;
   }

  


	// Choose which GPU to run on, change this on a multi-GPU system.
    cudaStatus = cudaSetDevice(0);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
        goto Error;
    }

	
	/*// n = p*q
	mul ( p, q, n );

	clearLongVariable (  _1 );
	_1[ MAX - 1 ] = 1;		// for storing digit "1"
	clearLongVariable ( p_minus_1 );
	clearLongVariable ( q_minus_1 );

	// phi = (p-1)*(q-1);
	sub ( p, _1, p_minus_1 );
	sub ( q, _1, q_minus_1 );
	mul ( p_minus_1, q_minus_1, phi );

	

	*/
	
	Me = 64;
	/*
	convertToNormalForm( "130ebebd67b16a9ab2c53a437badbf8f01a80c750095a7fcfe95742c3d5ed1abb318babc5cb5d9350fee4da65ee074f65e1758117e6945f0fcfc8137528053ce9d1da8618890dee24e5e0bf8c87795bb1d09eddd544640824ee0dd0ea9fd908d27b0f8a1ae5c37f3647fbf2f5795500ad76c195b3387d0458a8f51b701472301" , $N);      // modulus
	convertToNormalForm( "0123" , $input);  // message
	convertToNormalForm( "010001" , $e); // public exponenta
	convertToNormalForm( "12e8da920d4599458e84ec5ef1656161807f427d05eb79182b7418259d6f6c14364d1f5caf9130c8d9d9d6ea71d1bdbc87781a46a16bcb9e672814fed3b9c96ddffe0a1b0955ae68055c8f92fef518a04fc32a2ea8390e617cc5556a251f9ae9eee70a32e579cb3e9f298848a9b3aaf634f5930ffbf74473f7cb6c0cefee1751" , $d); // secret exponenta 
	*/
	
	
	convertToNormalForm( "025123" , $N);      // modulus
	
	convertToNormalForm( "01365D" , $e); // public exponenta
	convertToNormalForm( "0AD" , $d) ; // secret exponenta 
	

	convertPowerToArrayOfBits($e, e, &eBitsCount, &numberOfIterationsForE);
	convertPowerToArrayOfBits($d, d, &dBitsCount, &numberOfIterationsForD);

	for (int i = 0; i < MAX; i++) {
		ZERO_MEMORY[ i ] = 0; 
	}

	

	for(int i=0; i<THREAD_NUMBER;i++){
		Base[i] = 4294967296 - HHH[i];
		Base_[i] = 4294967296 - KKK[i];
	}


	

	generateStaticData();
	

	convertToNormalForm( argv[1] , $input);  // message

	for(int z=0; z< THREAD_NUMBER;z++) {
		mod_($input, Base[z], &input[z]);
		mod_($input, Base_[z], &input_[z]);
	}


	if( mode == 1 ){

	//encrypt
	iterationsCounter = 1;
	currentOperationIndex = 0;
	begin = clock();
	for (int z = MAX * INT_SIZE - eBitsCount; z < MAX * INT_SIZE; z++){

		MM(ABmodN[0], ABmodN[1], ABmodN[0], ABmodN[1], R, R_);
		MM(R, R_, SQR_M_MOD_N_RNS, SQR_M_MOD_N_RNS_, ABmodN[0], ABmodN[1]);

		if (e[z]){

			MM(ABmodN[0], ABmodN[1], input, input_, R, R_);
			MM(R, R_, SQR_M_MOD_N_RNS, SQR_M_MOD_N_RNS_, ABmodN[0], ABmodN[1]);

		} 
	}
	if(unit == 2)
		printf("\nEncryption ended in : %f s\n",(double)(clock() - begin) / CLOCKS_PER_SEC / 15);
	else
		printf("\nEncryption ended in : %f s\n",(double)(clock() - begin) / CLOCKS_PER_SEC);

	};
	
	if( mode == 2 ){

	
	//decrypt
	iterationsCounter = 1;
	currentOperationIndex = 1;
	begin = clock();
	for (int z = MAX * INT_SIZE - dBitsCount; z < MAX * INT_SIZE; z++){

		MM(ABmodN[0], ABmodN[1], ABmodN[0], ABmodN[1], R, R_);
		MM(R, R_, SQR_M_MOD_N_RNS, SQR_M_MOD_N_RNS_, ABmodN[0], ABmodN[1]);

		if (d[z]){

			MM(ABmodN[0], ABmodN[1], input, input_, R, R_);
			MM(R, R_, SQR_M_MOD_N_RNS, SQR_M_MOD_N_RNS_, ABmodN[0], ABmodN[1]);

		} 
	}

	if(unit == 2)
		printf("\nDecryption ended in : %f s\n",(double)(clock() - begin) / CLOCKS_PER_SEC / 15);
	else
		printf("\nDecryption ended in : %f s\n",(double)(clock() - begin) / CLOCKS_PER_SEC);
	


	}

	convertFromRNS(ABmodN[0], temp);
	convertFromNormalForm( output, temp );
	printf("%s\n", output);


	showProfilerInfo();

	/*

	const int arraySize = 5;
    const int a[arraySize] = { 1, 2, 3, 4, 5 };
    const int b[arraySize] = { 10, 20, 30, 40, 50 };
    int c[arraySize] = { 0 };
    
	
    // Add vectors in parallel.
   
	for (int i=0; i< 100000; i++) {
	cudaError_t cudaStatus = addWithCuda(c, a, b, arraySize);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "addWithCuda failed!");
        return 1;
    }
	}

    printf("{1,2,3,4,5} + {10,20,30,40,50} = {%d,%d,%d,%d,%d}\n",
        c[0], c[1], c[2], c[3], c[4]);

    // cudaDeviceReset must be called before exiting in order for profiling and
    // tracing tools such as Nsight and Visual Profiler to show complete traces.
   /* cudaStatus = cudaDeviceReset();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceReset failed!");
        return 1;
    }
	*/
Error:
	getch();
    return 0;
}








