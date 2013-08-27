#include "shifts.h"
#include "Constants.h"

// ������� 32-���� ����� ������ �� 1 �� � ��������� carry flag
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

// ������� 32-���� ����� ���� �� 1 �� � ��������� carry flag
int shl (unsigned int &input, int rightBit){
	
	int CF;
	CF = ( input & ( 1 << INT_SIZE - 1 ) ) == ( 1 << INT_SIZE - 1 );
	input <<= 1;             
	input &= 0xFFFFFFFE;	
	input |= rightBit;		
	return CF;
}


// ������� 64-���� ����� ������ �� 1 �� � ��������� carry flag
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

//������� "����� �����" ���� �� 1 �� � ��������� carry flag
int shiftToLeftVariable ( unsigned int *input ) {
	int CF = 0;
	for (int i = MAX - 1; i >= 0; i--) {
		CF = shl( input[ i ], CF );
	}
	return CF;
}

//������� "����� �����" ������ �� 1 �� � ��������� carry flag
int shiftToRightVariable ( unsigned int *input ) {
	int CF = 0;
	for (int i = 0; i < MAX; i++) {
		CF = shr( input[ i ], CF );
	}

	return CF;
}
