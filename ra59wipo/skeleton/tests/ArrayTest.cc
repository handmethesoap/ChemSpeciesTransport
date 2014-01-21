
#include "Array.hh"
#include "Debug.hh"

#include <iostream>

void copyTest()
{
   const int size = 10;

   // Fill single array...
   Array<real> arr (size);
   for ( int i = 0; i<size; ++i )
      arr(i) = i;

   // check if values where set correctly
   for ( int i = 0; i<size; ++i )
      CHECK( arr(i) == i );


   // create a copy of the array, and check if values are the same
   Array<real> arrCopy (arr);
   for ( int i = 0; i<size; ++i )
      CHECK( arrCopy(i) == i );


   // write new values to copy
   for ( int i = 0; i<size; ++i )
      arrCopy(i) = 2*i;

   // and check both arrays again
   for ( int i = 0; i<size; ++i )
      CHECK( arr(i) == i );
   for ( int i = 0; i<size; ++i )
      CHECK( arrCopy(i) == 2*i );

}

void  contiguousMemoryTest()
{
   Array<real> arr ( 10, 5 );

   size_t memDifference =  & ( arr (9,4 ) ) - & ( arr(0,0) );
   CHECK( memDifference == 10*5 - 1 );
}


int main( )
{
   std::cout << "Copy Test: ";
   copyTest();
   std::cout << "OK" << std::endl;

   std::cout << "Contiguous Memory Test: ";
   contiguousMemoryTest();
   std::cout << "OK" << std::endl;

   return 0;
}
