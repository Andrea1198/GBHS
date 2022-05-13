#!/bin/awk -f
{
   if ( match($0,"#") ) {print; next} 
   if ( NF != 4) next

   printf "%15.9f %15.9f \n", $1, $2+$4
}
