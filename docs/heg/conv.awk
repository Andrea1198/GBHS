#! /bin/awk -f

BEGIN{ ev2ry=13.6058 }

{

  if ( match($0,"#") ) {print; next}
  if ( NF != 2 ) {print; next}
 
  printf "%15.9f %15.9f \n", $1, $2/ev2ry

}
