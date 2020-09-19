/*
 * ExecuteLeuteron.C
 *
 *  Created on: 28 January 2020
 *      Author: Michael Jung
 */

#include "GetCorrelationsLeuteron.C"

int main(int argc, char* argv[]) {

  const char* file = argv[1];
  const char* prefix = (argv[2]) ? argv[2] : "";
  const char* addon = (argv[3]) ? argv[3] : ""; 

  GetCorrelationsLeuteron(file,"","");

  return 1;

}
