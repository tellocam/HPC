#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

void txt_gen(double p, double max_c){
  char file_name[64];
  char p_char[16] = (char)p + "_" + (char)max_c;   // number of processes for the filename
  p_char += _
  FILE *fp; // file pointer
  int i;

  // Create the file name using the variable
  snprintf(file_name, sizeof(file_name), "file_%s.txt", p_char);

  // Open the file for writing
  fp = fopen(file_name, "w");
  if (fp == NULL) {
    printf("Error opening file!\n");
    return 1;
  }

  // // Write some data to the file
  // fprintf(fp, "Variable: %s\n", var);
  // for (i = 0; i < 10; i++) {
  //     fprintf(fp, "Value: %d\n", i);
  }
  return 0;
}

int main(int argc, char *argv[]){



int count = 1;
int gentxt = 0;
for (int i=1; i<argc&&argv[i][0]=='-'; i++) {
    if (argv[i][1]=='c') i++, sscanf(argv[i],"%d",&count);
    if (argv[i][1]=='g') i++, sscanf(argv[i],"%d",&gentxt);
}

double test_vector[count];

for (int i = 0; i<count; i++){
    test_vector[i] = i;
}

for (int i = 0; i<count; i++){
    fprintf(stderr, "%.2f, %.2f \n", 
	test_vector[i], test_vector[i] );
}

}