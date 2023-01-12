#pragma GCC diagnostic ignored "-Wformat-truncation"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#pragma GCC diagnostic ignored "-Wformat-truncation"

int main(int argc, char *argv[]){



int count = 1;
int gentxt = 0;
for (int i=1; i<argc&&argv[i][0]=='-'; i++) {
    if (argv[i][1]=='c') i++, sscanf(argv[i],"%d",&count);
    if (argv[i][1]=='g') i++, sscanf(argv[i],"%d",&gentxt);
}
    FILE *fp; // file pointer
if (gentxt!=0){
    char file_name[128];
    char file_suffix[64];
    char p_char[16];
    char c_char[16];
    char uline[16] = "_";
    //snprintf(p_char, sizeof(p_char), "%X", count);
    sprintf(p_char, "%X", count);
    //snprintf(c_char, sizeof(c_char), "%X", count);
    sprintf(c_char, "%X", count);
    strcat(file_suffix, p_char);
    strcat(file_suffix, uline);
    strcat(file_suffix, c_char);
    //snprintf(file_name, sizeof(file_name), "Ex1_%s.txt", file_suffix);
    sprintf(file_name, "Ex1_%s.txt", file_suffix);
    fp = fopen(file_name, "w");
    if (fp == NULL) {
    printf("Error opening file!\n");
    }
    fprintf(stderr, "%s \n", file_name);
}
double test_vector[count];

for (int i = 0; i<count; i++){
    test_vector[i] = i;
}



for (int i = 0; i<count; i++){
    fprintf(stderr, "%.2f, %.2f \n", 
	test_vector[i], test_vector[i] );
    if (gentxt!=0){
        fprintf(fp, "%.2f, %.2f \n", test_vector[i], test_vector[i]);
    }

}
fclose(fp);
}