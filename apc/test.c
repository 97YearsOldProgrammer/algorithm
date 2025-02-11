#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char a[] = "abcd";

int main(void){

    for (int i = 0 ; i < strlen(a); i++)
    {
        printf("%c", a[i]);
    }
    printf("\n");
    return 0;
}