#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int main(int argc, char *argv[])
{
    FILE *file = fopen(argv[1], "r");

    char line[256];
    char *token;

    if (file == NULL)
    {
        printf("Invalid path. Can't find file!\n");
    }

    while( fgets( line, sizeof(line) , file) != NULL )
    {
        if ( line[0] == '%')     continue;
        
        printf("Line: %s", line);
        
        token = strtok(line, " \t\n");

        while ( token != NULL )
        {
            token = strtok(NULL, " \t\n");
            printf("%s\n", line);
        }

    }

    fclose(file);
    return 0;
}