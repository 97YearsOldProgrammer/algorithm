#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "model.h"


// emission probability //

void donor_parser(Lambda *l, char *filename)            // get emission probability for donor site
{
    FILE *file = fopen(filename, "r");

    char line[256];
    char *token;
    double p;                                           // probability we are going to store

    if (file == NULL)
    {
        printf("Can't find file for donor site emission probability!\n");
        return;
    }
    
    int c_line = -1;                                    // count of line

    while( fgets( line, sizeof(line) , file) != NULL )  // nest while loop to get elements
    {

        if ( line[0] == '%')     continue;              // skip the first line        

        c_line++;
        int c_token = -1;

        token = strtok(line, " \t\n");

        while ( token != NULL )
        {
            c_token ++;
            p = atof(token);                            // convert string into double
            l->B.dons[c_line][c_token] = p;             // 
            token = strtok(NULL, " \t\n");              // move to next element
        }
    }
    fclose(file);
}

void acceptor_parser(Lambda *l, char *filename)         // get emission probability for acceptor site
{
    FILE *file = fopen(filename, "r");

    char line[256];
    char *token;
    double p;                                           // probability we are gonna store

    if (file == NULL)
    { 
        printf("Can't find file for donor site emission probability!\n");
        return;
    }
    
    int c_line = -1;                                    // count of line

    while( fgets( line, sizeof(line) , file) != NULL )  // nest while loop to get elements
    {
        if ( line[0] == '%')     continue;              // skip the first line   

        c_line++;
        int c_token = -1;

        token = strtok(line, " \t\n");                  // get each probability

        while ( token != NULL )                         
        {
            c_token ++;
            p = atof(token);                            // convert string into double
            l->B.accs[c_line][c_token] = p;             // store the value
            token = strtok(NULL, " \t\n");              // move to next element
        }
    }
    fclose(file);
}

void exon_intron_parser(Lambda *l, char *filename)      // get emission probability for exon and intron
{
    FILE *file = fopen(filename, "r");

    char line[256];
    char *token;
    double p;                                           // probability we are gonna store

    if (file == NULL)
    {
        printf("Invalid emission probability for exon or intron. Can't find file!\n");
        return;
    }

    int c_line = -1;                                    // count of line

    while( fgets( line, sizeof(line) , file) != NULL )  // nest while loop to get elements
    {
        if ( line[0] == '%')     continue;              // skip the first line        

        int c_line ++;

        token = strtok(line, " \t\n");                  // 4 base pair from the file
        token = strtok(NULL, " \t\n");                  // get the probability 

        p = atof(token);
        l->B.exon[c_line] = p;                          // store the probability
    }
    fclose(file);
}

void explicit_duration_probability(explicit_duration *ed, char *filename, int digit)
{
    assert( 0 <= digit <= 1);

    FILE *file = fopen(filename, "r")

    char line[256];
    char *token;
    double p;

    if (file == NULL)
    {
        printf("Can't find the file for the explicit duration.");
        return;
    }

    int c_line = -1;

    while( fgets( line, sizeof(line) , file) != NULL)
    {
        if (line[0] == '%')     continue;

        int c_line ++;

        token = strtok(line, " \t\n");
        p = atof(token);

        if(digit == 1)  ed->ed_exon[c_line]   = p;
        else            ed->ed_intron[c_line] = p;
    }
    fclose(file);
}