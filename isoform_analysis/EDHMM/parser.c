#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "model.h"

// read the sequence from file

void read_sequence_file(const char *filename, Observed_events *info)
{
    FILE *file = fopen(filename, "r");
    
    if (file == NULL)
    {
        printf("Error: Cannot open sequence file %s\n", filename);
        return;
    }
    
    fseek(file, 0, SEEK_END);
    long file_size = ftell(file);
    rewind(file);

    char *buffer = (char*)malloc(file_size + 1);
    size_t read_size = fread(buffer, 1, file_size, file);

    buffer[read_size] = '\0';
    char *sequence = (char*)malloc(file_size + 1);

    size_t seq_index = 0;
    
    // Parse each line
    char *line = strtok(buffer, "\n");
    while (line != NULL)
    {
        
        // Extract valid DNA characters
        for (int i = 0; line[i] != '\0'; i++)
        {
            if (line[i] == 'A' || line[i] == 'C' || line[i] == 'G' || line[i] == 'T')
            {
                sequence[seq_index++] = line[i];
            }
        }
        
        line = strtok(NULL, "\n");
    }
    
    sequence[seq_index] = '\0';
        sequence = (char*)realloc(sequence, seq_index + 1);
    
    info->original_sequence = sequence;
    info->T = seq_index; 
    
    free(buffer);
    fclose(file);
}

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

void exon_intron_parser(Lambda *l, char *filename, int digit)      // get emission probability for exon and intron
{
    assert(digit == 0 || digit == 1);                   // 0  for exon, 1 for intron

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

        c_line ++;

        token = strtok(line, " \t\n");                  // 4 base pair from the file
        token = strtok(NULL, " \t\n");                  // get the probability 

        p = atof(token);
        if(digit == 0)  l->B.exon[c_line]   = p;        // store the prob for exon
        else            l->B.intron[c_line] = p;        // store the prob for intron
    }
    fclose(file);
}

// eplicit_duration //

void explicit_duration_probability(Explicit_duration *ed, char *filename, int digit)
{
    assert(digit == 0 || digit == 1);                   // 0 for exon, 1 for intron

    FILE *file = fopen(filename, "r");

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

        c_line ++;

        token = strtok(line, " \t\n");
        p = atof(token);

        if      (p != 0.0 && digit == 0) ed->min_len_exon = c_line;        // update min exon   len
        else if (p != 0.0 && digit == 1) ed->min_len_intron = c_line;      // update min intron len
        
        if(digit == 0)  ed->exon[c_line]   = p;
        else            ed->intron[c_line] = p;
    }

    if (digit == 0)     ed->max_len_exon   = c_line;
    else                ed->max_len_intron = c_line;

    fclose(file);
}