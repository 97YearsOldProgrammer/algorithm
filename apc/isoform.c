#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// initial setup donor and acceptor site

int dons_count = 0;
int accs_count = 0;

// gt; ag acceptor donor reader

void splice_site_reader(const char *seq)
{
    static int ds = 0;
    static int ac = 0;

    for (int i = 0; i < strlen(seq) - 1; i++)
    {

        if (seq[i] == 'G' && seq[i+1] == 'T')
        {
            dons[ds] = i+1;
            ds++;
            dons_count++;

        } 

        else if (seq[i] =='A' && seq[i+1] == 'G') 
        {
            accs[ac] = i+1;
            ac++;
            accs_count++;
        }
    }
}

// asign space to donor and acceptor site based on seq length

void da_array_assigner(const char *seq)
{
    int len = strlen(seq);
    
    if (len <= 500)
    {
        int dons[20];
        int accs[20];
    }
    else if (len > 500 && len <= 1500)
    {
        int dons[50];
        int accs[50];
    }
    else if (len > 1500 && len <= 3000)
    {
        int dons[100];
        int accs[100];
    }
}

// isoform generator

/*
void all_possible_isoforms(){

}
*/

char seq[] = "ACGTTGACGTAAGTAAAGCAGCGCCACGAGTAAGAGTAACCGTTTACC";

int main(void){

    da_array_assigner(seq);
    splice_site_reader(seq);

    for (int i = 0; i < dons_count; i++)
    {
        printf("%d\t", dons[i]);
    }

    printf("\n");

    for (int i = 0; i < accs_count; i++)
    {
        printf("%d\t", accs[i]);
    }

    printf("\n");
    return 0;
}



