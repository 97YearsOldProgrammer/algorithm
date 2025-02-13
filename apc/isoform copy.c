#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// ===================    global value  ==================== 

int *dons = NULL;
int *accs = NULL;
int dons_count = 0;
int accs_count = 0;

/*******************\
prototype being used
\*******************/

void splice_site_reader(const char *seq);
void da_array_assigner(const char *seq);
void pointer_printer(const int *start, const int *end);


// ====================  test sequence ==================== 

char seq[] = "ACGTTGACGTAAGTAAAGCAGCGCCACGAGTAAGAGTAACCGTTTACC";


int main(void){

    da_array_assigner(seq);
    splice_site_reader(seq);

    int *ds_start = dons;
    int *ac_start = accs;

    pointer_printer(ds_start, ds_start + dons_count - 1);
    printf("\n");
    pointer_printer(ac_start, ac_start + accs_count - 1);

    free(dons);
    free(accs);

    return 0;

}


// gt; ag acceptor donor reader

void splice_site_reader(const char *seq)

{
    dons_count = 0;
    accs_count = 0;

    size_t len = strlen(seq);

    for (int i = 0; i < len - 1; i++)
    {

        if (seq[i] == 'G' && seq[i+1] == 'T')
        {
            dons[dons_count] = i;
            dons_count++;

        } 

        else if (seq[i] =='A' && seq[i+1] == 'G') 
        {
            accs[accs_count] = i;
            accs_count++;
        }
    }
}

// asign space to donor and acceptor site based on seq length

void da_array_assigner(const char *seq)

{
    size_t len = strlen(seq);

    free(dons);
    free(accs);
    
    if (len <= 500)

    {
        dons = malloc ( 20 * sizeof(int) );
        accs = malloc ( 20 * sizeof(int) );
    }

    else if (len > 500 && len <= 1500)

    {
        dons = malloc ( 100 * sizeof(int) );
        accs = malloc ( 100 * sizeof(int) );
    }

    else if (len > 1500 && len <= 3000)

    {
        dons = malloc ( 150 * sizeof(int) );
        accs = malloc ( 150 * sizeof(int) );
    }

}

// print donor and acceptor array with pointer

void pointer_printer(const int *start, const int *end)

{

    for (const int *p = start; p <= end; p++ )

    {
        printf("%d\t", *p);
    }
    
}
