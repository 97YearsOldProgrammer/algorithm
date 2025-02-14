#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

// ===================    global value  ==================== 

int *dons = NULL;
int *accs = NULL;
int dons_count = 0;
int accs_count = 0;

// remember to set global value for min-intron and min-exon

int minin = 2;
int minex = 2;
int flank = 2;

/*******************\
prototype being used
\*******************/

void splice_site_reader(const char *seq);
void da_array_assigner(const char *seq);
void pointer_printer(const int *start, const int *end);
void all_isoform(const int *donor, const int *donor_end, const int *acceptor, const int *acceptor_end, const int spot);


// ====================  test sequence ==================== 

char seq[] = "ACGTTGACGTAAGTAAAGCAGCGCCACGAGTAAGAGTAACCGTTTACC";


int main(void){

    da_array_assigner(seq);
    splice_site_reader(seq);

    int *ds_start = dons;
    int *ac_start = accs;
    int *ds_end   = ds_start + dons_count - 1;
    int *ac_end   = ac_start + accs_count - 1;

    pointer_printer(ds_start, ds_end);

    printf("\n");

    pointer_printer(ac_start, ac_end);

    printf("\n"); 

    all_isoform(ds_start, ds_end, ac_start, ac_end , 0);

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

    for (int i = flank + minex ; i < len - 1 - flank - minex; i++)
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

/*******************************\
printer session for easier debug
\*******************************/

// print donor and acceptor array with pointer \\

void pointer_printer(const int *start, const int *end)
{

    for (const int *p = start; p <= end; p++ )
    {
        printf("%d\t", *p);
    }

    printf("\n");
    
}

// combinator


void all_isoform(const int *donor, const int *donor_end, const int *acceptor, const int *acceptor_end, int spot)
{

    assert(spot % 2 == 0);
    static int isoform[20];

    // whenever we are out of donor and acceptor, exit
    if ( donor > donor_end || acceptor > acceptor_end)
    {
        return;
    }

    // creating a nested loop for imitating the intron formation
    for (const int *p1 = donor; p1 <= donor_end; p1++ )
    {
        // check the middle exon
        if (spot != 0 && *p1 - isoform[spot - 1] < minex + 1)
        {
            continue;
        }

        // ==================================================================\\
        // add whatever features we want to continue the loop if bad isoform \\
        // ==================================================================\\

        isoform[spot] = *p1;
        

        // picking acceptors
        for (const int *p2 = acceptor; p2 <= acceptor_end; p2++)
        {
            // check the middle intron length
            if (*p2 - *p1 < minin - 1)
            {
                continue;
            }

            isoform[spot + 1] = *p2;
            
            pointer_printer( isoform, isoform + spot + 1 );

            // ====================================================================== \\
            // here we could add the mRNA function which used to store and check then \\
            // ====================================================================== \\            

            // recursion, continue picking
            all_isoform(p1, donor_end, p2, acceptor_end, spot + 2);
        }
    }   
    
}


