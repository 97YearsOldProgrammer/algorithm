#include "model.h"

void setup_initial_probability(Lambda *l)
{
    l->pi = calloc(HS, sizeof(double) );        // left-right HMM; only exon are 1
    l->pi[0] = 1;                               // initial probability of exon are 1
}

void setup_transition_probability(Lambda *l)
{
    l->A = malloc ( HS * sizeof(double*) );     // set up 13 spots in an array

    for (int i = 0 ; i < HS ; i++ )             // set up 13 x 13 matrix
    {
        l->A[i] = calloc(HS, sizeof(double) );  // assgin every element in this matrix to 0 initially
    }

    l->A[12][0] = 1.0;                          // a6 to exon is 1

    for ( int i = 1 ; i < 6 ; i++ )             // d(i) to d(i+1) is 1
    {
        l->A[i][i + 1] = 1.0;
    }

    for ( int i = 7 ; i < 12 ; i++ )
    {
        l->A[i][i + 1] = 1.0;                   // a(i) to a(i+1) is 1
    }
}
