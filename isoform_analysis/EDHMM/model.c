#include "model.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void setup_initial_probability(Lambda *l)                               // actually no longer needed
{
    l->pi = calloc(HS, sizeof(double) );                                // left-right HMM; only exon are 1
    l->pi[0] = 1;                                                       // initial probability of exon are 1
}

int power(int base, int exp)                                            // wtf, C don't have power for int
{
    int result = 1;

    for ( int i = 0 ; i < exp ; i++ )
    {
        result *= base;
    }
    return result;
}

int base4_to_int(int *array, int length) 
{
    int value = 0;
    
    for (int i = 0 ; i < length ; i++ )
    {
        value += array[i] * power(4, length - i - 1);
    }
    
    return value;
}

double total_prob(double *array, int length)
{
    double value = 1.0;

    for (int i = 0 ; i < length ; i ++)
    {
        value *= array[i];
    }

    return value
}

void initialize_donor_transition_matrix(Lambda *l, Apc *a, int depth)   // set the depth to 0 initially
{   
    if (depth == 5)                                                     // exit recursion; calculation start
    {
        int index = base4_to_int(a->position, 4);                       // this is where we plan to store that value
        double value = total_prob(a->prob, 4);                          // get total prob
        l->A.dons[index]  = value;                                      // store the value
        return;
    }

    for ( i = 0; i < 4 ; i++ )
    {
        double p = emission_prob[depth][i];                             // get the exon emission prob for current node
        a->prob[depth] = p;                                             // record the emission prob inside the array
        a->position[depth] = i;                                         // record each base pair each node choose
        initialize_donor_transition_matrix(l, a, depth + 1);            // send into next node
    }
}

void initialize_acceptor_transition_matrix(Lambda *l, Apc *a, int depth)// set the depth to 0 initially
{   
    if (depth == 6)                                                     // exit recursion; calculation start
    {
        int index = base4_to_int(a->position, 5);                       // this is where we plan to store that value
        double value = total_prob(a->prob, 5);                          // get total prob
        l->A.accs[index]  = value;                                      // store the value
        return;
    }

    for ( i = 0; i < 4 ; i++ )
    {
        double p = emission_prob[depth][i];                             // get the exon emission prob for current node
        a->prob[depth] = p;                                             // record the emission prob inside the array
        a->position[depth] = i;                                         // record each base pair each node choose
        initialize_acceptor_transition_matrix(l, a, depth + 1);         // send into next node
    }
}

void allocate_alpha(Observed_events *info, Forward_algorithm *alpha)    // assign data structure for alpha
{

    alpha->a = malloc ( info->T * sizeof(double*) );                    // 2x2; outer layer as length of sequence

    for (int i = 0 ; i < info->T ; i++ )
    {
        alpha->a[i] = calloc( HS , sizeof(double) );                    // inner layer as number of Hidden States
    }
}

void basis_of_forward_algo(Forward_algorithm *alpha, Explicit_duration *ed)   
{                                                                                               // left to right HMM; only exon need to assign
    int i = 0;
    while (ed->ed_exon[i] == 0.0)
    {

    }
    double initial_prob = log( l->pi[i] );                                                      // get initial probability
    double p;

    int current  = FLANK;                                                                       // start bps
    double prob_ed_exon = log( ed->ed_exon[0] );                                                // explicit duration at length at 1
    int before = current - 3;                                                                   // 3 bps before start
    int index = emission_probability_accessor(l, info , before , current );                     // index of explicit duration
    double prob_emission_exon = log( ed->ed_exon[index] );                                      // pointing to designated prob

    p = initial_prob +  prob_ed_exon + prob_emission_exon;                                      // everything in log space  
    alpha->a[0][i] = p;                                                                         // store as alpha

    // skip the rest of initial alpha for our model, since it's purely left-right-HMM
    // we don't consider any other prob for rest of hidden state as initial
}

int initialize_forward_algorithm(Forward_algorithm *alpha, Explicit_duration *ed)
{
    int i = 0;
    while ( ed->exon[i] == 0.0)
    {
        alpha->a[i][0] = 0.0;                                           // set invalid exon vlaue
        alpha->a[i][1] = 0.0;                                           // set invalid intron value ; actually missing 5 bps here
        i ++;    
    }
    return i;                                                         // we wanna know where is the valid start point
}

void forward_algorithm(Lambda *l, Forward_algorithm *alpha, Observed_events *info, Explicit_duration *ed)
{
    int beg = initialize_forward_algorithm(alpha, ed);                                          // which time t we wanna continue

    for ( t = 1 ; t < info->T ; t ++ )                                                          // iterate every element in the time scale
    {
        int real_t  = FLANK + t;                                                                // current bps is on flank + t
        int current = info->numerical_sequence[t];                                              // getting the current obs event
        int before  = info->numerical_sequence[t - 3];                                          // prepared for the exon or intron

        for ( i = 1 ; i < HS ; i ++ )
        {
             
        }
    }
}