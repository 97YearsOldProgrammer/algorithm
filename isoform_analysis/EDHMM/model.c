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

int base4_to_int(int *array, int beg, int length) 
{
    int value = 0;
    
    for (int i = beg ; i < length ; i++ )
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

    return value;
}

void initialize_donor_transition_matrix(Lambda *l, Apc *a, int depth)   // set the depth to 0 initially
{   
    if (depth == 5)                                                     // exit recursion; calculation start
    {
        int index = base4_to_int(a->position, 0, 5);                    // this is where we plan to store that value
        double value = total_prob(a->prob, 5);                          // get total prob
        l->A.dons[index]  = value;                                      // store the value
        return;
    }

    for ( int i = 0; i < 4 ; i++ )
    {
        double p = l->B.dons[depth][i];                                 // get the exon emission prob for current node
        a->prob[depth] = p;                                             // record the emission prob inside the array
        a->position[depth] = i;                                         // record each base pair each node choose
        initialize_donor_transition_matrix(l, a, depth + 1);            // send into next node
    }
}

void initialize_acceptor_transition_matrix(Lambda *l, Apc *a, int depth)// set the depth to 0 initially
{   
    if (depth == 6)                                                     // exit recursion; calculation start
    {
        int index = base4_to_int(a->position, 0, 6);                    // this is where we plan to store that value
        double value = total_prob(a->prob, 6);                          // get total prob
        l->A.accs[index]  = value;                                      // store the value
        return;
    }

    for ( int i = 0; i < 4 ; i++ )
    {
        double p = l->B.accs[depth][i];                                 // get the exon emission prob for current node
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

int initialize_forward_algorithm(Forward_algorithm *alpha, Explicit_duration *ed)
{
    int i = 0;
    while ( ed->exon[i] == 0.0)
    {
        alpha->a[i][0] = 0.0;                                                                   // set invalid exon vlaue
        alpha->a[i][1] = 0.0;                                                                   // set invalid intron value ; actually missing 5 bps here
        i ++;    
    }
    return i;                                                                                   // we wanna know where is the valid start point
}

double safe_log(double x)                                                                       // handle log 0
{
    const double epsilon = 1e-10;                                                               // log 0 gonna crush in C, name a value for it
    return log(x + epsilon);                                                                    // return huge negative log value
}

double log_sum_exp(double *logs, int n) 
{
    double max_log = logs[0];

    for (int i = 1; i < n ; i++)                                                                // get local maximum
    {
        if( logs[i] > max_log )         max_log = logs[i];                  
    }

    double sum = 0.0;

    for( int i = 0 ; i < n ; i++)       sum += exp(logs[i] - max_log);                          // log soft max trick

    return max_log + log(sum);       
}

void forward_algorithm(Lambda *l, Forward_algorithm *alpha, Observed_events *info, Explicit_duration *ed)
{
    int beg = initialize_forward_algorithm(alpha, ed);                                          // which time t we wanna continue

    for ( t = beg ; t < info->T ; t ++ )                                                        // iterate every element in the time scale
    {
        int current = info->numerical_sequence[t];                                              // getting the current obs event
        int before  = info->numerical_sequence[t - 3];                                          // prepared for the exon or intron
        int index   = base4_to_int(info->numerical_sequence, before, current);                  // value for emission prob for intron or exon                   

        for ( i = 0 ; i < HS ; i ++ )                                                           // the next position; 0 for exon, 1 for intron
        {
            int log_index = 0;                                                                  // prepare for log-space calcualtion

            if (i = 0)
            {
                int max_len        = ed->max_len_exon;                                          // max_len for exon state
                int transition_len = 5;                                                         // 5 bps for donor site 
            }
            else if (i = 1)
            {
                int max_len = ed->max_len_intron;                                               // max_len for intron state
                int transition_len = 6;                                                         // 6 bps for acceptor site
            }

            for ( j = 0 ; j < HS ; j ++ )                                                       // the situation of transition out
            {
                if (i = j) continue;                                                            // no self-transition for HSMM
                
                for ( int d = 1 ; d <= ed->max_len && t - d - transition_len + 1 >= 0 ; d ++ )
                {
                    // get product for emission prob

                    double emission_product = 1.0;                                              // what we want to calculate out here

                    for ( int s = t - d + 1 ; s <= t ;  s++ )
                    {
                        int index = base4_to_int(info->numerical_sequence, s - 3, 4);           // how we get index value from 4 base pair

                        double emission_prob;

                        if      ( i = 0 ) emission_prob = l->B.exon[index];                     // p emission exon
                        else if ( i = 1 ) emission_prob = l->B.intron[index];                   // p emission intron

                        emission_product *= emission_prob;                                      // update value
                    }

                    // get explicit duration prob

                    if      ( i = 0 ) ed_prob = ed->exon[d];
                    else if ( i = 1 ) ed_prob = ed->intron[d];

                    // get transition prob

                    double trans_prob;                                                          // transition prob

                    if ( j == 0 && i == 1)                                                      // from exon to intron
                    {
                        int index = base4_to_int(info->numerical_sequence , t - d - 4 , 5);     // 5bps motif
                        trans_prob = l->A.dons[index];
                    }
                    else if ( j == 1 && i == 0)                                                 // from intron to exon
                    {
                        int index = base4_to_int(info->numerical_sequence , t - d - 5 , 6);     // 6bps motif
                        trans_prob = l->A.accs[index];
                    }

                    // formal computation

                    double all = alpha->a[t - d][j] * trans_prob * duration_prob * emission_product;
                    if  ( all > 0 ) alpha->log_values[log_index] = safe_log(all); 
                    log_index ++;
                }
                
            }

            for ( int d = 1 ; d <= ed->max_len && t - transition_len - d + 1 >= 0 ; d ++ )      // for continue probability
            {
                // get emission product
                double emission_product = 1.0;

                for ( int s = t - d + 1 ; s <= t ;  s++ )
                {
                    int index = base4_to_int(info->numerical_sequence, s - 3, 4);               // how we get index value from 4 base pair

                    double emission_prob;

                    if      ( i = 0 ) emission_prob = l->B.exon[index];                         // p emission exon
                    else if ( i = 1 ) emission_prob = l->B.intron[index];                       // p emission intron

                    emission_product *= emission_prob;                                          // update value
                }

                // get explicit duration probability
                if      ( i = 0 ) ed_prob = ed->exon[d];
                else if ( i = 1 ) ed_prob = ed->intron[d];

                double all = alpha->a[t - d][i] * ed_prob * emission_product;
                if ( all > 0 ) alpha->log_values[log_index] = safe_log(all);
                log_index ++
            }

            alpha->a[t][i] = exp(log_sum_exp (alpha->log_values, log_index))
        }
    }
}