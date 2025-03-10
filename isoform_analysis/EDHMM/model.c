#include "model.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void numerical_transcription(Observed_events *info, const char *seq)
{
    // turns original sequence into int 
    size_t len = strlen(seq);

    info->T = len;
    info->numerical_sequence = malloc ( len * sizeof(int) );

    for( size_t i = 0; i < len; i++)
    {

        // A == 0 , C == 1, G == 2, T == 3  
        if      (seq[i] == 'A')     info->numerical_sequence[i] = 0;
        else if (seq[i] == 'C')     info->numerical_sequence[i] = 1;
        else if (seq[i] == 'G')     info->numerical_sequence[i] = 2;
        else if (seq[i] == 'T')     info->numerical_sequence[i] = 3;
    }
}

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
    
    for (int i = 0; i < length; i++)
    {
        value += array[beg + i] * power(4, length - i - 1);
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

    alpha->a = malloc ( (info->T - 2 * FLANK + 1) * sizeof(double*) );                    // 2x2; outer layer as length of sequence

    for (int i = 0 ; i < (info->T - 2 * FLANK + 1) ; i++ )
    {
        alpha->a[i] = calloc( HS , sizeof(double) );                    // inner layer as number of Hidden States
    }
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

void initial_forward_algorithm(Lambda *l, Explicit_duration *ed,  Forward_algorithm *alpha, Observed_events *info)
{
    for ( int i = 0 ; i < HS ; i ++ )
    {
        double initial_prob = l->pi[i];
        
        double ed_prob;
        
        if      ( i == 0 ) ed_prob = ed->exon[0];
        else if ( i == 1 ) ed_prob = ed->intron[0];

        int index = base4_to_int(info->numerical_sequence, FLANK - 3, 4);
        double emission_prob;

        if       ( i == 0 ) emission_prob = l->B.exon[index];
        else if  ( i == 1 ) emission_prob = l->B.intron[index];

        double all = initial_prob * ed_prob * emission_prob;

        if ( all > 0 )  alpha->a[0][i] = exp( log(all) );
        else            alpha->a[0][i] = 0.0;

    }

}

void forward_algorithm(Lambda *l, Forward_algorithm *alpha, Observed_events *info, Explicit_duration *ed)
{
    int terminal_region_intron = info->T - ed->min_len_exon - FLANK + 1;                        // no intron continue beyond this point
    int terminal_region_exon   = info->T - ed->min_len_exon - ed->max_len_exon - FLANK + 1;     // no exon to intron beyond this point 

    for ( int t = 0 ; t < info->T - 2 * FLANK + 1 ; t ++ )                                      // iterate every element in the time scale
    {
         
        for ( int i = 0 ; i < HS ; i ++ )                                                       // the next position; 0 for exon, 1 for intron
        {
            int log_index = 0;                                                                  // prepare for log-space calcualtion
            int max_len;

            // biological restrain
            if ( t >= terminal_region_intron && i == 1)                                         // can't have new intron after this
            {                                                                                   // cuz we have to end as exon
                alpha->a[t][i] = 0.0;
                continue;
            }

            if      ( i == 0 )    max_len = ed->max_len_exon;                                   // max_len for exon state
            else if ( i == 1 )    max_len = ed->max_len_intron;                                 // max_len for intron state

            for ( int j = 0 ; j < HS ; j ++ )                                                   // the situation of transition out
            {
                if (i == j) continue;                                                           // no self-transition for HSMM

                // biological restrain
                if (t >= terminal_region_exon && j == 0 && i == 1) continue;                    // no exon->intron transition beyond this point
                
                for ( int d = 1 ; d <= max_len && t - d >= 0 ; d ++ )                           // skip those unsatisfied len; cuz they all 0
                {
                    // get product for emission prob

                    double emission_product = 1.0;                                              // what we want to calculate out here

                    for ( int s = t - d + 1 ; s <= t ;  s++ )
                    {
                        int index = base4_to_int(info->numerical_sequence, s - 3 + FLANK, 4);   // how we get index value from 4 base pair

                        double emission_prob;

                        if      ( i == 0 ) emission_prob = l->B.exon[index];                    // p emission exon
                        else if ( i == 1 ) emission_prob = l->B.intron[index];                  // p emission intron

                        emission_product *= emission_prob;                                      // update value
                    }

                    // get explicit duration prob

                    double ed_prob;

                    if      ( i == 0 ) ed_prob = ed->exon[d];
                    else if ( i == 1 ) ed_prob = ed->intron[d];

                    // get transition prob

                    double trans_prob;                                                          // transition prob

                    if ( j == 0 && i == 1)                                                      // from exon to intron
                    {
                        int index = base4_to_int(info->numerical_sequence , t - d + FLANK , 5); // 5bps motif
                        trans_prob = l->A.dons[index];
                    }
                    else if ( j == 1 && i == 0)                                                 // from intron to exon
                    {
                        int index = base4_to_int(info->numerical_sequence , t-d- 6 + FLANK, 6); // 6bps motif
                        trans_prob = l->A.accs[index];
                    }

                    // formal computation

                    double all = alpha->a[t - d][j] * trans_prob * ed_prob * emission_product;
                    if  ( all > 0 ) l->log_values[log_index++] = safe_log(all); 
                }
            }

            for ( int d = 1 ; d <= max_len && t - d + 1 >= 0 ; d ++ )                           // for continue probability
            {
                // get emission product
                double emission_product = 1.0;

                for ( int s = t - d + 1 ; s <= t ;  s++ )
                {
                    int index = base4_to_int(info->numerical_sequence, s - 3 + FLANK, 4);       // how we get index value from 4 base pair

                    double emission_prob;

                    if      ( i == 0 ) emission_prob = l->B.exon[index];                        // p emission exon
                    else if ( i == 1 ) emission_prob = l->B.intron[index];                      // p emission intron
                    emission_product *= emission_prob;                                          // update value
                }

                // get explicit duration probability

                double ed_prob;

                if      ( i == 0 ) ed_prob = ed->exon[d];
                else if ( i == 1 ) ed_prob = ed->intron[d];

                double all = alpha->a[t - d][i] * ed_prob * emission_product;
                if ( all > 0 ) l->log_values[log_index++] = safe_log(all);
            }

            alpha->a[t][i] = exp( log_sum_exp (l->log_values, log_index) );
        }
    }
}

void allocate_beta(Observed_events *info, Backward_algorithm *beta)                             // assign data structure for backward algorithm
{

    beta->b = malloc ( (info->T - 2 * FLANK + 1) * sizeof(double*) );                           // 2x2; outer layer as length of sequence

    for (int i = 0 ; i < (info->T - 2 * FLANK + 1) ; i++ )
    {
        beta->b[i] = calloc( HS , sizeof(double) );                                             // inner layer as number of Hidden States
    }
}

void initial_backward_algorithm(Lambda *l, Backward_algorithm *beta, Observed_events *info)
{
    int last_bps = info->T - 2 * FLANK; 

    // constrain ; only possible to end as an exon 
    beta->b[last_bps][0] = 1.0;
    beta->b[last_bps][1] = 0.0;
}

void backward_algorithm(Lambda *l, Backward_algorithm *beta, Observed_events *info, Explicit_duration *ed)
{
    int start_region_intron = info->T - ed->min_len_exon - FLANK + 1;

    for ( int t = info->T - FLANK ; t >= FLANK ; t-- )
    {

        for ( int i = 0 ; i < HS ; i ++ )                                                       // the before position; 0 for exon, 1 for intron
        {
            int log_index = 0;                                                                  // prepare for log-space calcualtion
            int max_len;

            // biological restrain
            if ( t >= start_region_intron && i == 1)                                            // can't have new intron after this
            {                                                                                   // cuz we have to end as exon
                beta->b[t][i] = 0.0;
                continue;
            }

            if      ( i == 0 )    max_len = ed->max_len_exon;                                   // max_len for exon state
            else if ( i == 1 )    max_len = ed->max_len_intron;                                 // max_len for intron state

            for ( int j = 0 ; j < HS ; j ++ )                                                   // the situation of transition out
            {
                if (i == j) continue;                                                           // no self-transition for HSMM
                
                for ( int d = 1 ; d <= max_len && t + d >= 0 ; d ++ )                           // skip those unsatisfied len; cuz they all 0
                {

                    double emission_product = 1.0;                                              // what we want to calculate out here

                    for ( int s = t + 1 ; s <= t + d ;  s++ )
                    {
                        int index = base4_to_int(info->numerical_sequence, s - 3, 4);           // how we get index value from 4 base pair

                        double emission_prob;

                        if      ( i == 0 ) emission_prob = l->B.exon[index];                    // p emission exon
                        else if ( i == 1 ) emission_prob = l->B.intron[index];                  // p emission intron

                        emission_product *= emission_prob;                                      // update value
                    }

                    // get explicit duration prob

                    double ed_prob;

                    if      ( i == 0 ) ed_prob = ed->exon[d];
                    else if ( i == 1 ) ed_prob = ed->intron[d];

                    // get transition prob

                    double trans_prob;                                                          // transition prob

                    if ( j == 0 && i == 1)                                                      // from exon to intron
                    {
                        int index = base4_to_int(info->numerical_sequence , t - d , 5);         // 5bps motif
                        trans_prob = l->A.dons[index];
                    }
                    else if ( j == 1 && i == 0)                                                 // from intron to exon
                    {
                        int index = base4_to_int(info->numerical_sequence , t - d - 6 , 6);     // 6bps motif
                        trans_prob = l->A.accs[index];
                    }

                    // formal computation

                    double all = beta->b[t + d][j] * trans_prob * ed_prob * emission_product;
                    if  ( all > 0 ) l->log_values[log_index++] = safe_log(all); 
                }
            }

            for ( int d = 1 ; d <= max_len && t + d >= 0 ; d ++ )                               // for continue probability
            {
                // get emission product
                double emission_product = 1.0;

                for ( int s = t + 1 ; s <= t + d ;  s++ )
                {
                    int index = base4_to_int(info->numerical_sequence, s - 3, 4);               // how we get index value from 4 base pair

                    double emission_prob;

                    if      ( i == 0 ) emission_prob = l->B.exon[index];                        // p emission exon
                    else if ( i == 1 ) emission_prob = l->B.intron[index];                      // p emission intron
                    emission_product *= emission_prob;                                          // update value
                }

                // get explicit duration probability

                double ed_prob;

                if      ( i == 0 ) ed_prob = ed->exon[d];
                else if ( i == 1 ) ed_prob = ed->intron[d];

                double all = beta->b[t + d][i] * ed_prob * emission_product;
                if ( all > 0 ) l->log_values[log_index++] = safe_log(all);
            }

            beta->b[t][i] = exp( log_sum_exp (l->log_values, log_index) );
        }
    }
}

void free_alpha(Observed_events *info, Forward_algorithm *alpha)
{
    for (int i = 0 ; i < (info->T - 2 * FLANK + 1); i ++)
    {
        free( alpha->a[i] );
    }
    free( alpha->a ); 
}
void free_beta(Observed_events *info, Backward_algorithm *beta)
{
    for (int i = 0 ; i < (info->T - 2 * FLANK  + 1) ; i ++)
    {
        free( beta->b[i] );
    }
    free( beta->b );
}

// output //

// posterior prob //
void print_posterior_probabilities(Observed_events *info, Forward_algorithm *fw, Backward_algorithm *bw, int start_pos, int end_pos) 
{
    printf("\n=== Posterior Probabilities ===\n");
    printf("Position\tExon\t\tIntron\n");
    printf("----------------------------------------\n");
    
    // Calculate total observation probability for normalization
    double total_prob = 0.0;
    for (int i = 0; i < 2; i++) { // Just sum exon and intron states
        total_prob += fw->a[info->T - 2 * FLANK][i];
    }
    
    // Print posterior probabilities for the requested range
    for (int t = start_pos; t <= end_pos && t < info->T - 2 * FLANK + 1; t++) {
        printf("%d\t\t", t + FLANK);
        
        for (int i = 0; i < 2; i++) { // Just print for exon (0) and intron (1)
            double posterior = (fw->a[t][i] * bw->b[t][i]) / total_prob;
            printf("%.6f\t", posterior);
        }
        printf("\n");
    }
    printf("----------------------------------------\n");
}