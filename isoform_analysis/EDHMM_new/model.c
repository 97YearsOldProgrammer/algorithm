#include "model.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void numerical_transcription(Observed_events *info, const char *seq)
{
    printf("Start transforming original sequence into base4:\n");

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

    printf("\tWe get numerical sequence with Seq len: %d\n", info->T);
    printf("\t\u2713\n");
}

void setup_initial_probability(Lambda *l)                               // actually no longer needed
{
    printf("Start getting initial probability down:");
    l->pi = calloc(HS, sizeof(double) );                                // left-right HMM; only exon are 1
    l->pi[0] = 1;                                                       // initial probability of exon are 1
    printf("\t\u2713\n");
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
    int values = 0;
    int value;
    for (int i = 0; i < length; i++)
    {
        value  =  array[beg + i];
        values += value * power(4, length - i - 1);
    }
    
    return values;
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

void normalize_transition_prob(Lambda *l, int len, int dons_or_accs)
{
    printf("Start normalizing transition prob:");

    double sum = 0.0;
    if (dons_or_accs == 0)
    {
        for ( int i = 0 ; i < len ; i ++ )  sum += l->A.dons[i];
        for ( int i = 0 ; i < len ; i ++ )  l->A.dons[i] /= sum;
    }
    else if (dons_or_accs == 1)
    {
        for ( int i = 0 ; i < len ; i ++ )  sum += l->A.accs[i];
        for ( int i = 0 ; i < len ; i ++ )  l->A.accs[i] /= sum;
    }

    printf("\t\u2713\n");
}

double safe_log(double x)                                                                       // handle log 0
{
    if ( x == 0.0 )
    {
        const double epsilon = 1e-10;                                                           // log 0 gonna crush in C, name a value for it
        return log(x + epsilon);                                                                // return huge negative log value
    }
    else    return log(x);
}

double log_sum_exp(double *logs, int n) 
{
    double max = logs[0];

    for (int i = 1; i < n ; i++)                                                                // get local maximum
    {
        if( logs[i] > max )     max = logs[i];                  
    }

    double sum = 0.0;

    for( int i = 0 ; i < n ; i++)       sum += exp(logs[i] - max_log);                          // log soft max trick

    return max + log(sum);       
}

void allocate_alpha(Observed_events *info, Forward_algorithm *alpha , Explicit_duration *ed)                            
{
    printf("Start allocate memory for the forward algorithm:");

    int arary_size = info->T - 2 * FLANK;

    /*
        alpha->a[t][i]
        [t]: specific time among all observed events 
        [i]: [0] for exon ; [1] for intron

        each spot is storing a(t)(m, 1) ; based on 2006 implementation
        [m]: types of hidden state
    */

    alpha->a    = malloc ( ( arary_size ) * sizeof(double*) );         
    for (int i = 0 ; i < arary_size; i++ )     
        alpha->a[i] = calloc( HS , sizeof(double) );                                        
    
    /*
        alpha->basis[i][d]
        [i]: [0] for exon ; [1] for intron
        [d]: max duration for exon or intron ; depends on [i]

        each one assign 1D array for each t - 1 layer of all possible D computation
        based on forward algorithm; each at(m, d) partially depends one a(t-1)(m, d+1)
    */

    alpha->basis    = malloc( HS * sizeof(double*) );                   
    alpha->basis[0] = calloc( ed->max_len_exon, sizeof(double) );
    alpha->basis[1] = calloc( ed->max_len_intron, sizeof(double));

    printf("\t\u2713\n");
}

void basis_forward_algorithm(Lambda *l, Explicit_duration *ed,  Forward_algorithm *alpha, Observed_events *info)
{
    /*
        first loop: updating all possible a(1)(m, d)

        a(1)(m, d) = pi(m) * bm(o1) * pm(d)
        pi(m):  initial probability
        bm(o1): emission probability
        pm(d): explicit duration probability

        well we don't need initial prob for us; since uh, exon = 1.0; intron = 0.0;
    */

    int max_len;
    if ( info->T >= ed->max_len_exon)    max_len = ed->max_len_exon;
    else                                 max_len = info->T;

    double log_bm_sum = 0.0;

    for ( int d = 0 ; d < max_len ; d ++)
    { 

        int    index         = base4_to_int(info->numerical_sequence, d - 3 + FLANK, 4);
        double emission_prob = l->B.exon[index];

        log_bm_sum += log(emission_prob);

        /*
            explicit duration is the only term that gonna be 0.0
            if it's 0; means it's invalid length for exon
            directly assign it as 0.0
        */

        double total;
        double ed_prob = ed->exon[d];

        if (ed_prob == 0.0)
        {
            total = 0.0;
            alpha->basis[0][d] = total;
            continue;     
        }

        total = exp(log(ed_prob) + log_bm_sum);
        alpha->basis[0][d] = total;
    }

    /*
        another part of basis function
        update a(1)(m, 1);
        they all 0.0 cuz no intron or exon can be length of 0
        not generalize here; specific for our interests
    */

    alpha->a[0][0] = 0.0;
}

void forward_algorithm(Lambda *l, Forward_algorithm *alpha, Observed_events *info, Explicit_duration *ed, Viterbi_algorithm *vit)
{
    printf("Start computation for forward algorithm.\n");
    printf("\tStart working on total Seq len: %d, Flank size: %d, Real bps len without Flank %d\n", info->T, FLANK, info->T - 2 * FLANK);

    // area for computation//

    int terminal_region_intron = info->T - 2 * FLANK - ed->min_len_exon;

    for ( int t = 1 ; t < info->T - 2 * FLANK ; t ++ )
    {
        /*
            [tau]: the residual/remaining time for explicit duration
        */

        int tau = info->T - 2 * FLANK - t;

        for ( int i = 0 ; i < HS ; i ++ )
        {
            /*
                [bm]:    fancy way of saying emission probability

                [trans]: sum(n != m) a(t - 1)(m, 1) * transition prob
                    all possible transition to current hidden event    
            */

            double log_bm_sum = 0.0;
            double log_trans;
            double trans_prob;

            if ( i == 0 )
            {
                int index = base4_to_int(info->numerical_sequence , t + FLANK - 5, 6);
                trans_prob = l->A.accs[index];
            }
            else
            {
                int index = base4_to_int(info->numerical_sequence , t + FLANK , 5);
                trans_prob = l->A.dons[index];
            }
            
            // the previous state
            int j = (i == 0) ? 1 : 0;

            if (alpha->a[j][t - 1] == 0.0)   log_trans = 0.0;
            else if (trans_prob == 0.0)      log_trans = 0.0;
            else log_trans = log(trans_prob) + log(alpha->a[j][t - 1]);

            for ( int d = 0 ; d < tau ; d ++ )
            {   
                // calculate emission probability and update the sum (log space)
                int index = base4_to_int(info->numerical_sequence, d + FLANK + t - 3);
                double emission_prob = ( i == 0 ) ? l->B.exon[index] : l->B.intron[index];
                log_bm_sum += log(emission_prob);

                //  get explicit duration probability
                double ed_prob = ( i == 0 ) ? ed->exon[d] : ed->intron[d];

                /*
                    calculate final forward component alpha(t)(m, d)
                    formula: a(t - 1)(m, d + 1)( bm(Ot) ) + sum(n != m) a(t - 1)(n, 1) transition(n , 1) * bm(Ot) * pm(d)

                    for our case: only intron and exon
                    final formula: bm(Ot) ( a(t - 1)(m, d + 1) + a(t - 1)(n, 1) transition(n , 1) pm(d))
                */

                double total;

                l->log_values[0] = basis[i][d + 1]
                l->log_values[1] = (ed_prob != 0) exp(log_trans) + ed_prob : 0.0;
                total = log_sum_exp(l->log_values, 2);
                alpha->a[t][i] = exp( log(total) + log_bm_sum );
            }
        }
    }


    printf("\tComputation for forward algorithm finished. \n");
    printf("\n");
}

void free_alpha(Observed_events *info, Forward_algorithm *alpha)
{
    printf("Clearning up forward algorithm memory.\n");

    for (int i = 0 ; i < ( info->T - 2 * FLANK )    ; i ++ )       free( alpha->a[i] );
    free( alpha->a ); 
    printf("\tFinish up clearning alpha memory.\n");

    for (int i = 0 ; i < ( info->T - 2 * FLANK - 1) ; i ++ )       free( alpha->a_star[i] );
    free( alpha->a_star );
    printf("\tFinish up clearning alpha_star memory.\n");
}

void allocate_beta(Observed_events *info, Backward_algorithm *beta)                             // assign data structure for backward algorithm
{
    printf("Start allocate memory for the backward algorithm:\n");

    int array_size = info->T - 2 * FLANK;

    beta->b      = malloc ( ( array_size )    * sizeof(double*) );                                  
    for (int i = 0 ; i < array_size     ; i++ )   beta->b[i]      = calloc( HS , sizeof(double) );
    printf("\tFinish beta memory allocation\n"); 

    beta->b_star = malloc ( ( array_size - 1) * sizeof(double*) );
    for (int i = 0 ; i < array_size - 1 ; i++ )   beta->b_star[i] = calloc( HS , sizeof(double) );
    printf("\tFinish beta_star memory allocation\n"); 
    printf("\tFinished\n");
}

void initial_backward_algorithm(Lambda *l, Backward_algorithm *beta, Observed_events *info, Explicit_duration *ed)
{
    printf("Start initialize backward algorithm:");

    int last_t = info->T - 2 * FLANK - 1;                                                       // minus on cus start as 0
    beta->b[last_t][0] = 1.0;                                                                   // only consider end as an exon
    beta->b[last_t][1] = 0.0;                                                                   // intron ending is not valid

    // basis for exon
    for ( int t = last_t - 1 ; t >= 0 ; t --)
    {
        int log_index = 0; 

        // explicit duration probability
        int duration_len = last_t - t;                                                          // make sence if plug that back to condition for t
        double ed_prob = ed->exon[duration_len];                                            
        l->log_values[log_index++] = safe_log(ed_prob);                                         // send to log space
        
        // get emission prob
        double log_emission_sum = 1.0;

        for ( int s = t ; s <= last_t ; s ++ )
        {
            int index = base4_to_int(info->numerical_sequence, s + FLANK - 3, 4);
            double emission_prob = l->B.exon[index];
            log_emission_sum += safe_log(emission_prob);
        }
        l->log_values[log_index++] = log_emission_sum;

        // calculate up
        beta->b[t][0] = exp( ( log_sum_exp (l->log_values, log_index) ) );  
    }

    printf("\t\u2713\n");
}

void backward_algorithm(Lambda *l, Backward_algorithm *beta, Observed_events *info, Explicit_duration *ed)
{
    printf("Start backward algorihtm.\n");
    printf("\tStart working on total Seq len: %d, Flank size: %d, Real bps len without Flank %d\n", info->T, FLANK, info->T - 2 * FLANK);

    for ( int t = info->T - 2* FLANK - 2 ; t >= 0 ; t-- )                                       // -1 cuz array start at 0; -1 again since already set up last one
    {

        for ( int i = 0 ; i < HS ; i ++ )                                                       // the before position; 0 for exon, 1 for intron
        {
            if ( ( t > info->T - 2 * FLANK - 1 - ed->min_len_exon ) && ( i == 1 ) )
            {
                beta->b[t][i]      = 0.0;
                beta->b_star[t][i] = 0.0;
                continue;
            }

            int log_index = 0;                                                                  // prepare for log-space calcualtion
            int max_len = (i == 0) ? ed->max_len_intron : ed->max_len_exon ;                    // we are checking 

            for ( int j = 0 ; j < HS ; j ++ )                                                   // the situation of transition out
            {
                if (i == j) continue;                                                           // no self-transition for HSMM
                
                if (beta->b_star[t][j] == 0.0)                                                  // if that is not pre-calculated
                {
                    for ( int d = 1 ; d <= max_len; d ++ )                                          
                    {
                        if ( (t + d)  > ( info->T - 2 * FLANK - 1 ) )    break;                     // can not reach further than last beta 

                        double log_emission_sum = 1.0;                                              // what we want to calculate out here

                        for ( int s = t + 1 ; s <= t + d ;  s++ )
                        {
                            int index = base4_to_int(info->numerical_sequence, s - 3 + FLANK, 4);   // how we get index value from 4 base pair
                            double emission_prob = ( j == 0) ? l->B.exon[index] : l->B.intron[index]; 
                            log_emission_sum += safe_log(emission_prob);                        // update value
                        }

                        // get explicit duration prob
                        double ed_prob = ( i == 0) ? ed->intron[d - 1] : ed->exon[d - 1] ;          // different logic as forward

                        // formal computation
                        double all = safe_log( beta->b[t + d][j] ) + safe_log( ed_prob ) + log_emission_sum;
                        l->log_values[log_index++] = all;
                    }

                    beta->b_star[t][j] = exp( log_sum_exp (l->log_values, log_index) );             // update the value         
                }
            }

            // get transition prob
            double trans_prob;                                                                       // transition prob

            if ( i == 0)
            {
                int index = base4_to_int(info->numerical_sequence , t + FLANK + 4, 5);               // 5bps motif ; donor site transition prob
                trans_prob = l->A.dons[index];
            }
            else
            {
                int index = base4_to_int(info->numerical_sequence , t + FLANK - 5 , 6);              // 6bps motif ; acceptor site transiton prob
                trans_prob = l->A.accs[index];
            }

            l->log_values[log_index ++] = safe_log(trans_prob);
            beta->b[t][i] += exp( ( log_sum_exp (l->log_values, log_index) ) );   
        }
    }
    printf("\tComputation for backward algorithm finished.\n");
}

void free_beta(Observed_events *info, Backward_algorithm *beta)
{
    printf("Clearning up backward algorithm memory.\n");

    for (int i = 0 ; i < (info->T - 2 * FLANK )    ; i ++)     free( beta->b[i] );
    free( beta->b );
    printf("\tFinish up clearning beta memory.\n");

    for (int i = 0 ; i < (info->T - 2 * FLANK - 1) ; i++)      free( beta->b_star[i]);
    free( beta->b_star);
    printf("\tFinish up clearning beta_star memory.\n");
}

// output //
