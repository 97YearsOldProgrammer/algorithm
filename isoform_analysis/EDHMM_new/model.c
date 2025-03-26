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

double log_sum_exp(double *array, int n) 
{
    double max = array[0];

    for (int i = 1; i < n ; i++)
    {
        if( array[i] > max )     max = array[i];                  
    }

    double sum = 0.0;

    for( int i = 0 ; i < n ; i++)
    {
        sum += exp(array[i] - max);
    }

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

    double log_bm_sum_exon = 0.0;

    for ( int d = 0 ; d < max_len ; d ++)
    { 

        int    index         = base4_to_int(info->numerical_sequence, d - 3 + FLANK, 4);
        double emission_prob = l->B.exon[index];

        log_bm_sum_exon += log(emission_prob);

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

        total = exp(log(ed_prob) + log_bm_sum_exon);
        alpha->basis[0][d] = total;
    }
    
    alpha->a[0][0] = 0.0;

    /*
        this part is assign prob for intron
        since for this model; it;s kinda different
        we don't have initial intron probability
        btw the real intron initial probability comes after min len of exon
        which make sense
        so we make our first layer of calculation at point where passes min exon len

        [tau] = refers to remaining time for explicit duration
    */

    for ( int t = 0 ; t < ed->min_len_exon ; t ++ )
    {
        alpha->a[1][t] = 0.0;
    }

    int tau_intron = info->T - 2 * FLANK - ed->min_len_exon;
    double log_sum_intron = 0.0;

    for ( int d = 0 ; d < tau_intron ; d ++ )
    {
        int index = base4_to_int(info->numerical_sequence, d - 3 + FLANK + ed->min_len_exon);
        double emission_prob = l->B.intron[index];

        log_sum_intron += log(emission_prob);

        double total;
        double ed_prob = ed->intron[d];

        if (ed_prob == 0.0)
        {
            total = 0.0;
            alpha->basis[1][d] = total;
            continue;
        }

        total = exp(log(ed_prob) + log_bm_sum_exon);
        alpha->basis[0][d] = total;
    }
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

                [trans]: transition prob  

                [log_trans]: sum(n != m) a(t - 1)(m, 1) * transition prob
                    all possible transition to current hidden state
            */

            int initial_intron_len = ed->min_len_exon ;

            if ( t <= initial_intron_len && i == 1)
            {
                alpha->a[t][i] = 0.0;
                continue;
            }

            // the transition prob
            double trans_prob;

            // get the transition probability at this point
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
            
            double log_trans;

            // the previous state
            int j = (i == 0) ? 1 : 0;

            if (alpha->a[j][t - 1] == 0.0)   log_trans = 0.0;
            else if (trans_prob == 0.0)      log_trans = 0.0;
            else log_trans = log(trans_prob) + log(alpha->a[j][t - 1]);

            if ( i == 1 )   tau -= ed->min_len_exon;

            double log_bm = 0.0;

            for ( int d = 0 ; d < tau ; d ++ )
            {   
                // calculate emission probability and update the sum (log space)
                int index = base4_to_int(info->numerical_sequence, d + t - 3, 4 );

                double emission_prob = ( i == 0 ) ? l->B.exon[index] : l->B.intron[index];
                
                log_bm += log(emission_prob);

                //  get explicit duration probability
                double ed_prob = ( i == 0 ) ? ed->exon[d] : ed->intron[d];

                /*
                    calculate final forward component alpha(t)(m, d)
                    formula: a(t - 1)(m, d + 1)( bm(Ot) ) + sum(n != m) a(t - 1)(n, 1) transition(n , 1) * bm(Ot) * pm(d)

                    for our case: only intron and exon
                    final formula: bm(Ot) ( a(t - 1)(m, d + 1) + a(t - 1)(n, 1) transition(n , 1) pm(d))
                */

                double total;

                l->log_values[0] = alpha->basis[t - 1][d + 1];

                if (ed_prob == 0.0)     l->log_values[1] = 0.0;
                else l->log_values[1]   = exp ( log_trans + log(ed_prob) );

                total = log_sum_exp(l->log_values, 2);
                alpha->a[t][i] = exp( log(total) + log_bm );
            }
        }
    }

    printf("\tComputation for forward algorithm finished. \n");
    printf("\n");
}

void viterbi_basis(Viterbi_algorithm *vit, Forward_algorithm *alpha)
{
    double sum = 0.0;
    
    gamma.exon   = log_sum_exp(alpha->basis[0], 2);
    gamma.intron = log_sum_exp(alpha->basis[1], 2);
}

void free_alpha(Observed_events *info, Forward_algorithm *alpha, Explicit_duration *ed)
{
    printf("Clearing up forward algorithm memory:");
    
    int array_size = info->T - 2 * FLANK;
    
    for (int i = 0; i < array_size; i++) {
        free(alpha->a[i]);
    }

    free(alpha->a);
    free(alpha->basis[0]);
    free(alpha->basis[1]);
    free(alpha->basis);

    print("\tFinished\n");
}

void allocate_beta(Observed_events *info, Backward_algorithm *beta, Explicit_duration *ed)                             
{
    printf("Start allocate memory for the backward algorithm:\n");

    int arary_size = info->T - 2 * FLANK;
                                    
    /*
        β->basis[i][d]
        [i]: [0] for exon ; [1] for intron
        [d]: max duration for exon or intron ; depends on [i]

        each one assign 1D array for each t - 1 layer of all possible D computation
    */

    beta->basis    = malloc( HS * sizeof(double*) );                   
    beta->basis[0] = calloc( ed->max_len_exon, sizeof(double) );
    beta->basis[1] = calloc( ed->max_len_intron, sizeof(double));

    printf("\tFinished\n");
}

void initial_backward_algorithm(Lambda *l, Backward_algorithm *beta, Observed_events *info, Explicit_duration *ed)
{
    printf("Start initialize backward algorithm:");

    for ( int i = 0 ; i < HS ; i++ )
    {
        beta->basis[i][0] = 1.0;
    }

    printf("\tFinished\n");
}

void initial_viterbi_algorithm(Viterbi_algorithm *vit, Observed_events *info)
{
    printf("Start Initialize Viterbi Algorithm");

    /*
        given recursive viterbi formula: γ(t) = y(t+1)(m) + sum(n != m) ( ξ(t+1)(m, n) - ξ(t+1)(n, n) )
        for our model: degrade that formula
        γ(t) = y(t+1)(m) + ξ(t+1)(m, n) - ξ(t+1)(n, n)

        given definition of ξ
        ξ(t)(m, n) = α(t - 1)(m, 1) * amn * bn(Ot) * sum(d >= 1) (pn(d) β(t)(n, d))

        in terms of t + 1
        ξ(t+1)(m, n) = α(t)(m, 1) * amn * bn(Ot + 1) * sum(d >= 1) (pn(d) β(t+1)(n, d))
    */

    vit->xi = malloc( HS * sizeof(double) );
    vit->gamma = malloc( HS * sizeof(double) );
    vit->path = malloc( info->T * sizeof(int) );

    printf("\tFinished\n");
}

void argmax_viterbi(Viterbi_algorithm *vit)
{   
    int argmax;
    /*
        for our model: consider following viterbi formula

        γ(t)(exon)   = γ(t+1)(exon)   + ξ(exon, intron) - ξ(intron, exon)
        γ(t)(intron) = γ(t+1)(intron) + ξ(intron, exon) - ξ(exon, intron)
        
        deduct from 2006 paper
    */
   

}
void xi_calculation(Lambda *l, Forward_algorithm *alpha, Viterbi_algorithm *vit, Observed_events *info, double backward_sum, int t)
{
    /*
        [alpha_component]: α(t - 1)(m, 1)
        [trans_prob]: either donor or acceptor
    */

    for( int i = 0 ; i < HS ; i++ )
    {
        double alpha_component;
        alpha_component = alpha->a[t][i];

        double trans_prob;
        if ( i == 0 )
        {
            int index = base4_to_int(info->numerical_sequence, t , 5);
            trans_prob = l->A.dons[index];
        }
        else
        {
            int index = base4_to_int(info->numerical_sequence, t - 6, 6);
            trans_prob = l->A.accs[index];
        }
        
        double emission_prob;
        int index_emission = base4_to_int(info->numerical_sequence, t - 3, 4);
        emission_prob = (i == 0) ? l->B.intron[index_emission] : l->B.exon[index_emission];

        double total;

        if      (trans_prob == 0.0)         total = 0.0;
        else if (alpha_component == 0.0)    total = 0.0;
        else    total = exp( log(trans_prob) + log(alpha_component) + log(emission_prob) );

        vit->xi[i] = total;
    }

}

void backward_algorithm(Lambda *l, Backward_algorithm *beta, Observed_events *info, Explicit_duration *ed)
{
    printf("Start Backward Algorithm:");

    int start_bps = info->T - 2 * FLANK - 2;

    for ( int t = start_bps ; t >= 0 ; t-- )                                      // -1 cuz array start at 0; -1 again since already set up last one
    {
        int tau = info->T - 2 * FLANK - start_bps;

        for ( int i = 0 ; i < HS ; i ++ )                                                       // the before position; 0 for exon, 1 for intron
        {
            if ( t > start_bps - ed->min_len_exon + 1 && i == 1) 
            {
                double backward_sum = 0.0;

                continue;
            }

            /*
                calculate β(t)(m, 1) = sum(n != m) * a(mn) * bn(Ot+1) * sum(d>=1) pn(d) * beta(t+1)(n, d)
                the sum(d>=1) pn(d) * beta(t+1)(n, d) can be coupled with ξ calculation for viterbi algorithm

                [total]: β(t)(m, 1)
            */

            double total;

            double backward_sum;

            for ( int d = 0 ; d < tau ; d++ )
            {
                double ed_prob = ( i == 0 ) ? ed->exon[d] : ed->intron[d];
                
                if (ed_prob == 0.0)     l->log_values[d] = 0.0;
                else l->log_values[d] = exp( log( beta->basis[t + 1][d] ) + log (ed_prob) );
            }

            backward_sum = log_sum_exp(l->log_values, tau);

            int j = ( i == 0) ? 1 : 0;

            double trans_prob;

            if ( i == 0 )
            {   
                int index = base4_to_int(info->numerical_sequence, t + 1 , 5);
                trans_prob = l->A.dons[index];
            }
            else if ( i == 1 )
            {
                int index = base4_to_int(info->numerical_sequence, t - 5 , 6);
                trans_prob = l->A.accs[index];
            }

            double emission_prob;

            int index = base4_to_int(info->numerical_sequence, t - 2, 4);
            emission_prob = (i == 0) ? l->B.intron[index] : l->B.exon[index];

            if   (trans_prob == 0.0)   total = 0.0;
            else total = exp( log(trans_prob) + log(emission_prob) + log(backward_sum) );

            /*
                calculate β(t)(m, d) = bm(Ot+1) * β(t+1)(m, d - 1)
                    for all possible d > 1
            */
            
            int index = base4_to_int(info->numerical_sequence, t - 2, 4);
            double emission_prob = (i == 0) ? l->B.exon[index] : l->B.intron[index];

            for( int d = 1 ; d < tau ; d++ )
            {   
                double previous_node;
                previous_node = beta->basis[i][d - 1];

                if (previous_node == 0.0)   beta->basis[i][d] = 0.0;
                else    beta->basis[i][d] = exp( log(emission_prob) + log(previous_node) );
            }

            beta->basis[i][0] = total;    
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
