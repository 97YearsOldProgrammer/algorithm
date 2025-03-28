#include "model.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

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
    printf("\tFinished\n");
    printf("\n");
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
    double value = 0.0;

    for (int i = 0 ; i < length ; i ++)
    {
        if (array[i] == 0.0)
        {
            value = 0.0;
            break;
        }

        value += log(array[i]);
    }

    if (value != 0.0)   value = exp(value);

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

        double value;
        for ( int i = 0 ; i < len ; i ++ )
        {
            value = l->A.dons[i];
            if (value == 0.0)       continue;
            else                    l->A.dons[i] = exp ( log(value) - log(sum) );
            
        }
    }
    else if (dons_or_accs == 1)
    {
        for ( int i = 0 ; i < len ; i ++ )  sum += l->A.accs[i];

        double value;
        for ( int i = 0 ; i < len ; i ++ )
        {
            value = l->A.accs[i];

            if (value == 0.0)       continue;
            else                    l->A.accs[i] = exp( (log(value) - log(sum) ) );
        }
    }

    printf("\tFinished\n");
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

    printf("\tFinished\n");
}

void basis_forward_algorithm(Lambda *l, Explicit_duration *ed,  Forward_algorithm *alpha, Observed_events *info)
{
    printf("Start forward algorithm basis calculation:");

    /*
        first loop: updating all possible a(1)(m, d)

        a(1)(m, d) = pi(m) * bm(o1) * pm(d)
        pi(m):  initial probability
        bm(o1): emission probability
        pm(d): explicit duration probability
        [tau]: used for stimulate different duration

        well we don't need initial prob for us; since uh, exon = 1.0; intron = 0.0;
    */

    int tau_exon;
    int comp_len_exon = info->T - 2 * FLANK;

    if  (comp_len_exon > ed->max_len_exon)      tau_exon = ed->max_len_exon;
    else                                        tau_exon = comp_len_exon;

    int    index         = base4_to_int(info->numerical_sequence, FLANK - 3, 4);
    double emission_prob = l->B.exon[index];

    /*
        explicit duration is the only term that gonna be 0.0
        if it's 0; means it's invalid length for exon
        directly assign it as 0.0
    */

    double total;
    double ed_prob;

    for( int d = 0 ; d < tau_exon ; d ++)
    {
        ed_prob = ed->exon[d];

        if (ed_prob == 0.0)
        {
            alpha->basis[0][d] = 0.0;
        }
        else
        {
            total = exp( log(ed_prob) + log(emission_prob) );
            alpha->basis[0][d] = total;
        }
    }

    /*
        this part is assign prob for intron
        since for this model; it;s kinda different
        we don't have initial intron probability
        btw the real intron initial probability comes after min len of exon
        which make sense
        so we make our first layer of calculation at point where passes min exon len
    */

    int tau_intron;
    int comp_len_intron = comp_len_exon - ed->min_len_exon * 2;

    if ( comp_len_intron > ed->max_len_intron)      tau_intron = ed->max_len_intron;
    else                                            tau_intron = comp_len_intron;

    index = base4_to_int(info->numerical_sequence, FLANK + ed->min_len_exon - 3 , 4);
    emission_prob = l->B.intron[index];

    for( int d = 0 ; d < tau_intron ; d ++)
    {
        ed_prob = ed->intron[d];

        if (ed_prob == 0.0)
        {
            alpha->basis[1][d] = 0.0;
        }
        else
        {
            total = exp( log(ed_prob) + log(emission_prob) );
            alpha->basis[1][d] = total;
        }
    }

    printf("\tFinished\n");
}

void forward_algorithm(Lambda *l, Forward_algorithm *alpha, Observed_events *info, Explicit_duration *ed)
{
    printf("Start computation for forward algorithm:");

    /*
        recall
        [tau]: the residual/remaining time for explicit duration
    */

    int tau = info->T - 2 * FLANK;

    for ( int t = 1 ; t < info->T - 2 * FLANK ; t ++ )
    {
        tau --;

        for ( int i = 0 ; i < HS ; i ++ )
        {
            int tau_modified;
            if ( i == 1)    tau_modified = tau - ed->min_len_exon;
            else            tau_modified = tau;

            /*
                boundary condition

                before first min exon exit; no intron possible to generate
                update only alpha[t][i]
                which is all time interval before min_len_exon for α(intron, 1)

                and also for intron near the end which min_len_exon exist
            */

            if ( ( tau_modified <= 0 || t > info->T - ed->min_len_exon ) && i == 1)
            {
                alpha->a[t][i] = 0.0;
                continue;
            }

            /*
                [bm]:    fancy way of saying emission probability

                [trans]: transition prob  

                [alpha_trans]: sum(n != m) a(t - 1)(m, 1) * transition prob
                    all possible transition to current hidden state
            */

            double trans_prob;

            if ( i == 0 )
            {
                int index = base4_to_int(info->numerical_sequence , t + FLANK - 6, 6);
                trans_prob = l->A.accs[index];
            }
            else
            {
                int index = base4_to_int(info->numerical_sequence , t + FLANK , 5);
                trans_prob = l->A.dons[index];
            }
            
            double alpha_trans;
            double previous_node;
            int j = (i == 0) ? 1 : 0;

            previous_node = alpha->a[t - 1][j];

            if      (previous_node == 0.0)        alpha_trans = 0.0;
            else if (trans_prob == 0.0)           alpha_trans = 0.0;
            else                                  alpha_trans = exp( log(trans_prob) + log(previous_node) );

            double emission_prob;

            int index = base4_to_int(info->numerical_sequence, t + FLANK - 3, 4);
            emission_prob = ( i == 0 ) ? l->B.exon[index] : l->B.intron[index];

            double ed_prob;
            double total;  

            /*
                since mostly explicit duration when d = 0 is 0
                α(t)(m, 1) = α(t - 1)(m, 2) * bm(ot)
            */
            
            previous_node = alpha->basis[i][1];

            if ( previous_node == 0.0 )     total = 0.0;
            else                            total = exp( log( previous_node ) + log(emission_prob) );

            alpha->a[t][i]     = total;
            alpha->basis[i][0] = total;

            for ( int d = 1 ; d < tau_modified ; d ++ )
            {   
                /*
                    calculate final forward component alpha(t)(m, d)
                    formula: a(t - 1)(m, d + 1)( bm(Ot) ) + sum(n != m) a(t - 1)(n, 1) transition(n , 1) * bm(Ot) * pm(d)

                    for our case: only intron and exon
                    final formula: bm(Ot) ( a(t - 1)(m, d + 1) + a(t - 1)(n, 1) transition(n , 1) pm(d))
                */

                previous_node = alpha->basis[i][d + 1];
                if ( previous_node == 0.0 )     l->log_values[0] = 0.0;
                else                            l->log_values[0] = previous_node;

                ed_prob = ( i == 0 ) ? ed->exon[d] : ed->intron[d];
                if   (ed_prob == 0.0)   l->log_values[1] = 0.0;
                else                    l->log_values[1] = exp ( log(alpha_trans) + log(ed_prob) );

                total = log_sum_exp(l->log_values, 2);
                alpha->basis[i][d] = exp( log(total) + log(emission_prob) );
            }
        }
    }

    printf("\tFinished\n");
}

void free_alpha(Observed_events *info, Forward_algorithm *alpha)
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

    printf("\tFinished\n");
}

void allocate_viterbi(Viterbi_algorithm *vit, Observed_events *info)
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

    vit->xi    = calloc( HS , sizeof(double) );
    vit->gamma = malloc( HS * sizeof(double) );
    vit->path  = malloc( (info->T - 2 * FLANK) * sizeof(int) );

    printf("\tFinished\n");
}

void argmax_viterbi(Viterbi_algorithm *vit, int t)
{   
    int argmax;

    /*
        for our model: consider following viterbi formula

        γ(t)(exon)   = γ(t+1)(exon)   + ξ(exon, intron) - ξ(intron, exon)
        γ(t)(intron) = γ(t+1)(intron) + ξ(intron, exon) - ξ(exon, intron)
        
        deduct from 2006 paper
    */

    vit->gamma[0] += vit->xi[0] - vit->xi[1];
    vit->gamma[1] += vit->xi[1] - vit->xi[0];

    /*
        not set for equal conditon here; i don't believe that would happen
        or maybe it would happen; lets see
    */

    if ( vit->gamma[0] > vit->gamma[1] )    argmax = 0;
    else                                    argmax = 1;

    vit->path[t] = argmax;
}

void xi_calculation(Lambda *l, Forward_algorithm *alpha, Viterbi_algorithm *vit, Observed_events *info, double backward_sum, int t, int type)
{
    assert(type == 0 || type == 1);

    /*
        input parameter

        [backward_sum]: sum(d >= 1) pn(d)*β(n, d)
            that was compute for ξ(m, n)
        [t]: time when it computed
            which means which t - 1 gamma we wanna compute
        [type]: either 0 or 1
            so we know which ξ(m, n) we shall update in vit

        computation notation

        [alpha_component]: α(t - 1)(m, 1)
        [trans_prob]: either donor or acceptor
    */


    double alpha_component;
    alpha_component = alpha->a[t][type];

    double trans_prob;
    if (type == 0)
    {
        int index = base4_to_int(info->numerical_sequence, t + 1 + FLANK, 5);
        trans_prob = l->A.dons[index];
    }
    else
    {
        int index = base4_to_int(info->numerical_sequence, t - 6 + FLANK, 6);
        trans_prob = l->A.accs[index];
    }
        
    double emission_prob;
    int index_emission = base4_to_int(info->numerical_sequence, t - 3 + FLANK, 4);
    emission_prob = (type == 0) ? l->B.intron[index_emission] : l->B.exon[index_emission];

    double total;
    if      (trans_prob == 0.0)         total = 0.0;
    else if (alpha_component == 0.0)    total = 0.0;
    else    total = exp( log(trans_prob) + log(alpha_component) + log(emission_prob) + log(backward_sum) );

    vit->xi[type] = total;
}

void viterbi_basis(Viterbi_algorithm *vit, Forward_algorithm *alpha)
{
    printf("Start assign basis for Viterbi Algorithm:");

    /*
        γ(t)(m) = sum(d>=1) α(t)(m, d)
        at final t; there is only one α(t)(m, 1) existed
        based on definition for d
        (if forward algorithm is compute without boundary issue)

        γ(t)(1) = shall be 0; intron uh is invalid in that position
    */

    double gamma_exon;
    double gamma_intron;

    gamma_exon   = alpha->basis[0][0];
    gamma_intron = alpha->basis[1][0];

    vit->gamma[0] = gamma_exon;
    vit->gamma[1] = gamma_intron;
    
    printf("\tFinished\n");
}

void allocate_beta(Backward_algorithm *beta, Explicit_duration *ed)                             
{
    printf("Start allocate memory for the backward algorithm:");
                                    
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

void initial_backward_algorithm(Backward_algorithm *beta)
{
    printf("Start initialize backward algorithm:");

    for ( int i = 0 ; i < HS ; i++ )
    {
        beta->basis[i][0] = 1.0;
    }

    printf("\tFinished\n");
}

void backward_algorithm(Lambda *l, Backward_algorithm *beta, Observed_events *info, Explicit_duration *ed, Viterbi_algorithm *vit, Forward_algorithm *alpha)
{
    printf("Start Backward Algorithm:");

    int start_bps = info->T - 2 * FLANK - 2;
    int tau = 0;

    for ( int t = start_bps ; t >= 0 ; t-- )                                      // -1 cuz array start at 0; -1 again since already set up last one
    {
        argmax_viterbi(vit, t + 1);

        tau ++;

        for ( int i = 0 ; i < HS ; i ++ )                                                       // the before position; 0 for exon, 1 for intron
        {
            int tau_modified;

            if ( i == 0 )   tau_modified = tau;
            else            tau_modified = tau - ed->min_len_exon;

            if       ( i == 0 && tau >= ed->max_len_exon)                        tau_modified = ed->max_len_exon;
            else if  ( i == 1 && tau >= ed->max_len_intron - ed->min_len_exon)   tau_modified = ed->max_len_intron;
    
            /*
                boundary condition

                [final bps]: total_bps - 2 * FLANK (- 1) as in array form
                [last intron range]: final bps - min_exon (-1) in array form

                [first intron range]: FLANK + min_exon

                in terms of backward algorithm; purpose is calculate summation of all possible backward transition
                with in those boundary; directly assign such calculation to 0; represent they are inreasonable nodes for network
            */

            if ( tau_modified <= 0 && i == 1) 
            {
                double backward_sum = 0.0;
                xi_calculation(l, alpha, vit, info, backward_sum, t, 0);
                continue;
            }

            else if ( t < FLANK + ed->min_len_exon && i == 1)
            {
                double backward_sum = 0.0;
                xi_calculation(l, alpha, vit, info, backward_sum, t, 0);
                continue;
            }

            /*
                calculate β(t)(m, 1) = sum(n != m) * a(mn) * bn(Ot+1) * sum(d>=1) pn(d) * beta(t+1)(n, d)
                the sum(d>=1) pn(d) * beta(t+1)(n, d) can be coupled with ξ calculation for viterbi algorithm

                [total]: β(t)(m, 1)
                [j]: conjudated hidden state
                    if exon(0) it's intron(1); vice versa
            */

            int j = ( i == 0) ? 1 : 0;

            double backward_sum;

            for ( int d = 0 ; d < tau_modified ; d++ )
            {
                double ed_prob = ( j == 0 ) ? ed->exon[d] : ed->intron[d];
                
                if (ed_prob == 0.0)     l->log_values[d] = 0.0;
                else l->log_values[d] = exp( log( beta->basis[j][d] ) + log (ed_prob) );
            }

            backward_sum = log_sum_exp(l->log_values, tau_modified);
            xi_calculation(l, alpha, vit, info, backward_sum, t, j);

            double trans_prob;

            if ( i == 0 )
            {   
                int index = base4_to_int(info->numerical_sequence, t + 1 + FLANK , 5);
                trans_prob = l->A.dons[index];
            }
            else
            {
                int index = base4_to_int(info->numerical_sequence, t - 6 + FLANK , 6);
                trans_prob = l->A.accs[index];
            }

            int index = base4_to_int(info->numerical_sequence, t - 2 + FLANK, 4);
            double emission_prob = (j == 0) ? l->B.exon[index] : l->B.intron[index];

            double total;

            if   (trans_prob == 0.0)   total = 0.0;
            if   (backward_sum == 0.0) total = 0.0;
            else total = exp( log(trans_prob) + log(emission_prob) + log(backward_sum) );

            /*
                calculate β(t)(m, d) = bm(Ot+1) * β(t+1)(m, d - 1)
                    for all possible d > 1
            */
            
            index = base4_to_int(info->numerical_sequence, t - 2 + FLANK, 4);
            emission_prob = (i == 0) ? l->B.exon[index] : l->B.intron[index];

            /*
                logic here is linear array for all layer of network
                and each time passed there would be one node lost for each layer
                every computation based on the previous node
                so we cannot directly update that after calculation
            */

            double waitlist_node;
            waitlist_node = total;

            for( int d = 1 ; d < tau_modified ; d++ )
            {   
                double previous_node;
                previous_node = beta->basis[i][d - 1];

                double  update_node;
                if      (previous_node == 0.0)  update_node = 0.0;
                else    update_node = exp( log(emission_prob) + log(previous_node) );

                beta->basis[i][d - 1] = waitlist_node;
                waitlist_node = update_node;
            } 
        }
    }

    argmax_viterbi(vit, 0);

    printf("\tFinished.\n");
}

void free_beta(Backward_algorithm *beta)
{
    printf("Clearning up backward algorithm memory:");

    free(beta->basis[0]);
    free(beta->basis[1]);
    free(beta->basis);

    printf("\tFinished\n");
}

void free_viterbi(Viterbi_algorithm *vit)
{
    printf("Clearning up viterbi algorithm memory:");

    free(vit->path);
    free(vit->gamma);
    free(vit->xi);

    printf("\tFinished\n");
}

void viterbi_path_test(Viterbi_algorithm *vit, Observed_events *info, Explicit_duration *ed)
{
    printf("\n");
    printf("Start Viterbi Check:");

    int state = vit->path[0];

    if ( state == 0 )   printf("Exon");
    else                printf("Intron");

    int bps = 1 + FLANK;

    printf("\t%i", bps);

    for ( int i = 1 ; i < info->T - 2 * FLANK ; i++ )
    {
        bps ++;

        if  (vit->path[i] == state)
        {
            state = vit->path[i];
            continue;
        }
        else
        {
            printf("\t%d\n", bps - 1);
            state = vit->path[i];
            if ( state == 0)    printf("Exon\t%i", bps);
            else                printf("Intron\t%i", bps);
        }
    }
    printf("\t%i\n", bps);
}

