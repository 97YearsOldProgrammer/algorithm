#include "model.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void numerical_transcription(Observed_events *info, const char *seq)
{
    printf("Start transforming original sequence into base4:");

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

    printf("\t\u2713\n");
    printf("We get numerical sequence with Seq len: %d\n\n", len);
}

void setup_initial_probability(Lambda *l)                               // actually no longer needed
{
    printf("Start getting initial probability down");
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

void allocate_alpha(Observed_events *info, Forward_algorithm *alpha)    // assign data structure for alpha
{
    printf("Start allocate memory for the forward algorithm:");
    alpha->a = malloc ( (info->T - 2 * FLANK + 1) * sizeof(double*) );                    // 2x2; outer layer as length of sequence

    for (int i = 0 ; i < (info->T - 2 * FLANK + 1) ; i++ )
    {
        alpha->a[i] = calloc( HS , sizeof(double) );                    // inner layer as number of Hidden States
    }
    printf("\t\u2713\n");
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
    printf("Start initialize forward algorithm:");
   alpha->a[0][0] = 1.0;
   alpha->a[0][1] = 0.0;
    printf("\t\u2713\n");
}

void forward_algorithm(Lambda *l, Forward_algorithm *alpha, Observed_events *info, Explicit_duration *ed)
{
    int terminal_region_intron = info->T - ed->min_len_exon - FLANK + 1;                          // no intron continue beyond this point
    int terminal_region_exon   = info->T - ed->min_len_exon - ed->min_len_intron - FLANK + 1;     // no exon to intron beyond this point 

    if ( terminal_region_exon   <= 0 )  printf("Invalid data input. Check FLANK size. Seq len. Min Exon\n");
    else                                printf("Start forwrad algorithm with terminal intron region %d\n", terminal_region_intron);
    if ( terminal_region_intron <= 0 )  printf("Invalid data input. Check FLANK size. Seq len. Min Intron & Exon\n");   
    else                                printf("Start forward algorithm with terminal exon region %d\n", terminal_region_exon);

    int real_bps = info->T - 2 * FLANK + 1;
    printf("Start working on total Seq len: %d, Flank size: %d, Real bps len without Flank %d\n", info->T + 1, FLANK, real_bps);

    for ( int t = 0 ; t < info->T - 2 * FLANK + 1 ; t ++ )                                      // iterate every element in the time scale
    {
         
        for ( int i = 0 ; i < HS ; i ++ )                                                       // the next position; 0 for exon, 1 for intron
        {

            // biological restrain
            if ( t >= terminal_region_intron && i == 1)                                         // can't have new intron after this
            {                                                                                   // cuz we have to end as exon
                alpha->a[t][i] = 0.0;
                continue;
            }

            int log_index = 0;                                                                  // prepare for log-space calcualtion
            int max_len;

            if      ( i == 0 )    max_len = ed->max_len_exon;                                   // max_len for exon state
            else if ( i == 1 )    max_len = ed->max_len_intron;                                 // max_len for intron state

            for ( int j = 0 ; j < HS ; j ++ )                                                   // the situation of transition out
            {
                if (i == j) continue;                                                           // no self-transition for HSMM

                // biological restrain
                if (t >= terminal_region_exon && j == 0 && i == 1) continue;                    // no exon->intron transition beyond this point
                
                for ( int d = 1 ; d <= max_len ; d ++ )                                         // skip those unsatisfied len; cuz they all 0
                {
                    if ( t < d )   break;                                                       // not making sence if d larger than s

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

                    if      ( all > 0 ) l->log_values[log_index++] = all;
                    else                l->log_values[log_index++] = safe_log(all); 
                }
            }

            for ( int d = 1 ; d <= max_len; d ++ )                                              // for continue probability
            {
                if ( t < d )    break;   

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

                if ( all > 0 ) l->log_values[log_index++] = all;
                else           l->log_values[log_index++] = safe_log(all);
            }

            if ( log_index == 0 )   alpha->a[t][i] = 0.0;
            else                    alpha->a[t][i] = exp( log_sum_exp (l->log_values, log_index) );
        }
    }
}

void allocate_beta(Observed_events *info, Backward_algorithm *beta)                             // assign data structure for backward algorithm
{
    printf("Start allocate memory for the backward algorithm:");

    beta->b = malloc ( (info->T - 2 * FLANK + 1) * sizeof(double*) );                           // 2x2; outer layer as length of sequence

    for (int i = 0 ; i < (info->T - 2 * FLANK + 1) ; i++ )
    {
        beta->b[i] = calloc( HS , sizeof(double) );                                             // inner layer as number of Hidden States
    }
    printf("\t\u2713\n");
}

void initial_backward_algorithm(Lambda *l, Backward_algorithm *beta, Observed_events *info)
{
    printf("Start initialize backward algorithm:");
    int last_bps = info->T - 2 * FLANK; 

    // constrain ; only possible to end as an exon 
    beta->b[last_bps][0] = 1.0;
    beta->b[last_bps][1] = 0.0;
    printf("\t\u2713\n");
}

void backward_algorithm(Lambda *l, Backward_algorithm *beta, Observed_events *info, Explicit_duration *ed)
{
    int start_region_intron = info->T - ed->min_len_exon - FLANK + 1;

    if ( start_region_intron < FLANK || start_region_intron < 0 )   printf("Invalid data. Check Seq_len. FLANK. Min_exon.\n");
    else                                                            printf("Start backward algorithm. No exon to intron transiton before %d.\n", start_region_intron);

    for ( int t = info->T - 2* FLANK - 1 ; t >= 0 ; t-- )
    {

        for ( int i = 0 ; i < HS ; i ++ )                                                       // the before position; 0 for exon, 1 for intron
        {

            // biological restrain
            if ( t >= start_region_intron && i == 1)                                            // can't have new intron after this
            {                                                                                   // cuz we have to end as exon
                beta->b[t][i] = 0.0;
                continue;
            }

            int log_index = 0;                                                                  // prepare for log-space calcualtion
            int max_len;

            if      ( i == 0 )    max_len = ed->max_len_exon;                                   // max_len for exon state
            else if ( i == 1 )    max_len = ed->max_len_intron;                                 // max_len for intron state

            for ( int j = 0 ; j < HS ; j ++ )                                                   // the situation of transition out
            {
                if (i == j) continue;                                                           // no self-transition for HSMM
                
                for ( int d = 1 ; d <= max_len; d ++ )                                          // skip those unsatisfied len; cuz they all 0
                {
                    if ( (t + d)  > (info->T - 2 * FLANK) )    break;                           // can durate longer than current point       

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
                        int index = base4_to_int(info->numerical_sequence , t - 4 + d + FLANK , 5);     // 5bps motif
                        trans_prob = l->A.dons[index];
                    }
                    else if ( j == 1 && i == 0)                                                 // from intron to exon
                    {
                        int index = base4_to_int(info->numerical_sequence , t + d + FLANK, 6);          // 6bps motif
                        trans_prob = l->A.accs[index];
                    }

                    // formal computation

                    double all = beta->b[t + d][j] * trans_prob * ed_prob * emission_product;
                    if  ( all > 0 ) l->log_values[log_index++] = all;
                    else            l->log_values[log_index++] = safe_log(all);
                }
            }

            for ( int d = 1 ; d <= max_len && t + d >= 0 ; d ++ )                               // for continue probability
            {
                if ( t + d >= info->T - 2 * FLANK + 1)  break;
                
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

                if ( all > 0 ) l->log_values[log_index++] = all;
                else           l->log_values[log_index++] = safe_log(all);
            }

            if (log_index == 0) beta->b[t][i] = 0.0;
            else                beta->b[t][i] = exp( log_sum_exp (l->log_values, log_index) );   
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
void print_posterior_probabilities(Observed_events *info, Forward_algorithm *fw, Backward_algorithm *bw, int start_pos, int end_pos) {
    printf("DEBUG: Entering print_posterior_probabilities function\n");
    printf("DEBUG: info->T = %d, FLANK = %d\n", info->T, FLANK);
    printf("DEBUG: start_pos = %d, end_pos = %d\n", start_pos, end_pos);
    printf("DEBUG: fw = %p, bw = %p\n", (void*)fw, (void*)bw);
    
    // Basic sanity checks first
    if (info == NULL) {
        printf("ERROR: info is NULL\n");
        return;
    }
    if (fw == NULL) {
        printf("ERROR: fw is NULL\n");
        return;
    }
    if (bw == NULL) {
        printf("ERROR: bw is NULL\n");
        return;
    }
    
    // Check matrix pointers
    printf("DEBUG: Checking fw->a = %p\n", (void*)fw->a);
    if (fw->a == NULL) {
        printf("ERROR: fw->a is NULL\n");
        return;
    }
    
    printf("DEBUG: Checking bw->b = %p\n", (void*)bw->b);
    if (bw->b == NULL) {
        printf("ERROR: bw->b is NULL\n");
        return;
    }
    
    // Check sequence length
    if (info->T <= 2 * FLANK) {
        printf("ERROR: Sequence too short for analysis (T=%d, FLANK=%d)\n", info->T, FLANK);
        return;
    }
    
    printf("\n=== Posterior Probabilities ===\n");
    printf("Position\tExon\t\tIntron\n");
    printf("----------------------------------------\n");
    
    // Calculate safe boundaries
    int max_pos = info->T - 2 * FLANK;
    printf("DEBUG: max_pos = %d\n", max_pos);
    
    // Ensure start_pos and end_pos are within bounds
    if (start_pos < 0) {
        printf("DEBUG: Adjusted start_pos from %d to 0\n", start_pos);
        start_pos = 0;
    }
    if (end_pos > max_pos) {
        printf("DEBUG: Adjusted end_pos from %d to %d\n", end_pos, max_pos);
        end_pos = max_pos;
    }
    
    printf("DEBUG: Final range: start_pos=%d, end_pos=%d\n", start_pos, end_pos);
    
    // Simpler approach with conservative checks
    for (int t = start_pos; t <= end_pos; t++) {
        printf("DEBUG: Processing position t = %d\n", t);
        
        // Check bounds thoroughly
        if (t < 0) {
            printf("DEBUG: t < 0, skipping\n");
            continue;
        }
        if (t >= info->T - 2 * FLANK + 1) {
            printf("DEBUG: t >= info->T - 2 * FLANK + 1, skipping\n");
            break;
        }
        
        // Check array access
        printf("DEBUG: Checking fw->a[%d]\n", t);
        if (fw->a[t] == NULL) {
            printf("ERROR: fw->a[%d] is NULL\n", t);
            continue;
        }
        
        printf("DEBUG: Checking bw->b[%d]\n", t);
        if (bw->b[t] == NULL) {
            printf("ERROR: bw->b[%d] is NULL\n", t);
            continue;
        }
        
        // Output position
        printf("%d\t\t", t + FLANK);
        
        // Calculate normalization - double check each step
        double pos_total = 0.0;
        for (int i = 0; i < 2; i++) {
            printf("DEBUG: Accessing fw->a[%d][%d] = %f\n", t, i, fw->a[t][i]);
            printf("DEBUG: Accessing bw->b[%d][%d] = %f\n", t, i, bw->b[t][i]);
            
            double fw_val = fw->a[t][i];
            double bw_val = bw->b[t][i];
            double product = fw_val * bw_val;
            
            printf("DEBUG: Product = %f\n", product);
            pos_total += product;
        }
        
        printf("DEBUG: pos_total = %f\n", pos_total);
        
        // Handle potential division by zero
        if (pos_total < 1e-10) {
            printf("DEBUG: pos_total too small, printing zeros\n");
            printf("0.000000\t0.000000\n");
            continue;
        }
        
        // Calculate and print posteriors
        for (int i = 0; i < 2; i++) {
            double posterior = (fw->a[t][i] * bw->b[t][i]) / pos_total;
            printf("%.6f\t", posterior);
        }
        printf("\n");
    }
    printf("----------------------------------------\n");
    printf("DEBUG: Successfully completed print_posterior_probabilities function\n");
}