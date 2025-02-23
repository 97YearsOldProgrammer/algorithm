// constructing Hidden Markov Model for isoform analysis
// given numerical definition of HMM 

// HMM = {N, M, A, B, pi}
// necessary condition we need to calcualte
// for isoform analysis, we don't need to give out the initial state of pi
// since the starting point is settle down either A T C or G

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>

// defining parameter for markov model
#define OS 4     // number of observed state ； OS == N
#define HS 3     // number of hidden state   ； HS == M
#define Ot 4     // given probability distribution at T = Ot

typedef struct
{
    int T;                  // T  is length for observed events
    int *numerical_seq;     // used for easier calculation
}   Hidden_markov_model;

typedef struct
{                         // lambda = {A, B, pi}
    double A[HS][HS];     // A  is the transition matrix for observed state
    double B[HS][OS];     // B  is the emission matrix for hidden state
    double pi[HS];        // pi is the initial probability for hidden state
}   Lambda;

typedef struct
{
    double **alpha;       // alpha for forward algorithm
}   Fw_algo;

typedef struct
{
    double **beta;        // beta for backward algorithm
}   Bw_algo;

/*******************\
prototype being used
\*******************/

void numerical_transcription(Hidden_markov_model *hmm, const char *seq);

// prevent overflow
double log_sum_exp(double *logs, int n);

// forward algorithm
void allocate_alpha(Hidden_markov_model *hmm, Fw_algo *fw);
void basis_of_forward_algorithm(Hidden_markov_model *hmm, Lambda *l, Fw_algo *fw);
void forward_algorithm(Hidden_markov_model *hmm, Lambda *l, Fw_algo *fw);

// backward algorithm
void allocate_beta(Hidden_markov_model *hmm, Bw_algo *bw);
void basis_of_backward_algorithm(Hidden_markov_model *hmm, Lambda *l, Bw_algo *bw);
void backward_algorithm(Hidden_markov_model *hmm, Lambda *l, Bw_algo *bw);

void free_alpha(Hidden_markov_model *hmm, Fw_algo *fw);
void free_beta(Hidden_markov_model *hmm, Bw_algo *bw);

// ====================  test sequence ==================== 

char seq[] = "ACGTTTTGCGT";

// ====================  execution     ====================

int main(){

    Hidden_markov_model hmm;

    Lambda l = {
        .A =
        {   // given transition matrix
            // ds   ac    etc     
            {0.70, 0.15, 0.15}, // ds
            {0.15, 0.70, 0.15}, // ac
            {0.05, 0.05, 0.90}, // etc
        },

        .B =
        {   // Emission matrix
            // A     C      G     T
            { 0.10, 0.10, 0.40, 0.40}, // ds
            { 0.40, 0.10, 0.40, 0.10}, // ac
            { 0.30, 0.20, 0.20, 0.30}, // etc
        },
            // ds     ac   etc
        .pi = {0.01, 0.01, 0.98}       // initial probability
    };

    Fw_algo fw;
    Bw_algo bw;

    numerical_transcription(&hmm, seq);
    allocate_alpha(&hmm, &fw);
    allocate_beta(&hmm, &bw);

    // forward algorithm
    basis_of_forward_algorithm(&hmm, &l, &fw);
    basis_of_backward_algorithm(&hmm, &l, &bw);

    forward_algorithm(&hmm, &l, &fw);
    backward_algorithm(&hmm, &l, &bw);

    double log_total = log_sum_exp(fw.alpha[hmm.T-1], HS);

    printf(" Total sequence probability (linear): %.3e\n", exp(log_total) );

    printf(" \nState probabilities at t = %d:\n", Ot );
    
    for(int i = 0; i < HS; i++) 
    {
        double log_gamma = fw.alpha[Ot][i] + bw.beta[Ot][i] - log_total;
        printf("State %d: %.4f\n", i, exp(log_gamma));
    }

    // free memory
    free_alpha(&hmm, &fw);
    free_beta(&hmm, &bw);
    free(hmm.numerical_seq);

    return 0;
}

void numerical_transcription(Hidden_markov_model *hmm, const char *seq)
{

    // turns original sequence into int 
    size_t len = strlen(seq);

    hmm->T = len;
    hmm->numerical_seq = malloc ( len * sizeof(int) );

    for( int i = 0; i < len; i++)
    {

        // A == 0 , C == 1, G == 2, T == 3  
        if      (seq[i] == 'A')     hmm->numerical_seq[i] = 0;
        else if (seq[i] == 'C')     hmm->numerical_seq[i] = 1;
        else if (seq[i] == 'G')     hmm->numerical_seq[i] = 2;
        else if (seq[i] == 'T')     hmm->numerical_seq[i] = 3;
    }

}

void allocate_alpha(Hidden_markov_model *hmm, Fw_algo *fw)
{

    assert(Ot < hmm->T && "Wrong Ot input: Ot exceeds sequence length");

    fw->alpha = malloc ( hmm->T * sizeof(double*) );

    for (int i = 0 ; i < hmm->T ; i++ )
    {
        fw->alpha[i] = calloc(HS, sizeof(double) );
    }
}

void basis_of_forward_algorithm(Hidden_markov_model *hmm, Lambda *l, Fw_algo *fw)
{
    int obs = hmm->numerical_seq[0];

    for (int i = 0; i < HS; i++)
    {
        fw->alpha[0][i] = log(l->pi[i]) + log(l->B[i][obs]);;
    }
}

double log_sum_exp(double *logs, int n) 
{
    double max_log = logs[0];

    for (int i = 1; i < n ; i++)
    {
        if( logs[i] > max_log )
        {
            max_log = logs[i];
        }
    }

    double sum = 0.0;

    for( int i = 0; i < n; i++)
        sum += exp(logs[i] - max_log);

    return max_log + log(sum);       
}

void forward_algorithm(Hidden_markov_model *hmm, Lambda *l, Fw_algo *fw)
{
    for (int t = 1; t < hmm->T; t++)
    {
        int obs = hmm->numerical_seq[t];    // getting observed value

        for (int i = 0; i < HS; i ++)
        {
            double log_fw[HS];

            for (int j = 0; j < HS; j++)
            {
                log_fw[j] = fw->alpha[t - 1][j] + log(l->A[j][i]);
            }

            double p = log_sum_exp(log_fw, HS);

            fw->alpha[t][i] = p + log( l->B[i][obs] );
        }
    }
}
void free_alpha(Hidden_markov_model *hmm, Fw_algo *fw)
{
    for (int i = 0; i < hmm->T; i ++)
    {
        free( fw->alpha[i] );
    }
    free(fw->alpha);
}

void allocate_beta(Hidden_markov_model *hmm, Bw_algo *bw)
{

    bw->beta = malloc ( hmm->T * sizeof(double*) );

    for (int i = 0 ; i < hmm->T ; i++ )
    {
        bw->beta[i] = calloc(HS, sizeof(double) );
    }
}

void basis_of_backward_algorithm(Hidden_markov_model *hmm, Lambda *l, Bw_algo *bw)
{
    int last_idx = hmm->T - 1;

    for (int i = 0; i < HS; i++)
    {
        bw->beta[last_idx][i] = 0.0;  // set them to 1.0; but log(1) = 0.0
    }
}

void backward_algorithm(Hidden_markov_model *hmm, Lambda *l, Bw_algo *bw)
{
    for (int t = (hmm->T - 2); t >= 0; t--)
    {
        int obs = hmm->numerical_seq[t + 1];    // getting observed value

        for (int i = 0; i < HS; i ++)
        {
            double log_bw[HS];

            for (int j = 0; j < HS; j++)
            {
                log_bw[j] = log(l->A[i][j]) + log(l->B[j][obs]) + bw->beta[t+1][j];
            }

            bw->beta[t][i] = log_sum_exp(log_bw, HS);

        }
    }
}

void free_beta(Hidden_markov_model *hmm, Bw_algo *bw)
{
    for (int i = 0; i < hmm->T; i ++)
    {
        free( bw->beta[i] );
    }
    free(bw->beta); 
}