#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <cblas.h>


// constructing Hidden Markov Model for isoform analysis
// given numerical definition of HMM 

// HMM = {N, M, A, B, pi}
// necessary condition we need to calcualte
// for isoform analysis, we don't need to give out the initial state of pi
// since the starting point is settle down either A T C or G

// global index explaination
// A = 1
// T = 2
// C = 3
// G = 4

typedef struct{
    int N;
    int M;
    int numerical_seq;
} Hidden_markov_model;

// note if we are not making everything into the log space or whatever turns them into the int
// the struct of lamba shall be as double

// lamba = {A, B, pi}
typedef struct{
    int A;
    int B;
    int alpha;
    int pi;
} Lamba;

// N  = number of observed state
// M  = number of hidden state
void set_parameter(Hidden_markov_model *hmm, const int N, const int M)
{
    hmm.N->N;
    hmm.M->M;
    hmm.A[4][4]; // if we are consider the base pair as observed events
    hmm.B[N][M]; // this is bug for compile, make it dynamic allocation later on
    hmm.alpha[4]; // this gonna be either 0 0 1 0 or  0 1 0 0 or 1 0 0 0 or 0 0 0 1
}

int A[4][4];

int B[4][N];

// setting up the alpha for forward-backward algorithm
// given numerical definiton of alpha
// alpha = P(sequence of event, end state of i | given condition lamba)

// first layer:  given the position of the state, which is up to T
// middle layer: given the observed state; which is up to N
// third layer:  given the hidden state;   which is up to N

int alpha[T][N][M];

// A  = transition matrix for observed state
// B  = matrix for given state, having hidden state M
// pi = 1D array for probability of given start state
// outer layer: given base pair
// inner layer: the probability of next base pair be

int main(void){

    free(hmm.numerical_seq);
}

// turns original sequence into int 
// A == 1 , C == 2, G == 3, T == 4

void numerical_transcription(Hidden_markov_model *hmm, const char *seq)
{
    
    size_t len = strlen(seq);
    hmm->numerical_seq = malloc ( len * sizeof(int));

    for( i = 0; i < len; i++)
    {
        if      (seq[i] == 'A'){
            hmm.numerical_seq[i] = 1
        }
        else if (seq[i] == 'C'){
            hmm.numerical_seq[i] = 2
        }
        else if (seq[i] == 'G'){
            hmm.numerical_seq[i] = 3
        }
        else if (seq[i] == 'T'){
            hmm.numerical_seq[i] = 4
        }
    }

}

void transition_matrix(Lamba *l)
{
    // current transition between each observed base pair is set to 0.25
    // edit directly if we wanna change soemthing 
    // or later for automatic update

    l->A = double A[4][4] = {

    //    A     C     G     T

        {0.25, 0.25, 0.25, 0.25}, // A
        {0.25, 0.25, 0.25, 0.25}, // C
        {0.25, 0.25, 0.25, 0.25}, // G
        {0.25, 0.25, 0.25, 0.25}  // T

    }

}

void starting_matrix(Lamba *l)
{
    l->pi = double pi[4] = {

    //    A     C     G     T

        0.25, 0.25, 0.25, 0.25

    }
}

void hidden_state_matrix(Lamba *l)
{
    // outer layer: given state X as base pair
    // inner layer: number of hidden state associated with
    // i assume there is 3 hidden state rn

    l->B = double B[4][3] = {
    
    //     ds    as    etc
        { 0.33, 0.33, 0.33 }, // A
        { 0.33, 0.33, 0.33 }, // C
        { 0.33, 0.33, 0.33 }, // G
        { 0.33, 0.33, 0.33 }  // T

    }
}

// parameter; given a hidden state we want to find probability on
// which is const int *hs; pointer of that

// this is calculating first 
void basis_of_forward_algorithm(Lamba *l)
{
    // alpha(current observed state) = P(observed state) * P(hidden state | observed state)
    // data structure [current position][observed state][hidden state]

    // this part require dynamic allocation
    // since we cannot tell how much position we need unless we are given a sequence 
    
    l->alpha[0][]
}

void forward_algorithm(Hidden_markov_model *HMM, Lamba *l, int *hs, int *alpha, const int *A, const *B, int depth)
{

}
