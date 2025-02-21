#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>

// constructing Hidden Markov Model for isoform analysis
// given numerical definition of HMM 

// HMM = {N, M, A, B, pi}
// necessary condition we need to calcualte

// global index explaination
// A = 1
// T = 2
// C = 3
// G = 4

typedef struct{
    int N;
    int M;
    int *pA;
    int *pB;
    int *palpha;    
} Hidden_markov_model;

// note if we are not making everything into the log space or whatever turns them into the int
// the struct of lamba shall be as double

// lamba = {A, B, pi}
typedef struct{
    int A;
    int B;
    int alpha;
} lamba;

// N  = number of observed state
// M  = number of hidden state
void set_parameter(Hidden_markov_model *hmm, const int N, const int M)
{
    hmm.N->N;
    hmm.M->M;
    hmm.A[4][4]; // if we are consider the base pair as observed events
    hmm.B[N][M]; // this is bug for compile, make it dynamic allocation later on
    hmm.alpha[3]
}

// A  = transition matrix for observed state
// B  = matrix for given state, having hidden state M
// pi = 1D array for probability of given start state
// outer layer: given base pair
// inner layer: the probability of next base pair be
void 

int A[4][4];

for i in range()

// outer layer: given state X as base pair
// inner layer: number of hidden state associated with

int B[4][N];

// setting up the alpha for forward-backward algorithm
// given numerical definiton of alpha
// alpha = P(sequence of event, end state of i | given condition lamba)

// first layer:  given the position of the state, which is up to T
// middle layer: given the observed state; which is up to N
// third layer:  given the hidden state;   which is up to N

int alpha[T][N][M];


void forward_algorithm(Hidden_markov_model *HMM, int *alpha, const int *A, const *B, int depth)
{

}