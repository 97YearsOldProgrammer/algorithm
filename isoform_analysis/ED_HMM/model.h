#ifndef HMM_MODEL
#define HMM_MODEL

#define OS 4
#define HS 13                   // 1 (exon) + 1 (intron) + 5 (donor site) + 6(acceptor site)
#define OT                      // used for print function; see which probability spot you wanna see

typedef struct Ot               // observed events with length T
{
    int T;                      // overall length for sequence
    int numerical_sequence;     // transcribe from base pair to digits
};

typedef struct Lambda
{
    double A[][];               // the transition probability
    double pi[];                // the initial probability
};

typedef struct B
{
    double dons[][];          // the emission probability for donor sites
    double accs[][];          // the emission probability for acceptor sites
    double exon[][];          // the emission probability for exon
    double intron[][];        // the emission probability for intron
};

#endif