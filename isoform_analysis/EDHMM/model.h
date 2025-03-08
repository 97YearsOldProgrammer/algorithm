#ifndef HMM_MODEL
#define HMM_MODEL

#define OS 4
#define HS 13                           // 1 (exon) + 1 (intron) + 5 (donor site) + 6(acceptor site)
#define OT                              // used for print function; see which probability spot you wanna see

typedef struct                          // observed events with length T
{
    int T;                              // overall length for sequence
    int numerical_sequence;             // transcribe from base pair to digits
} Observed_events;

typedef struct
{
    double dons[5][4];                  // the emission probability for donor sites
    double intron[6][4];                // the emission probability for intron
    double accs[256];                   // the emission probability for acceptor sites
    double exon[256];                   // the emission probability for exon
} Emission_matrix;

typedef struct                          // degrade the sequence of conventional transition prob from donor 1-5 acceptor 1-6
{
    double dons[1024];                  // enumerating exon->intron ; aka donor site series
    double accs[4096];                  // enumerating intron->exon ; aka acceptor site series                       
} Transition_matrix;

typedef struct 
{
    double prob[6];                     // for apc algorithm to calculate transition prob
    int position[6];                    // for apc algorithm to get index to store in transition matrix
} Apc;


typedef struct
{
    Transition_matrix A;                // the transition probability
    Emission_matrix B;                  // the pre-defined emission probibility data strcuture
    double pi;                          // the initial probability
} Lambda;

typedef struct
{
    double exon[1000];                  // the ed probability for exon
    double intron[1000];                // the ed probability for intron
} Explicit_duration;

typedef struct
{
    double  a;                          // alpha component for forward algorithm
    double log_fw[HS];                  // prepared for log softmax trick   
} Forward_algorithm;

// declared function
void donor_parser(Lambda *l, char *filename);
void acceptor_parser(Lambda *l, char *filename);
void exon_intron_parser(Lambda *l, char *filename);
void explicit_duration_probability(explicit_duration *ed, char *filename, int digit);

#endif