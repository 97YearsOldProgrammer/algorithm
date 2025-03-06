#ifndef HMM_MODEL
#define HMM_MODEL

#define OS 4
#define HS 13                   // 1 (exon) + 1 (intron) + 5 (donor site) + 6(acceptor site)
#define OT                      // used for print function; see which probability spot you wanna see

typedef struct                  // observed events with length T
{
    int T;                      // overall length for sequence
    int numerical_sequence;     // transcribe from base pair to digits
} Ot;

typedef struct
{
    double dons[5][4];          // the emission probability for donor sites
    double intron[6][4];        // the emission probability for intron
    double accs[256];           // the emission probability for acceptor sites
    double exon[256];           // the emission probability for exon
} emission_matrix;

// transition probability matrix illustration
// exon donor_sites1...donor_sites5 intron acceptor_site1...acceptor_sites6

// numerical abbreviation;
// exon = 0;
// ds_sites1-5 = 1-5;
// ac_sites1-6 = 6-11;
// intron = 12;

typedef struct
{
    double A;                   // the transition probability
    emission_matrix B;          // the pre-defined emission probibility data strcuture
    double pi;                  // the initial probability
} Lambda;

typedef struct
{
    double ed_exon[1000];             // the ed probability for exon
    double ed_intron[1000];           // the ed probability for intron
} explicit_duration;

// declared function
void donor_parser(Lambda *l, char *filename);
void acceptor_parser(Lambda *l, char *filename);
void exon_intron_parser(Lambda *l, char *filename);
void explicit_duration_probability(explicit_duration *ed, char *filename, int digit);

#endif