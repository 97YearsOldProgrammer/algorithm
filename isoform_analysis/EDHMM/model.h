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

typedef struct emission_matrix
{
    double dons;                // the emission probability for donor sites
    double accs;                // the emission probability for acceptor sites
    double exon;                // the emission probability for exon
    double intron;              // the emission probability for intron
};

// transition probability matrix illustration
// exon donor_sites1...donor_sites5 intron acceptor_site1...acceptor_sites6

// numerical abbreviation;
// exon = 0;
// ds_sites1-5 = 1-5;
// ac_sites1-6 = 6-11;
// intron = 12;

typedef struct Lambda
{
    double A;                   // the transition probability
    emission_matrix B;          // the pre-defined emission probibility data strcuture
    double pi;                  // the initial probability
};

typedef struct explicit_duration
{
    double ed_exon;             // the ed probability for exon
    double ed_intron;           // the ed probability for intron
}

#endif