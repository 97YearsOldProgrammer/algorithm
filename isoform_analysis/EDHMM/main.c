#include "model.h"

int main(int argc, char *argv[])
{
    if (argc < 7)
    {
        printf("Check the requirement for EDHMM. Missing files.\n");
        printf("Proper Usage:\n");
        printf("acc.pwm don.pwm exon_emission intron_emission");
    }

    // argv section for command-line inputs
    char *don_emission;
    char *acc_emission;
    char *exon_emission;
    char *intron_emission;
    char *Ped_exon;
    char *Ped_intron;

    don_emission    = argv[1];
    acc_emission    = argv[2];
    exon_emission   = argv[3];
    intron_emission = argv[4];
    Ped_exon        = argv[5];
    Ped_intron      = argv[6];

    Lamdba l;
    explicit_duration ed;

    // initialize datas
    void donor_parser(&l, don_emission);                            // donor emission prob
    void acceptor_parser(&l, acc_emission);                         // acceptor emission prob
    void exon_intron_parser(&l, exon_emission);                     // exon emission prob
    void exon_intron_parser(&l, intron_emission);                   // intron emission prob
    void explicit_duration_probability(&ed, Ped_exon, int 0);       // exon ed prob
    void explicit_duration_probability(&ed, Ped_intron, int 1);     // intron ed prob




}