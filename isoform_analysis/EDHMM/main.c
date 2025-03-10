#include "model.h"

int main(int argc, char *argv[])
{
    if (argc < 7)
    {
        printf("Check the requirement for EDHMM. Missing files.\n");
        printf("Proper Usage:\n");
        printf("acc.pwm don.pwm exon_emission intron_emission");
    }

    // argv section for command-line inputs //
    char *don_emission;
    char *acc_emission;
    char *exon_emission;
    char *intron_emission;
    char *Ped_exon;
    char *Ped_intron;

    // handle argv value //
    don_emission    = argv[1];
    acc_emission    = argv[2];
    exon_emission   = argv[3];
    intron_emission = argv[4];
    Ped_exon        = argv[5];
    Ped_intron      = argv[6];

    // data structure //
    Observed_events info;
    Apc apc;
    Lambda l;
    Explicit_duration ed;
    Forward_algorithm fw;
    Backward_algorithm bw;

    // get sequence //
    read_sequence_file(const char *filename, Observed_events *info);
    numerical_transcription(Observed_events *info, const char *seq);
    
    // initialize datas //
    donor_parser(&l, don_emission);                            // donor emission prob
    acceptor_parser(&l, acc_emission);                         // acceptor emission prob
    exon_intron_parser(&l, exon_emission, 0);                  // exon emission prob
    exon_intron_parser(&l, intron_emission, 1);                // intron emission prob
    explicit_duration_probability(&ed, Ped_exon,   0);         // exon ed prob
    explicit_duration_probability(&ed, Ped_intron, 1);         // intron ed prob

    // initialize computation //
    setup_initial_probability(&l);                             // setup pi 
    initialize_donor_transition_matrix(&l, &apc, 0);           // setup transition prob for exon to intron
    initialize_acceptor_transition_matrix(&l, &apc, 0);        // setup transition prob for intron to exon

    // initialize algorihtm //
    initial_forward_algorithm(&l, &ed, &fw, &info);            // set up alpha 0
    initial_backward_algorithm(&l, &bw, &info);                // set up beta t

    // initialize memory //
    allocate_alpha(&info, &fw);                                 // allocate forward  algorithm
    allocate_beta(&info, &bw);                                  // allocate backward algorithm

    // forward and backward algo
    forward_algorithm(&l, &fw, &info, &ed);                     
    backward_algorithm(&l, &bw, &info, &ed);

    // free memory
    free_alpha(&info, &fw);
    free_beta(&info, &bw);

}