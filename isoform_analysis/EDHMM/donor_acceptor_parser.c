#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "model.h"


// Function to load a PWM from a file
PWM* load_pwm_from_file(const char *filename) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error opening file: %s\n", filename);
        return NULL;
    }
    
    char line[256];
    
    // Read the header line
    if (!fgets(line, sizeof(line), file)) {
        fprintf(stderr, "Error reading header from %s\n", filename);
        fclose(file);
        return NULL;
    }
    
    // Parse the header to get the number of rows
    char model_path[256];
    int rows;
    if (sscanf(line, "%% PWM %s %d", model_path, &rows) != 2) {
        fprintf(stderr, "Invalid header format in %s\n", filename);
        fclose(file);
        return NULL;
    }
    
    // Extract the model name from the path (after the last slash or the whole string)
    char *model_name = strrchr(model_path, '/');
    if (model_name) {
        model_name++; // Skip the slash
    } else {
        model_name = model_path;
    }
    
    // Create a PWM with 4 columns (A, C, G, T)
    PWM *pwm = create_pwm(model_name, rows, 4);
    if (!pwm) {
        fclose(file);
        return NULL;
    }
    
    // Read the matrix values
    for (int i = 0; i < rows; i++) {
        if (!fgets(line, sizeof(line), file)) {
            fprintf(stderr, "Error reading line %d from %s\n", i+1, filename);
            free_pwm(pwm);
            fclose(file);
            return NULL;
        }
        
        // Parse the four values on this line
        if (sscanf(line, "%lf %lf %lf %lf", 
                  &pwm->matrix[i][0], &pwm->matrix[i][1], 
                  &pwm->matrix[i][2], &pwm->matrix[i][3]) != 4) {
            fprintf(stderr, "Invalid data format at line %d in %s\n", i+1, filename);
            free_pwm(pwm);
            fclose(file);
            return NULL;
        }
    }
    
    fclose(file);
    return pwm;
}

// Function to print a PWM for verification
void print_pwm(const PWM *pwm) {
    if (!pwm) return;
    
    printf("PWM: %s (%d positions x %d nucleotides)\n", pwm->name, pwm->rows, pwm->cols);
    printf("      A      C      G      T\n");
    printf("-----------------------------------\n");
    
    for (int i = 0; i < pwm->rows; i++) {
        printf("%d | ", i+1);
        for (int j = 0; j < pwm->cols; j++) {
            printf("%6.4f ", pwm->matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

// Function to initialize emission probabilities from PWMs
void setup_emission_probability(Lambda *lambda, PWM *donor_pwm, PWM *acceptor_pwm) {
    // Allocate memory for the emission matrix
    lambda->B = (struct emission_matrix*)malloc(sizeof(struct emission_matrix));
    if (!lambda->B) {
        fprintf(stderr, "Memory allocation failed for emission matrix\n");
        return;
    }
    
    // Initialize donor site emission probabilities
    // For donor sites, we need to set up emission probabilities for positions d1-d5
    // Each position corresponds to a row in the donor PWM
    lambda->B->dons = donor_pwm;  // Store the entire PWM for access during emission calculations
    
    // Initialize acceptor site emission probabilities
    // For acceptor sites, we need to set up emission probabilities for positions a1-a6
    // Each position corresponds to a row in the acceptor PWM
    lambda->B->accs = acceptor_pwm;  // Store the entire PWM for access during emission calculations
    
    // Set default emission probabilities for exon and intron
    // These are typically simpler and might be based on nucleotide frequencies
    // For simplicity, we'll use uniform distributions (0.25 for each nucleotide)
    lambda->B->exon = 0.25;    // Placeholder - in reality, this would be a more complex model
    lambda->B->intron = 0.25;  // Placeholder - in reality, this would be a more complex model
    
    printf("Emission probabilities initialized from PWMs\n");
}

// Function to setup explicit duration probabilities
void setup_explicit_duration(Lambda *lambda) {
    // Allocate memory for explicit duration structure
    lambda->ed = (struct explicit_duration*)malloc(sizeof(struct explicit_duration));
    if (!lambda->ed) {
        fprintf(stderr, "Memory allocation failed for explicit duration\n");
        return;
    }
    
    // Set default values for explicit duration probabilities
    // These would typically be derived from biological data
    // For simplicity, we'll use placeholder values
    lambda->ed->ed_exon = 0.8;    // Probability of staying in exon state
    lambda->ed->ed_intron = 0.9;  // Probability of staying in intron state
    
    printf("Explicit duration probabilities initialized\n");
}

// Function to initialize the complete EDHMM model
Lambda* initialize_edhmm_model(const char *donor_file, const char *acceptor_file) {
    // Load PWMs
    PWM *donor_pwm = load_pwm_from_file(donor_file);
    if (!donor_pwm) {
        fprintf(stderr, "Failed to load donor PWM\n");
        return NULL;
    }
    
    PWM *acceptor_pwm = load_pwm_from_file(acceptor_file);
    if (!acceptor_pwm) {
        fprintf(stderr, "Failed to load acceptor PWM\n");
        free_pwm(donor_pwm);
        return NULL;
    }
    
    // Create and initialize Lambda structure
    Lambda *lambda = (Lambda*)malloc(sizeof(Lambda));
    if (!lambda) {
        fprintf(stderr, "Memory allocation failed for Lambda\n");
        free_pwm(donor_pwm);
        free_pwm(acceptor_pwm);
        return NULL;
    }
    
    // Initialize model components
    setup_initial_probability(lambda);
    setup_transition_probability(lambda);
    setup_emission_probability(lambda, donor_pwm, acceptor_pwm);
    setup_explicit_duration(lambda);
    
    printf("EDHMM model initialization complete\n");
    
    return lambda;
}

// Main function for testing
int main(int argc, char *argv[]) {
    const char *donor_file = "isoform_analysis/models/don.pwm";
    const char *acceptor_file = "isoform_analysis/models/acc.pwm";
    
    // Initialize the EDHMM model
    Lambda *model = initialize_edhmm_model(donor_file, acceptor_file);
    if (!model) {
        fprintf(stderr, "Failed to initialize EDHMM model\n");
        return 1;
    }
    
    // Print model information for verification
    printf("Successfully initialized EDHMM model with emission probabilities from PWMs\n");
    
    // Clean up
    // Note: In a real application, you would pass this model to other functions
    // rather than freeing it here
    free(model->pi);
    for (int i = 0; i < HS; i++) {
        free(model->A[i]);
    }
    free(model->A);
    free(model->B);
    free(model->ed);
    free(model);
    
    return 0;
}