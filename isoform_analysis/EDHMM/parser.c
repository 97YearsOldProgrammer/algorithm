#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Structure to represent a Position Weight Matrix (PWM)
typedef struct {
    char *name;        // Name of the model
    int rows;          // Number of positions (rows)
    int cols;          // Number of nucleotides (columns) - typically 4 for A,C,G,T
    double **matrix;   // 2D array to store the probability values
} PWM;

// Function to allocate memory for a PWM
PWM* create_pwm(const char *name, int rows, int cols) {
    PWM *pwm = (PWM*)malloc(sizeof(PWM));
    if (!pwm) {
        fprintf(stderr, "Memory allocation failed for PWM struct\n");
        return NULL;
    }
    
    pwm->name = strdup(name);
    pwm->rows = rows;
    pwm->cols = cols;
    
    // Allocate memory for the 2D matrix
    pwm->matrix = (double**)malloc(rows * sizeof(double*));
    if (!pwm->matrix) {
        fprintf(stderr, "Memory allocation failed for matrix rows\n");
        free(pwm->name);
        free(pwm);
        return NULL;
    }
    
    for (int i = 0; i < rows; i++) {
        pwm->matrix[i] = (double*)malloc(cols * sizeof(double));
        if (!pwm->matrix[i]) {
            fprintf(stderr, "Memory allocation failed for matrix column %d\n", i);
            // Free previously allocated memory
            for (int j = 0; j < i; j++) {
                free(pwm->matrix[j]);
            }
            free(pwm->matrix);
            free(pwm->name);
            free(pwm);
            return NULL;
        }
    }
    
    return pwm;
}

// Function to free memory allocated for a PWM
void free_pwm(PWM *pwm) {
    if (pwm) {
        if (pwm->matrix) {
            for (int i = 0; i < pwm->rows; i++) {
                free(pwm->matrix[i]);
            }
            free(pwm->matrix);
        }
        free(pwm->name);
        free(pwm);
    }
}

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

int main(int argc, char argv[]) {
    // Load the acceptor PWM
    PWM *acc_pwm = load_pwm_from_file("isoform_analysis/models/acc.pwm");

    if (!acc_pwm) {
        fprintf(stderr, "Failed to load acceptor PWM\n");
        return 1;
    }
    
    // Load the donor PWM
    PWM *don_pwm = load_pwm_from_file("isoform_analysis/models/don.pwm");
    if (!don_pwm) {
        fprintf(stderr, "Failed to load donor PWM\n");
        free_pwm(acc_pwm);
        return 1;
    }
    
    // Print the loaded PWMs for verification
    printf("Successfully loaded PWM models:\n\n");
    print_pwm(acc_pwm);
    print_pwm(don_pwm);
    
    // Clean up
    free_pwm(acc_pwm);
    free_pwm(don_pwm);
    
    return 0;
}