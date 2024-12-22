#include <stdio.h>
#include <stdlib.h>
#include <math.h> // Per fabs()

#define VECTOR_SIZE 256
#define EPSILON 0.001 // Soglia per confronto fino alla terza cifra decimale

typedef struct {
    int header1;
    int header2;
    double data[VECTOR_SIZE];
} BinaryFileContent;

void printFileContent(const BinaryFileContent* content) {
    printf("Header1: %d\n", content->header1);
    printf("Header2: %d\n", content->header2);
    printf("Data (double):\n");
    for (int i = 0; i < VECTOR_SIZE; i++) {
        printf(" %.3f", content->data[i]);
        if ((i + 1) % 8 == 0) printf("\n");
    }
    printf("\n");
}

void compareFiles(const BinaryFileContent* file1, const BinaryFileContent* file2, int* identicalCount, int* differentCount) {
    *identicalCount = 0;
    *differentCount = 0;

    for (int i = 0; i < VECTOR_SIZE; i++) {
        if (fabs(file1->data[i] - file2->data[i]) <= EPSILON) {
            (*identicalCount)++;
        } else {
            (*differentCount)++;
        }
    }
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Uso: %s <file1> <file2>\n", argv[0]);
        return EXIT_FAILURE;
    }

    const char* file1Name = argv[1];
    const char* file2Name = argv[2];

    FILE* file1 = fopen(file1Name, "rb");
    FILE* file2 = fopen(file2Name, "rb");
    if (!file1 || !file2) {
        perror("Errore nell'apertura dei file");
        return EXIT_FAILURE;
    }

    BinaryFileContent content1;
    BinaryFileContent content2;

    // Leggi il contenuto del primo file
    if (fread(&content1, sizeof(BinaryFileContent), 1, file1) != 1) {
        fprintf(stderr, "Errore nella lettura del primo file.\n");
        fclose(file1);
        fclose(file2);
        return EXIT_FAILURE;
    }

    // Leggi il contenuto del secondo file
    if (fread(&content2, sizeof(BinaryFileContent), 1, file2) != 1) {
        fprintf(stderr, "Errore nella lettura del secondo file.\n");
        fclose(file1);
        fclose(file2);
        return EXIT_FAILURE;
    }

    fclose(file1);
    fclose(file2);

    // Stampa i contenuti dei file
    printf("Contenuto del primo file:\n");
    printFileContent(&content1);

    printf("\nContenuto del secondo file:\n");
    printFileContent(&content2);

    // Confronta i file
    int identicalCount = 0, differentCount = 0;
    compareFiles(&content1, &content2, &identicalCount, &differentCount);

    // Risultati
    printf("\nRisultati del confronto:\n");
    printf("Numeri identici: %d\n", identicalCount);
    printf("Numeri diversi: %d\n", differentCount);

    return EXIT_SUCCESS;
}
