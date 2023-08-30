#include <stdio.h>
#include <stdlib.h>
#include <unistd.h> 
#include <time.h> 
#include <math.h>

#define NUM_STEPS 10
#define NUM_FISH 5

// Declare structure for fish, holding coordinates (for now)
typedef struct _fish {
    int x;
    int y;
    double weight; 
} FISH;

double obj_func (FISH* fishes) 
{
    double total_sum;
    int pre_root_val;
    double post_root_val;

    for (int i = 0; i < NUM_FISH; i++)
    {
        pre_root_val = pow((double)fishes[i].x, 2) + pow((double)fishes[i].y, 2);
        post_root_val = sqrt((double) pre_root_val);

        total_sum += post_root_val;
    }

    return total_sum;
}

void main(int argc, char* argv[])
{
    FISH *fishes;

    fishes = (FISH*) malloc(NUM_FISH * sizeof(FISH));

    // fishes->x=1;
    // fishes->y=1;

    // (fishes+1)->x=2;
    // (fishes+1)->y=2;

    // printf("x coord of 1st fish=%d\n", fishes->x);
    // printf("x coord of 2nd fish=%d\n", (fishes+1)->x); 


    // free(fishes);

    srand(time(NULL));

    clock_t begin = clock();
    // Generate positions for the fish
    for (int i = 0; i < NUM_FISH; i++)
    {
        int x_rand_num = rand() % 201 - 100;
        int y_rand_num = rand() % 201 - 100;

        fishes[i].x = x_rand_num;
        fishes[i].y = y_rand_num;

        printf("Fish #%d coordinates: (%d, %d)\n", i+1, fishes[i].x, fishes[i].y);
    }

    for (int i = 0; i < NUM_STEPS; i++)
    {
        double total_sum = obj_func(fishes);

        printf("Objective function at Step %d: %.2f\n", i+1, total_sum);
    }

    clock_t end = clock();

    double time_spent = (double) (end-begin) / CLOCKS_PER_SEC;

    printf("Time spent: %.2f\n", time_spent); 
}