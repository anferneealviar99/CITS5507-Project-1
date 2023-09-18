#include <stdio.h>
#include <stdlib.h>
#include <unistd.h> 
#include <time.h> 
#include <math.h>

#define NUM_STEPS 10
#define NUM_FISH 5
#define FISH_INIT_WEIGHT 15

// Declare structure for fish, holding coordinates (for now)
typedef struct _fish {
    int prev_x;
    int prev_y;
    int x;
    int y;
    double prev_f_i;
    double f_i;
    double delta_f_i;
    double prev_weight;
    double weight;  
} FISH;

//TODO Collective Action

void swim(FISH fish)
{
    int random_x_int = rand() % 201;
    int random_y_int = rand() % 201;

    double random_x_movement = (random_x_int / 1000.0) - 0.1;
    double random_y_movement = (random_y_int / 1000.0) - 0.1;

    fish.prev_x = fish.x;
    fish.prev_y = fish.y;

    fish.x = fish.x + random_x_movement;
    fish.y = fish.y + random_y_movement; 
}

double weight_function_calc (FISH* fishes)
{
    // double weight_func_val;

    // Calculate the quotient of a fish's change in objective function 
    // and the maximum change in objective function when accounting for all fish 

    double maxDelta = 0.0;

    // Calculate the maximum change in the objective function
    for (int i = 0; i < NUM_FISH; i++) {
        fishes[i].delta_f_i = fishes[i].f_i - fishes[i].prev_f_i;
        
        if (fishes[i].delta_f_i > maxDelta) 
        {
            maxDelta = fishes[i].delta_f_i;
        }
    }

    // Calculate the weight function value for all fish
    double weight_func_val = 0.0;
    for (int i = 0; i < NUM_FISH; i++) 
    {
        fishes[i].weight += fishes[i].delta_f_i / maxDelta;
        weight_func_val += fishes[i].weight;
    }

    return weight_func_val;
}

double calc_euc_dist (FISH fish)
{
    return sqrt(pow((double)fish.x, 2) + pow((double)fish.y, 2));
}

double obj_func (FISH* fishes) 
{
    double total_sum = 0;
    // int pre_root_val;
    // double post_root_val;

    for (int i = 0; i < NUM_FISH; i++)
    {
        // total_sum += calc_euc_dist(fishes[i]);

        // Set the current f_i to the previous f_i
        fishes[i].prev_f_i = fishes[i].f_i;

        // Get the distance of the current fish from the center
        double distance = (double)(fishes[i].x * fishes[i].x + fishes[i].y * fishes[i].y);

        // Get the objective function of the current fish
        double obj_func_val = sqrt(distance);

        // Set the current objective function for the fish
        fishes[i].f_i = obj_func_val;

        // Add the current f_i to the overall objective function
        total_sum += fishes[i].f_i;
    }
    
    return total_sum;
}

// void main(int argc, char* argv[])
int main(int argc, char* argv[])
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
        fishes[i].weight = FISH_INIT_WEIGHT; 
        fishes[i].f_i = calc_euc_dist(fishes[i]);

        printf("Fish #%d coordinates: (%d, %d)\n", i+1, fishes[i].x, fishes[i].y);
    }

    double weight_func = 0;
    for (int i = 0; i < NUM_STEPS; i++)
    {
        double total_sum = obj_func(fishes);

        if (i == 0)
        {
            // Weight function is random value at the very first step
            
            weight_func = rand() % 5 - 1;
            
            for (int j = 0; j < NUM_FISH; j++)
            {
                swim(fishes[j]);
            }
            
        }
        else
        {
            // Calculate the objective functions to get the weight function 
            weight_func = weight_function_calc(fishes);   
        }
    
        printf("Objective function at Step %d: %.2f\n", i+1, total_sum);
    }

    clock_t end = clock();

    double time_spent = (double) (end-begin) / CLOCKS_PER_SEC;

    free(fishes);

    printf("Time spent: %.2f\n", time_spent); 
}