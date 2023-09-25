#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <unistd.h> 
#include <time.h> 
#include <math.h>

#define NUM_STEPS 10
#define NUM_FISH 10000000
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
double CollectiveAction(FISH* fishes, int num_fish, double total_obj_func) {
    
    double total_distance_times_weight = 0.0;

    #pragma omp parallel
    {
        #pragma omp for
        for (int i = 0; i < num_fish; i++) 
        {

            double weight = fishes[i].weight;
            double distance = sqrt(fishes[i].x * fishes[i].x + fishes[i].y * fishes[i].y);

            total_distance_times_weight += (distance * weight);

        }
    }

    // Calculate the barycenter
    double barycenter = total_distance_times_weight / total_obj_func;
}

// void swim(FISH fish)
// {
//     int random_x_int = rand() % 201;
//     int random_y_int = rand() % 201;

//     double random_x_movement = (random_x_int / 1000.0) - 0.1;
//     double random_y_movement = (random_y_int / 1000.0) - 0.1;

//     fish.prev_x = fish.x;
//     fish.prev_y = fish.y;

//     // printf("Fish previous coordinates: %d, %d\n", fish.prev_x, fish.prev_y);

//     fish.x = fish.x + random_x_movement;
//     fish.y = fish.y + random_y_movement; 

//     // printf("Fish current coordinates: %d, %d\n", fish.x, fish.y);
// }

// Function for parallelized fish swimming
void swim(FISH* fishes, int num_fish) {
    #pragma omp parallel 
    {
        #pragma omp for
        for (int i = 0; i < num_fish; i++) 
        {   
            int random_x_int = rand() % 201;
            int random_y_int = rand() % 201;
                
            double random_x_movement = (random_x_int / 1000.0) - 0.1;
            double random_y_movement = (random_y_int / 1000.0) - 0.1;
                
            fishes[i].prev_x = fishes[i].x;
            fishes[i].prev_y = fishes[i].y;
                
            fishes[i].x = fishes[i].x + random_x_movement;
            fishes[i].y = fishes[i].y + random_y_movement;
             
        }
    }
}

// void weight_function (FISH* fishes)
// {
//     // double weight_func_val;

//     // Calculate the quotient of a fish's change in objective function 
//     // and the maximum change in objective function when accounting for all fish 

//     double maxDelta = 0.0;

//     // Calculate the maximum change in the objective function
//     for (int i = 0; i < NUM_FISH; i++) {
//         fishes[i].delta_f_i = fishes[i].f_i - fishes[i].prev_f_i;
        
//         if (fishes[i].delta_f_i > maxDelta) 
//         {
//             maxDelta = fishes[i].delta_f_i;
//         }
//     }

//     // printf("Max delta: %f\n", maxDelta);

//     // Calculate the weight function value for all fish
//     double weight_func_val = 0.0;
//     for (int i = 0; i < NUM_FISH; i++) 
//     {
//         fishes[i].weight += (fishes[i].delta_f_i / maxDelta);
//     }
// }

// Function for parallelized weight_function
void weight_function(FISH* fishes, int num_fish) {
    double maxDelta = 0.0;

    #pragma omp parallel
    { 
        #pragma omp for // reduction(max:maxDelta)
        for (int i = 0; i < num_fish; i++) {
            fishes[i].delta_f_i = fishes[i].f_i - fishes[i].prev_f_i;

            if (fishes[i].delta_f_i > maxDelta) {
                maxDelta = fishes[i].delta_f_i;
            }
        }
    // printf("Max delta: %f\n", maxDelta);
    }

    #pragma omp parallel 
    {   
        #pragma omp for
        for (int i = 0; i < num_fish; i++) {
            fishes[i].weight += (fishes[i].delta_f_i / maxDelta);
        }
    }
}

double calc_euc_dist (FISH fish)
{
    return sqrt(pow((double)fish.x, 2) + pow((double)fish.y, 2));
}

double obj_func (FISH* fishes) 
{
    // Thinking we don't need this for this function; we just iterate through the objective functions first
    double total_sum = 0;
    // int pre_root_val;
    // double post_root_val;

    #pragma omp parallel
    {
        #pragma omp for 
        for (int i = 0; i < NUM_FISH; i++)
        {
            // total_sum += calc_euc_dist(fishes[i]);

            // Set the current f_i to the previous f_i
            fishes[i].prev_f_i = fishes[i].f_i;
            // printf("Fish #%d previous objective function: %f\n", i+1, fishes[i].prev_f_i);

            // Get the distance of the current fish from the center
            double distance = (double)(fishes[i].x * fishes[i].x + fishes[i].y * fishes[i].y);

            // Get the objective function of the current fish
            double obj_func_val = sqrt(distance);
            // printf("Fish #%d new objective function: %f\n", i+1, obj_func_val);
            // Store previous objective function value
            fishes[i].prev_f_i = fishes[i].f_i;

            // Set the current objective function for the fish
            fishes[i].f_i = obj_func_val;

            // Add the current f_i to the overall objective function
            total_sum += fishes[i].f_i;
        }
    }   

    return total_sum;
}

// void main(int argc, char* argv[])
int main(int argc, char* argv[])
{
    FISH *fishes;

    fishes = (FISH*) malloc(NUM_FISH * sizeof(FISH));

    double start = omp_get_wtime();
    
    // Generate positions for the fish

    #pragma omp parallel
    {
        #pragma omp for
        for (int i = 0; i < NUM_FISH; i++)
        {
            int x_rand_num = rand() % 201 - 100;
            int y_rand_num = rand() % 201 - 100;

            fishes[i].x = x_rand_num;
            fishes[i].y = y_rand_num;
            fishes[i].weight = FISH_INIT_WEIGHT; 
            fishes[i].f_i = calc_euc_dist(fishes[i]);

            // printf("Fish #%d coordinates: (%d, %d)\n", i+1, fishes[i].x, fishes[i].y);
        }
    }

    double total_sum;
    
    // for collective action
    double barycenter_x;
    double barycenter_y;
    // ------

    for (int i = 0; i < NUM_STEPS; i++)
    {
        
        if (i == 0)
        {
            // printf("Calculating the weight function...\n");
            
            // Weight function is random value at the very first step
            total_sum = obj_func(fishes);

            #pragma omp parallel 
            {
                #pragma omp for 
                for (int j = 0; j < NUM_FISH; j++)
                {
                    double current_weight = fishes[j].weight;
                    double weight_func = rand() % 5 - 1;
                    fishes[j].weight += weight_func;

                }
            }
            
            swim(fishes, NUM_FISH);
        }
        else
        {
            // ADD ALL WEIGHTS FOR FISH 
            total_sum = obj_func(fishes);
            // weight_function(fishes); 
            weight_function(fishes, NUM_FISH);  

            swim(fishes, NUM_FISH);

            double barycentre = CollectiveAction(fishes, NUM_FISH, total_sum);
            
            // printf("Barycentre at Step %d: %f\n", i+1, barycentre);
        }
        
        // printf("Objective function at Step %d: %f\n", i+1, total_sum);
        
    }

    double end = omp_get_wtime();

    double time_taken= end-start;

    free(fishes);

    printf("Time spent: %.2f\n", time_taken); 
}