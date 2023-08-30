#include <stdio.h>
#include <stdlib.h>
#include <unistd.h> 

#define NUM_STEPS 10
#define NUM_FISH 100

// Declare structure for fish, holding coordinates (for now)
typedef struct _fish {
    int x;
    int y;
} FISH;

void main(int argc, char* argv[])
{
    FISH *fishes;

    fishes = (FISH*) malloc(NUM_FISH * sizeof(FISH));

    fishes->x=1;
    fishes->y=1;

    (fishes+1)->x=2;
    (fishes+1)->y=2;

    printf("x coord of 1st fish=%d\n", fishes->x);
    printf("x coord of 2nd fish=%d\n", (fishes+1)->x); 


    free(fishes);
}