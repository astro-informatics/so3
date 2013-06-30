#include <assert.h>
#include <stdbool.h>
#include <stdio.h>

#include "../so3_types.h"
#include "../so3_sampling.h"

static void test_sampling_elmn2ind();

int main() {
    test_sampling_elmn2ind();
    printf("All unit tests passed!\n");
    return 0;
}

void test_sampling_elmn2ind()
{
    int ind;
    // Test padded storage with n-order 0, -1, 1, -2, 2, ...
    // We pass in a nonsensical N to make sure implementation does
    // not depend on it.
    so3_sampling_elmn2ind(&ind, 0, 0, 0, 3, -1, SO3_STORE_ZERO_FIRST_PAD);
    assert( ind == 0 && 
            "Element (0,0,0) should be first." );

    so3_sampling_elmn2ind(&ind, 2, 1, 0, 3, -1, SO3_STORE_ZERO_FIRST_PAD);
    assert( ind == 7 && 
            "Element (2,1,0) in wrong position for L = 3." );

    so3_sampling_elmn2ind(&ind, 2, 1, -1, 3, -1, SO3_STORE_ZERO_FIRST_PAD);
    assert( ind == 16 && 
            "Element (2,1,-1) in wrong position for L = 3." );

    so3_sampling_elmn2ind(&ind, 2, 1, 1, 3, -1, SO3_STORE_ZERO_FIRST_PAD);
    assert( ind == 25 && 
            "Element (2,1,1) in wrong position for L = 3." );

    // Test compact storage with n-order 0, -1, 1, -2, 2, ...
    // We pass in a nonsensical N to make sure implementation does
    // not depend on it.
    so3_sampling_elmn2ind(&ind, 0, 0, 0, 3, -1, SO3_STORE_ZERO_FIRST_COMPACT);
    assert( ind == 0 && 
            "Element (0,0,0) should be first." );

    so3_sampling_elmn2ind(&ind, 2, 1, 0, 3, -1, SO3_STORE_ZERO_FIRST_COMPACT);
    assert( ind == 7 && 
            "Element (2,1,0) in wrong position for L = 3." );

    so3_sampling_elmn2ind(&ind, 2, 1, -1, 3, -1, SO3_STORE_ZERO_FIRST_COMPACT);
    assert( ind == 15 && 
            "Element (2,1,-1) in wrong position for L = 3." );

    so3_sampling_elmn2ind(&ind, 2, 1, 1, 3, -1, SO3_STORE_ZERO_FIRST_COMPACT);
    assert( ind == 23 && 
            "Element (2,1,1) in wrong position for L = 3." );
   
    so3_sampling_elmn2ind(&ind, 2, 1, -2, 3, -1, SO3_STORE_ZERO_FIRST_COMPACT);
    assert( ind == 28 && 
            "Element (2,1,-2) in wrong position for L = 3." );

    so3_sampling_elmn2ind(&ind, 2, 1, 2, 3, -1, SO3_STORE_ZERO_FIRST_COMPACT);
    assert( ind == 33 && 
            "Element (2,1,2) in wrong position for L = 3." );

    // Test padded storage with n-order ..., -2, -1, 0, 1, 2, ...
    so3_sampling_elmn2ind(&ind, 0, 0, 0, 3, 3, SO3_STORE_NEG_FIRST_PAD);
    assert( ind == 18 && 
            "Element (0,0,0) in wrong position for L = N = 3." );

    so3_sampling_elmn2ind(&ind, 2, -2, -2, 3, 3, SO3_STORE_NEG_FIRST_PAD);
    assert( ind == 4 && 
            "Element (2,-2,-2) in wrong position for L = N = 3." );

    so3_sampling_elmn2ind(&ind, 2, 1, 0, 3, 3, SO3_STORE_NEG_FIRST_PAD);
    assert( ind == 25 && 
            "Element (2,1,0) in wrong position for L = N = 3." );

    so3_sampling_elmn2ind(&ind, 2, 1, -1, 3, 3, SO3_STORE_NEG_FIRST_PAD);
    assert( ind == 16 && 
            "Element (2,1,-1) in wrong position for L = N = 3." );

    so3_sampling_elmn2ind(&ind, 2, 1, 1, 3, 3, SO3_STORE_NEG_FIRST_PAD);
    assert( ind == 34 && 
            "Element (2,1,1) in wrong position for L = N = 3." );

    // Test compact storage with n-order ..., -2, -1, 0, 1, 2, ...
    so3_sampling_elmn2ind(&ind, 0, 0, 0, 3, 3, SO3_STORE_NEG_FIRST_COMPACT);
    assert( ind == 13 && 
            "Element (0,0,0) in wrong position for L = N = 3." );

    so3_sampling_elmn2ind(&ind, 2, -2, -2, 3, 3, SO3_STORE_NEG_FIRST_COMPACT);
    assert( ind == 0 && 
            "Element (2,-2,-2) in wrong position for L = N = 3." );

    so3_sampling_elmn2ind(&ind, 2, 1, 0, 3, 3, SO3_STORE_NEG_FIRST_COMPACT);
    assert( ind == 20 && 
            "Element (2,1,0) in wrong position for L = N = 3." );

    so3_sampling_elmn2ind(&ind, 2, 1, -1, 3, 3, SO3_STORE_NEG_FIRST_COMPACT);
    assert( ind == 11 && 
            "Element (2,1,-1) in wrong position for L = N = 3." );

    so3_sampling_elmn2ind(&ind, 2, 1, 1, 3, 3, SO3_STORE_NEG_FIRST_COMPACT);
    assert( ind == 28 && 
            "Element (2,1,1) in wrong position for L = N = 3." );
   
    so3_sampling_elmn2ind(&ind, 2, 1, -2, 3, 3, SO3_STORE_NEG_FIRST_COMPACT);
    assert( ind == 3 && 
            "Element (2,1,-2) in wrong position for L = N = 3." );

    so3_sampling_elmn2ind(&ind, 2, 1, 2, 3, 3, SO3_STORE_NEG_FIRST_COMPACT);
    assert( ind == 33 && 
            "Element (2,1,2) in wrong position for L = N = 3." );
}
