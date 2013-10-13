#include <assert.h>
#include <stdbool.h>
#include <stdio.h>

#include "../so3_types.h"
#include "../so3_sampling.h"

static void test_sampling_elmn2ind();
static void test_sampling_ind2elmn();
static void test_sampling_elmn2ind_real();
static void test_sampling_ind2elmn_real();

int main() {
    test_sampling_elmn2ind();
    test_sampling_ind2elmn();
    test_sampling_elmn2ind_real();
    test_sampling_ind2elmn_real();
    printf("All unit tests passed!\n");
    return 0;
}

void test_sampling_elmn2ind()
{
    int ind;
    // Test padded storage with n-order 0, -1, 1, -2, 2, ...
    // That is, for L = N = 2, the flmn are layed out as follows.
    // Each cell in the m? row corresponds to one coefficient in memory.
    // It contains an underscore for invalid el-m-n combinations
    // (which will be zeroes in memory) or the actual m-index for values
    // that will really be stored.
    //
    // n  |       0       |      -1       |       1       |
    // el | 0 |     1     | 0 |     1     | 0 |     1     |
    // m? | 0 |-1 | 0 | 1 | _ |-1 | 0 | 1 | _ |-1 | 0 | 1 |
    //
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
    // That is, for L = N = 2, the flmn are layed out as follows.
    // Each cell of the m-row corresponds to one coefficient in memory.
    //
    // n  |       0       |    -1     |     1     |
    // el | 0 |     1     |     1     |     1     |
    // m  | 0 |-1 | 0 | 1 |-1 | 0 | 1 |-1 | 0 | 1 |
    //
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
    // That is, for L = N = 2, the flmn are layed out as follows.
    // Each cell in the m? row corresponds to one coefficient in memory.
    // It contains an underscore for invalid el-m-n combinations
    // (which will be zeroes in memory) or the actual m-index for values
    // that will really be stored.
    //
    // n  |      -1       |       0       |       1       |
    // el | 0 |     1     | 0 |     1     | 0 |     1     |
    // m? | _ |-1 | 0 | 1 | 0 |-1 | 0 | 1 | _ |-1 | 0 | 1 |
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
    // That is, for L = N = 2, the flmn are layed out as follows.
    // Each cell of the m-row corresponds to one coefficient in memory.
    //
    // n  |    -1     |       0       |     1     |
    // el |     1     | 0 |     1     |     1     |
    // m  |-1 | 0 | 1 | 0 |-1 | 0 | 1 |-1 | 0 | 1 |
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

void test_sampling_ind2elmn()
{
    int el, m, n;
    // Test padded storage with n-order 0, -1, 1, -2, 2, ...
    // That is, for L = N = 2, the flmn are layed out as follows.
    // Each cell in the m? row corresponds to one coefficient in memory.
    // It contains an underscore for invalid el-m-n combinations
    // (which will be zeroes in memory) or the actual m-index for values
    // that will really be stored.
    //
    // n  |       0       |      -1       |       1       |
    // el | 0 |     1     | 0 |     1     | 0 |     1     |
    // m? | 0 |-1 | 0 | 1 | _ |-1 | 0 | 1 | _ |-1 | 0 | 1 |
    //
    // We pass in a nonsensical N to make sure implementation does
    // not depend on it.
    so3_sampling_ind2elmn(&el, &m, &n, 0, 3, -1, SO3_STORE_ZERO_FIRST_PAD);
    assert( el == 0 && m == 0 && n == 0 &&
            "Index 0 yields wrong indices (el,m,n)." );

    so3_sampling_ind2elmn(&el, &m, &n, 7, 3, -1, SO3_STORE_ZERO_FIRST_PAD);
    assert( el == 2 && m == 1 && n == 0 &&
            "Index 7 yields wrong indices (el,m,n)." );

    so3_sampling_ind2elmn(&el, &m, &n, 16, 3, -1, SO3_STORE_ZERO_FIRST_PAD);
    assert( el == 2 && m == 1 && n == -1 &&
            "Index 16 yields wrong indices (el,m,n)." );

    so3_sampling_ind2elmn(&el, &m, &n, 25, 3, -1, SO3_STORE_ZERO_FIRST_PAD);
    assert( el == 2 && m == 1 && n == 1 &&
            "Index 25 yields wrong indices (el,m,n)." );

    // Test compact storage with n-order 0, -1, 1, -2, 2, ...
    // That is, for L = N = 2, the flmn are layed out as follows.
    // Each cell of the m-row corresponds to one coefficient in memory.
    //
    // n  |       0       |    -1     |     1     |
    // el | 0 |     1     |     1     |     1     |
    // m  | 0 |-1 | 0 | 1 |-1 | 0 | 1 |-1 | 0 | 1 |
    //
    // We pass in a nonsensical N to make sure implementation does
    // not depend on it.
    so3_sampling_ind2elmn(&el, &m, &n, 0, 3, -1, SO3_STORE_ZERO_FIRST_COMPACT);
    assert( el == 0 && m == 0 && n == 0 &&
            "Index 0 yields wrong indices (el,m,n)." );

    so3_sampling_ind2elmn(&el, &m, &n, 7, 3, -1, SO3_STORE_ZERO_FIRST_COMPACT);
    assert( el == 2 && m == 1 && n == 0 &&
            "Index 7 yields wrong indices (el,m,n)." );

    so3_sampling_ind2elmn(&el, &m, &n, 15, 3, -1, SO3_STORE_ZERO_FIRST_COMPACT);
    assert( el == 2 && m == 1 && n == -1 &&
            "Index 15 yields wrong indices (el,m,n)." );

    so3_sampling_ind2elmn(&el, &m, &n, 23, 3, -1, SO3_STORE_ZERO_FIRST_COMPACT);
    assert( el == 2 && m == 1 && n == 1 &&
            "Index 23 yields wrong indices (el,m,n)." );

    so3_sampling_ind2elmn(&el, &m, &n, 28, 3, -1, SO3_STORE_ZERO_FIRST_COMPACT);
    assert( el == 2 && m == 1 && n == -2 &&
            "Index 28 yields wrong indices (el,m,n)." );

    so3_sampling_ind2elmn(&el, &m, &n, 33, 3, -1, SO3_STORE_ZERO_FIRST_COMPACT);
    assert( el == 2 && m == 1 && n == 2 &&
            "Index 33 yields wrong indices (el,m,n)." );

    // Test padded storage with n-order ..., -2, -1, 0, 1, 2, ...
    // That is, for L = N = 2, the flmn are layed out as follows.
    // Each cell in the m? row corresponds to one coefficient in memory.
    // It contains an underscore for invalid el-m-n combinations
    // (which will be zeroes in memory) or the actual m-index for values
    // that will really be stored.
    //
    // n  |      -1       |       0       |       1       |
    // el | 0 |     1     | 0 |     1     | 0 |     1     |
    // m? | _ |-1 | 0 | 1 | 0 |-1 | 0 | 1 | _ |-1 | 0 | 1 |
    so3_sampling_ind2elmn(&el, &m, &n, 18, 3, 3, SO3_STORE_NEG_FIRST_PAD);
    assert( el == 0 && m == 0 && n == 0 &&
            "Index 18 yields wrong indices (el,m,n)." );

    so3_sampling_ind2elmn(&el, &m, &n, 4, 3, 3, SO3_STORE_NEG_FIRST_PAD);
    assert( el == 2 && m == -2 && n == -2 &&
            "Index 4 yields wrong indices (el,m,n)." );

    so3_sampling_ind2elmn(&el, &m, &n, 25, 3, 3, SO3_STORE_NEG_FIRST_PAD);
    assert( el == 2 && m == 1 && n == 0 &&
            "Index 25 yields wrong indices (el,m,n)." );

    so3_sampling_ind2elmn(&el, &m, &n, 16, 3, 3, SO3_STORE_NEG_FIRST_PAD);
    assert( el == 2 && m == 1 && n == -1 &&
            "Index 16 yields wrong indices (el,m,n)." );

    so3_sampling_ind2elmn(&el, &m, &n, 34, 3, 3, SO3_STORE_NEG_FIRST_PAD);
    assert( el == 2 && m == 1 && n == 1 &&
            "Index 34 yields wrong indices (el,m,n)." );

    // Test compact storage with n-order ..., -2, -1, 0, 1, 2, ...
    // That is, for L = N = 2, the flmn are layed out as follows.
    // Each cell of the m-row corresponds to one coefficient in memory.
    //
    // n  |    -1     |       0       |     1     |
    // el |     1     | 0 |     1     |     1     |
    // m  |-1 | 0 | 1 | 0 |-1 | 0 | 1 |-1 | 0 | 1 |
    so3_sampling_ind2elmn(&el, &m, &n, 13, 3, 3, SO3_STORE_NEG_FIRST_COMPACT);
    assert( el == 0 && m == 0 && n == 0 &&
            "Index 13 yields wrong indices (el,m,n)." );

    so3_sampling_ind2elmn(&el, &m, &n, 0, 3, 3, SO3_STORE_NEG_FIRST_COMPACT);
    assert( el == 2 && m == -2 && n == -2 &&
            "Index 0 yields wrong indices (el,m,n)." );

    so3_sampling_ind2elmn(&el, &m, &n, 20, 3, 3, SO3_STORE_NEG_FIRST_COMPACT);
    assert( el == 2 && m == 1 && n == 0 &&
            "Index 20 yields wrong indices (el,m,n)." );

    so3_sampling_ind2elmn(&el, &m, &n, 11, 3, 3, SO3_STORE_NEG_FIRST_COMPACT);
    assert( el == 2 && m == 1 && n == -1 &&
            "Index 11 yields wrong indices (el,m,n)." );

    so3_sampling_ind2elmn(&el, &m, &n, 28, 3, 3, SO3_STORE_NEG_FIRST_COMPACT);
    assert( el == 2 && m == 1 && n == 1 &&
            "Index 28 yields wrong indices (el,m,n)." );

    so3_sampling_ind2elmn(&el, &m, &n, 3, 3, 3, SO3_STORE_NEG_FIRST_COMPACT);
    assert( el == 2 && m == 1 && n == -2 &&
            "Index 3 yields wrong indices (el,m,n)." );

    so3_sampling_ind2elmn(&el, &m, &n, 33, 3, 3, SO3_STORE_NEG_FIRST_COMPACT);
    assert( el == 2 && m == 1 && n == 2 &&
            "Index 33 yields wrong indices (el,m,n)." );
}

void test_sampling_elmn2ind_real()
{
    int ind;
    // Test padded storage for a real signal (n-order 0, 1, 2, ...)
    // That is, for L = N = 2, the flmn are layed out as follows.
    // Each cell in the m? row corresponds to one coefficient in memory.
    // It contains an underscore for invalid el-m-n combinations
    // (which will be zeroes in memory) or the actual m-index for values
    // that will really be stored.
    //
    // n  |       0       |       1       |
    // el | 0 |     1     | 0 |     1     |
    // m? | 0 |-1 | 0 | 1 | _ |-1 | 0 | 1 |

    so3_sampling_elmn2ind_real(&ind, 0, 0, 0, 3, 3, SO3_STORE_ZERO_FIRST_PAD);
    assert( ind == 0 &&
            "Element (0,0,0) should be first." );

    so3_sampling_elmn2ind_real(&ind, 2, 1, 0, 3, 3, SO3_STORE_ZERO_FIRST_PAD);
    assert( ind == 7 &&
            "Element (2,1,0) in wrong position for L = 3." );

    so3_sampling_elmn2ind_real(&ind, 2, 1, 1, 3, 3, SO3_STORE_ZERO_FIRST_PAD);
    assert( ind == 16 &&
            "Element (2,1,1) in wrong position for L = 3." );

    // Test compact storage for a real signal (n-order 0, 1, 2, ...)
    // That is, for L = N = 2, the flmn are layed out as follows.
    // Each cell of the m-row corresponds to one coefficient in memory.
    //
    // n  |       0       |     1     |
    // el | 0 |     1     |     1     |
    // m  | 0 |-1 | 0 | 1 |-1 | 0 | 1 |
    //
    // We pass in a nonsensical N to make sure implementation does
    // not depend on it.
    so3_sampling_elmn2ind_real(&ind, 0, 0, 0, 3, 3, SO3_STORE_ZERO_FIRST_COMPACT);
    assert( ind == 0 &&
            "Element (0,0,0) should be first." );

    so3_sampling_elmn2ind_real(&ind, 2, 1, 0, 3, 3, SO3_STORE_ZERO_FIRST_COMPACT);
    assert( ind == 7 &&
            "Element (2,1,0) in wrong position for L = 3." );

    so3_sampling_elmn2ind_real(&ind, 2, 1, 1, 3, 3, SO3_STORE_ZERO_FIRST_COMPACT);
    assert( ind == 15 &&
            "Element (2,1,1) in wrong position for L = 3." );

    so3_sampling_elmn2ind_real(&ind, 2, 1, 2, 3, 3, SO3_STORE_ZERO_FIRST_COMPACT);
    assert( ind == 20 &&
            "Element (2,1,2) in wrong position for L = 3." );
}

void test_sampling_ind2elmn_real()
{
    int el, m, n;
    // Test padded storage for a real signal (n-order 0, 1, 2, ...)
    // That is, for L = N = 2, the flmn are layed out as follows.
    // Each cell in the m? row corresponds to one coefficient in memory.
    // It contains an underscore for invalid el-m-n combinations
    // (which will be zeroes in memory) or the actual m-index for values
    // that will really be stored.
    //
    // n  |       0       |       1       |
    // el | 0 |     1     | 0 |     1     |
    // m? | 0 |-1 | 0 | 1 | _ |-1 | 0 | 1 |
    //
    // We pass in a nonsensical N to make sure implementation does
    // not depend on it.
    so3_sampling_ind2elmn_real(&el, &m, &n, 0, 3, 3, SO3_STORE_ZERO_FIRST_PAD);
    assert( el == 0 && m == 0 && n == 0 &&
            "Index 0 yields wrong indices (el,m,n)." );

    so3_sampling_ind2elmn_real(&el, &m, &n, 7, 3, 3, SO3_STORE_ZERO_FIRST_PAD);
    assert( el == 2 && m == 1 && n == 0 &&
            "Index 7 yields wrong indices (el,m,n)." );

    so3_sampling_ind2elmn_real(&el, &m, &n, 16, 3, 3, SO3_STORE_ZERO_FIRST_PAD);
    assert( el == 2 && m == 1 && n == 1 &&
            "Index 16 yields wrong indices (el,m,n)." );

    // Test compact storage for a real signal (n-order 0, 1, 2, ...)
    // That is, for L = N = 2, the flmn are layed out as follows.
    // Each cell of the m-row corresponds to one coefficient in memory.
    //
    // n  |       0       |     1     |
    // el | 0 |     1     |     1     |
    // m  | 0 |-1 | 0 | 1 |-1 | 0 | 1 |
    //
    // We pass in a nonsensical N to make sure implementation does
    // not depend on it.
    so3_sampling_ind2elmn_real(&el, &m, &n, 0, 3, 3, SO3_STORE_ZERO_FIRST_COMPACT);
    assert( el == 0 && m == 0 && n == 0 &&
            "Index 0 yields wrong indices (el,m,n)." );

    so3_sampling_ind2elmn_real(&el, &m, &n, 7, 3, 3, SO3_STORE_ZERO_FIRST_COMPACT);
    assert( el == 2 && m == 1 && n == 0 &&
            "Index 7 yields wrong indices (el,m,n)." );

    so3_sampling_ind2elmn_real(&el, &m, &n, 15, 3, 3, SO3_STORE_ZERO_FIRST_COMPACT);
    assert( el == 2 && m == 1 && n == 1 &&
            "Index 15 yields wrong indices (el,m,n)." );

    so3_sampling_ind2elmn_real(&el, &m, &n, 20, 3, 3, SO3_STORE_ZERO_FIRST_COMPACT);
    assert( el == 2 && m == 1 && n == 2 &&
            "Index 23 yields wrong indices (el,m,n)." );
}
