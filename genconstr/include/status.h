/* status.h: define keccak status operations */

#ifndef KECCAKSOLVER_STATUS_H
#define KECCAKSOLVER_STATUS_H

#include <algorithm>
#include <iostream>
#include <vector>
#include <assert.h>
#include <string.h>

#define keccak_var_num 800

/* macros: kidx
 * usage: given the index of variables in Keccak status array, computes
 *      its index.
 */
#define kidx(i, j, k) \
    ((((5 * i) + j) * 32) + k)

#define cbinom2(n) \
    ( ((n) * ((n)+1)) / 2)

#define deg2midx0(vnum) \
    (cbinom2(vnum) + (vnum))

#define deg2midx1(vnum, var1_idx) \
    (cbinom2(vnum) + (var1_idx))

#define deg2midx2(var1_idx, var2_idx) \
    ((var2_idx) > (var1_idx) ? (cbinom2(var2_idx) + (var1_idx)) : (cbinom2(var1_idx) + (var2_idx)))

#define extractbit(enc, idx) \
    ((enc >> idx) & 1)

#define ROR32(X, i) \
    (X.rotate(i))

// number of variables in quad and linear part
#define keccak_quad_size \
    ((keccak_var_num + 1) * keccak_var_num / 2)
#define keccak_lin_size \
    (keccak_var_num)
#define keccak_polynomial_size \
    (1 + keccak_lin_size + keccak_quad_size)
#define keccak_linear_system_size \
    (keccak_lin_size + 1)
#define keccak_status_linear_system_size \
    (keccak_linear_system_size * 32)

namespace Keccak {

/* convert big endian to little endian */
    inline uint32_t bitswap_32(uint32_t enc) {
        uint32_t ret = enc;
        uint32_t count = 31;

        enc >>= 1;
        while (enc) {
            ret <<= 1;
            ret |= enc & 1;
            enc >>= 1;
            count--;
        }
        ret <<= count;

        return ret;
    }

    typedef struct VariableRef {

        /* polynomial array */
        bool *polynomial;

        /* polynomial order; max:2 */
        uint64_t poly_ord;

        VariableRef();

        VariableRef(const VariableRef &src);

        ~VariableRef();

        /* add a monomial of order one */
        void add_variable_ord1(const int i, const int j, const int k);

        void add_variable_ord1(const int i);

        /* add a monomial of order two */
        void add_variable_ord2(const int i_1, const int j_1, const int k_1,
                               const int i_2, const int j_2, const int k_2);

        void add_variable_ord2(const int i_1, const int i_2);

        /* update the order of polynomial after operations */
        void update_poly_ext_ord();

        /* set all coef in polynomial to 0 */
        void clear();

        /*
         * using known variable values to eliminate the stored polynomial
         * input:
         *    determined_system: integer array with length keccak_var_num
         *                       0,1: known value
         *                       -1: unknown flag
         */
        void elimination(int *determined_system);

        /* misc pointers */
        bool *quadratic_part();

        bool *linear_part();

        bool *const_part();

        /* for debug use, print the polynomial in deg 2 grlex order */
        void display(int sidx = 0, int lidx = keccak_polynomial_size);

        /* assign by copy */
        const VariableRef &operator=(const VariableRef &other);

        /* add two variables together */
        VariableRef operator+(const VariableRef &other);

        /* add this variable with an constant */
        VariableRef operator+(const int c);

        /* multiply two variables together */
        VariableRef operator*(VariableRef &other);

        /* multiply a variable with a constant */
        VariableRef operator*(const int c);

    } variable_ref_t;

    typedef struct StatusEpi32 {

        /* variable array */
        variable_ref_t var32[32];

        /* dynamic constraint array */
        std::vector<variable_ref_t> constr;

        StatusEpi32();

        StatusEpi32(const StatusEpi32 &src);

        ~StatusEpi32();

        /* assign by copy */
        const StatusEpi32 &operator=(const StatusEpi32 &other);

        /* assign constance */
        StatusEpi32 operator=(const uint32_t c);

        /* rotate permutation by c */
        StatusEpi32 rotate(const uint32_t c);

        /* add two status rows together */
        StatusEpi32 operator^(const StatusEpi32 &other);

        /* add a row with a constance */
        StatusEpi32 operator^(const uint32_t c);

        /* multiply two row */
        StatusEpi32 operator&(StatusEpi32 &other);

        /* multiply a row with a constance */
        StatusEpi32 operator&(const uint32_t c);

        /* set variable level constraint
         * i.e. A[i_0][j_0] == A[i_1][j_1]
         * store the result as an equation, with the lhs
         *    is represented by a variable_ref_t, and rhs
         *    equals to zero
         * */
        void set_equiv_constraint(const uint32_t c);

        /* set variable versus constance constraint
         * i.e. A[i_0][j_0] == 0xFFFFFFFF
         * store the result as an equation, with the lhs
         *    is represented by a variable_ref_t, and rhs
         *    equals to zero
         * */
        void set_equiv_constraint(const StatusEpi32 &other);

        /* get the number of linear and quadratic constraints */
        uint32_t get_linear_constraint_num();

        uint32_t get_quad_constraint_num();

        /* copy the linear and quadratic equations to destination
         *    returns the number of constraints respectively
         * */
        uint32_t dump_linear_constraint(bool *system);

        uint32_t dump_quad_constraint(bool *system);

        /* assign each variable a index
         * call this function when startup
         * */
        void init_system(int x_index, int y_index);

        /* for debug use, print requested variable in deg 2 grlex order */
        void display_var(uint64_t vindex);

        void display_var_linpart(uint64_t vindex);

        void display_var_quadpart(uint64_t vindex);

        void display_all_var_quadpart();

        void display_all_var();

        void display_all_var_linpart();

        /* for debug use, print constraints in deg 2 grlex order */
        void display_constr();
    } status_epi32_t;

    typedef status_epi32_t Status;

    void init_status(status_epi32_t status[5][5]);

    void display_status(status_epi32_t status[5][5]);

    void display_status_quadpart(status_epi32_t status[5][5]);

    void display_status_linpart(status_epi32_t status[5][5]);

    void status_set_val(Status &status, uint32_t value);
}
#endif //KECCAKSOLVER_STATUS_H
