/* status.h: define keccak status operations */

#include "status.h"

namespace Keccak {

VariableRef::VariableRef() {
    polynomial = (bool *) malloc(sizeof(bool) * keccak_polynomial_size);
    clear();
}

VariableRef::VariableRef(const Keccak::VariableRef &src) {
    polynomial = (bool *) malloc(sizeof(bool) * keccak_polynomial_size);
    clear();
    std::copy(src.polynomial, src.polynomial + keccak_polynomial_size, polynomial);
    this->poly_ord = src.poly_ord;
}

VariableRef::~VariableRef() {
    if (polynomial != nullptr)
        free(polynomial);
}

void VariableRef::add_variable_ord1(const int i, const int j, const int k) {
    int var_index = deg2midx1(keccak_var_num, kidx(i, j, k));
    polynomial[var_index] ^= 1;
    update_poly_ext_ord();
}

void VariableRef::add_variable_ord1(const int i) {
    polynomial[keccak_quad_size + i] ^= 1;
    update_poly_ext_ord();
}

void VariableRef::add_variable_ord2(const int i_1, const int j_1, const int k_1,
                                    const int i_2, const int j_2, const int k_2) {
    int var_index = deg2midx2(kidx(i_1, j_1, k_1), kidx(i_2, j_2, k_2));
    polynomial[var_index] ^= 1;
    update_poly_ext_ord();
}

void VariableRef::add_variable_ord2(const int i_1, const int i_2) {
    int var_index = deg2midx2(i_1, i_2);
    polynomial[var_index] ^= 1;
    update_poly_ext_ord();
}

void VariableRef::update_poly_ext_ord() {
    this->poly_ord = 0;

    // check quad part
    for (int i = 0; i < keccak_quad_size; i++) {
        if (polynomial[i]) {
            this->poly_ord = 2;
            return;
        }
    }

    // check linear part
    for (int i = keccak_quad_size; i < keccak_quad_size + keccak_lin_size; i++) {
        if (polynomial[i]) {
            this->poly_ord = 1;
            return;
        }
    }
}

void VariableRef::clear() {
    for (int i = 0; i < keccak_polynomial_size; i++)
        this->polynomial[i] = 0;
    this->poly_ord = 0;
}

void VariableRef::elimination(int *determined_system) {
    for (int i = 0; i < keccak_var_num; i++) {
        if (determined_system[i] == -1)
            continue;

        uint64_t subpoly_start_index = (uint64_t) cbinom2(i);
        uint64_t subpoly_end_index = (uint64_t) cbinom2(i + 1);
        if (determined_system[i] == 1) {
            for (uint64_t j = subpoly_start_index; j < subpoly_end_index; j++) {
                if (this->polynomial[j]) {
                    uint64_t varindex = i + j - subpoly_start_index;
                    this->add_variable_ord1(varindex);
                    this->polynomial[j] = 0;
                }
            }
        } else {
            for (uint64_t j = subpoly_start_index; j < subpoly_end_index; j++)
                this->polynomial[j] = 0;
        }
    }
    this->update_poly_ext_ord();
}

bool *VariableRef::quadratic_part() {
    return &polynomial[0];
}

bool *VariableRef::linear_part() {
    return &polynomial[keccak_quad_size];
}

bool *VariableRef::const_part() {
    return &polynomial[keccak_quad_size + keccak_lin_size];
}

const VariableRef &VariableRef::operator=(const Keccak::VariableRef &other) {
    if (this == &other) return *this;
    std::copy(other.polynomial, other.polynomial + keccak_polynomial_size, this->polynomial);
    this->poly_ord = other.poly_ord;
    return *this;
}

VariableRef VariableRef::operator+(const Keccak::VariableRef &other) {
    VariableRef v;
    for (uint64_t i = 0; i < keccak_polynomial_size; i++)
        v.polynomial[i] = this->polynomial[i] ^ other.polynomial[i];
    v.update_poly_ext_ord();
    return v;
}

VariableRef VariableRef::operator+(const int c) {
    VariableRef v;
    std::copy(this->polynomial, this->polynomial + keccak_polynomial_size, v.polynomial);
    v.polynomial[keccak_polynomial_size - 1] ^= c;
    v.update_poly_ext_ord();
    return v;
}

VariableRef VariableRef::operator*(Keccak::VariableRef &other) {
    // only handle quadratic product
    assert(this->poly_ord + other.poly_ord <= 2);

    VariableRef v;
    if (this->poly_ord == 0 && (*(this->const_part()))) { // const 1 mul const
        v = other;
    } else if (this->poly_ord == 1) { // linear mul other
        bool *other_linear_part = other.linear_part();
        bool *other_const_part = other.const_part();
        bool *this_linear_part = this->linear_part();
        bool *this_const_part = this->const_part();

        int thisIdx;
        int otherIdx;

        for (thisIdx = 0; thisIdx < keccak_lin_size; thisIdx++)
            if (this_linear_part[thisIdx])
                for (otherIdx = 0; otherIdx < keccak_lin_size; otherIdx++)
                    if (other_linear_part[otherIdx])
                        v.add_variable_ord2(thisIdx, otherIdx);

        if (*this_const_part)
            for (otherIdx = 0; otherIdx < keccak_lin_size; otherIdx++)
                if (other_linear_part[otherIdx])
                    v.add_variable_ord1(otherIdx);

        if (*other_const_part)
            for (thisIdx = 0; thisIdx < keccak_lin_size; thisIdx++)
                if (this_linear_part[thisIdx])
                    v.add_variable_ord1(thisIdx);

        if (*other_const_part && *this_const_part)
            v = v + 1;

    } else if (this->poly_ord == 2) { // quadratic mul other
        if (*other.const_part())
            v = *this;
    }

    v.update_poly_ext_ord();

    return v;
}

VariableRef VariableRef::operator*(const int c) {
    VariableRef v;
    v.clear();
    if (c != 0)
        v = *this;
    return v;
}

void VariableRef::display(int sidx, int lidx) {
    for (int i = sidx; i < lidx; i++)
        std::cout << (int) polynomial[i];
    std::cout << std::endl;
}

StatusEpi32::StatusEpi32() {
}

StatusEpi32::StatusEpi32(const Keccak::StatusEpi32 &src) {
    *this = src;
    this->constr = src.constr;
}

StatusEpi32::~StatusEpi32() {
}

const StatusEpi32 &StatusEpi32::operator=(const Keccak::StatusEpi32 &other) {
    if (this == &other) return *this;
//    this->constr = other.constr;
    std::copy(other.var32, other.var32 + 32, this->var32);
    return *this;
}

StatusEpi32 StatusEpi32::operator=(const uint32_t c) {
    StatusEpi32 s;
    for (int i = 0; i < 32; i++)
        *(s.var32[i].const_part()) = extractbit(c, i);
    return s;
}

// TODO: design efficient procedures to speedup rotate operations
StatusEpi32 StatusEpi32::rotate(const uint32_t c) {
    StatusEpi32 s;
    s = *this;
    for (int i = 0; i < 32; i++)
        s.var32[i] = this->var32[(i + c) % 32];
    return s;
}

StatusEpi32 StatusEpi32::operator^(const Keccak::StatusEpi32 &other) {
    StatusEpi32 s;
    s = *this;
    for (int i = 0; i < 32; i++)
        s.var32[i] = this->var32[i] + other.var32[i];
    return s;
}

StatusEpi32 StatusEpi32::operator^(const uint32_t c) {
//    uint32_t enc = bitswap_32(c);
    StatusEpi32 s;
    s = *this;
    for (int i = 0; i < 32; i++)
        s.var32[i] = this->var32[i] + extractbit(c, i);
    return s;
}

StatusEpi32 StatusEpi32::operator&(Keccak::StatusEpi32 &other) {
    StatusEpi32 s;
    s = *this;
    for (int i = 0; i < 32; i++)
        s.var32[i] = this->var32[i] * other.var32[i];
    return s;
}

StatusEpi32 StatusEpi32::operator&(const uint32_t c) {
//    uint32_t enc = bitswap_32(c);
    StatusEpi32 s;
    s = *this;
    for (int i = 0; i < 32; i++) {
        s.var32[i] = this->var32[i] * extractbit(c, i);
    }
    return s;
}

void StatusEpi32::set_equiv_constraint(const uint32_t c) {
//    uint32_t enc = bitswap_32(c);
    for (int i = 0; i < 32; i++) {
        VariableRef constrRef = var32[i] + extractbit(c, i);
        if (constrRef.poly_ord != 0)
            constr.push_back(constrRef);
    }
}
void StatusEpi32::set_equiv_constraint(const Keccak::StatusEpi32 &other) {
    for (int i = 0; i < 32; i++) {
        VariableRef vRef = var32[i] + other.var32[i];
        if (vRef.poly_ord != 0)
            constr.push_back(vRef);
    }
}

uint32_t StatusEpi32::get_linear_constraint_num() {
    uint32_t ret = 0;
    for (uint32_t i = 0; i < constr.size(); i++)
        if (constr[i].poly_ord == 1)
            ret++;
    return ret;
}

uint32_t StatusEpi32::get_quad_constraint_num() {
    uint32_t ret = 0;
    for (uint32_t i = 0; i < constr.size(); i++)
        if (constr[i].poly_ord == 2)
            ret++;
    return ret;
}

uint32_t StatusEpi32::dump_linear_constraint(bool *system) {
    uint32_t ret = 0;
    uint32_t offset = 0;
    for (uint32_t i = 0; i < constr.size(); i++)
        if (constr[i].poly_ord == 1) {
            std::copy(constr[i].linear_part(), constr[i].linear_part() + keccak_linear_system_size, system + offset);
            offset += keccak_linear_system_size;
            ret++;
        }
    return ret;
}

uint32_t StatusEpi32::dump_quad_constraint(bool *system) {
    uint32_t ret = 0;
    uint32_t offset = 0;
    for (uint32_t i = 0; i < constr.size(); i++)
        if (constr[i].poly_ord == 2) {
            std::copy(constr[i].quadratic_part(), constr[i].quadratic_part() + keccak_polynomial_size, system + offset);
            offset += keccak_polynomial_size;
            ret++;
        }
    return ret;
}

void StatusEpi32::init_system(int x_index, int y_index) {
    for (uint64_t i = 0; i < 32; i++) {
        this->var32[i].clear();
        this->var32[i].add_variable_ord1(x_index, y_index, i);
    }
}

void StatusEpi32::display_var(uint64_t vindex) {
    var32[vindex].display();
}

void StatusEpi32::display_var_linpart(uint64_t vindex) {
    var32[vindex].display(keccak_quad_size, keccak_polynomial_size);
}

void StatusEpi32::display_var_quadpart(uint64_t vindex) {
    var32[vindex].display(0, keccak_quad_size);
}

void StatusEpi32::display_all_var() {
    for (int i = 0; i < 32; i++)
        var32[i].display();
}

void StatusEpi32::display_all_var_linpart() {
    for (int i = 0; i < 32; i++)
        var32[i].display(keccak_quad_size, keccak_polynomial_size);
}

void StatusEpi32::display_all_var_quadpart() {
    for (int i = 0; i < 32; i++)
        var32[i].display(0, keccak_quad_size);
}

void StatusEpi32::display_constr() {
    for (uint64_t i = 0; i < constr.size(); i++)
        constr[i].display(keccak_quad_size, keccak_polynomial_size);
}

void init_status(status_epi32_t status[5][5]) {
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            status[i][j].init_system(i, j);
        }
    }
}

void display_status(status_epi32_t status[5][5]) {
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            printf("Display status[%d][%d]:\n", i, j);
            for (int k = 0; k < 32; k++)
                status[i][j].display_var(k);
            printf("\n\n");
        }
    }
}

void display_status_quadpart(status_epi32_t status[5][5]) {
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            printf("Display status[%d][%d]:\n", i, j);
            for (int k = 0; k < 32; k++)
                status[i][j].display_var_quadpart(k);
            printf("\n\n");
        }
    }
}

void display_status_linpart(status_epi32_t status[5][5]) {
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            printf("Display status[%d][%d]:\n", i, j);
            for (int k = 0; k < 32; k++)
                status[i][j].display_var_linpart(k);
            printf("\n\n");
        }
    }
}

void status_set_val(Status &status, uint32_t value) {
    for (int i = 0; i < 32; i++) {
        *(status.var32[i].const_part()) = extractbit(value, i);
        status.var32[i].update_poly_ext_ord();
    }
}
}