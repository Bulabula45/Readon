// hashK.hpp

#ifndef HASHK_H
#define HASHK_H

#include <cstdint>
#include <string>

uint64_t hash64shift(uint64_t key);
std::string hash_64i(uint64_t key);



struct uint64_16t
{
    uint64_t a;
    uint16_t b;
};


struct align_res
{
    int t1_id;
    int t2_id;
    int hits;
    int second_hits;
    bool valid;
    
    align_res(int t1, int t2, int h, int sh, bool v)
        : t1_id(t1), t2_id(t2), hits(h), second_hits(sh), valid(v) {}

    
    align_res() : t1_id(0), t2_id(0), hits(0), second_hits(0), valid(false) {}
};

struct res_record
{
    align_res* align_res_array;
    bool flag;

    res_record(bool f = false) : align_res_array(new align_res[48]), flag(f) {}

    ~res_record() {
        delete[] align_res_array;
    }
};


struct f_align_res
{
    int q_id;
    int t1_id;
    int t2_id;
    int hits;
    int second_hits;
    int f_id;
    bool strand;
};

const uint64_t charToHash[] = {0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3,};
const char hashToChar[] = {'A', 'C', 'G', 'T'};
const uint64_t charToHash_r[] = {3, 0, 2, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,};


const uint64_t MIN_B_SIZE = 0x3FFF;
const uint64_t mask4 = 0xF, mask12 = 0xFFF, mask5 = 0b11111, mask11=0b11111111111;
const uint64_t mask20 = 0xFFFFF;
const off_t file_size_MAX = 1ULL << 30;
#endif // HASHK_H
