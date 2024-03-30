#include "search.hpp"
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <iostream>
#include <mutex>

using std::cout;
using std::endl;

inline int getHighFreq(int *lst, int l, int &up, int &down, int &second_hits, bool mode = false)
{
    std::unordered_map<int, int> freqMap;

    for (int i = 0; i < l; ++i)
    {
        freqMap[lst[i]]++;
    }

    std::vector<std::pair<int, int>> freqPairs;
    for (const auto &entry : freqMap)
    {
        freqPairs.push_back(std::make_pair(entry.first, entry.second));
    }

    std::sort(freqPairs.begin(), freqPairs.end(), [](const auto &a, const auto &b)
              { return a.second > b.second; });

    if (mode)
    {
        cout << "freqPairs: " << freqPairs.size() << endl;
        for (size_t i = 0; i < freqPairs.size(); ++i)
        {
            cout << freqPairs[i].first << " " << freqPairs[i].second << endl;
        }
    }

    if (freqPairs.size() < 2)
        return 0;

    else if (freqPairs.size() == 2)
    {
        if (freqPairs[0].first != 0)
        {
            up = freqPairs[0].first;
            down = freqPairs[1].first;
            second_hits = freqPairs[1].second;
            return freqPairs[0].second;
        }
        else
        {
            up = freqPairs[1].first;
            down = freqPairs[0].first;
            second_hits = freqPairs[0].second;
            return freqPairs[1].second;
        }
    }

    else
    {
        if (freqPairs[0].first == 0)
        {
            up = freqPairs[1].first;
            down = freqPairs[2].first;
            second_hits = freqPairs[2].second;
            return freqPairs[1].second;
        }
        else if (freqPairs[1].first == 0)
        {
            up = freqPairs[0].first;
            down = freqPairs[2].first;
            second_hits = freqPairs[2].second;
            return freqPairs[0].second;
        }
        else
        {
            up = freqPairs[0].first;
            down = freqPairs[1].first;
            second_hits = freqPairs[1].second;
            return freqPairs[0].second;
        }
    }
}

void searchRefs(uint64_t **miQQuery, int *query_mi_cnts, int query_cnt, const char *reFnames[], const int recSize[],
                int p_start, int p_end,
                res_record *res, res_record *res_r)
{

    const int MAX_DIST = 300000;
    const int Thr_hits = kmer_l > 20 ? 6 : 10;

    for (int f_id = p_start; f_id < p_end; f_id++)
    {
        const char *ref_name = reFnames[f_id];
        int rec_size = recSize[f_id];
        uint64_16t *miRef = new uint64_16t[rec_size * int(2250 / window_s)];
        uint64_16t *miRef_r = new uint64_16t[rec_size * int(2250 / window_s)];
        uint32_t *mi_ptr = new uint32_t[int(rec_size * 0.2)];

        int fd = open(ref_name, O_RDONLY);
        int len = lseek(fd, 0, SEEK_END);
        char *mbuf = (char *)mmap(NULL, len, PROT_READ, MAP_PRIVATE, fd, 0);
        char *ss = mbuf;
        char *ee = mbuf + len;
        char *seq_s = ss, *seq_e = ss;
        int seq_count = 0;
        int mi_count = 0;
        int ref_s = 0, ref_e = 0, ls_e = 0;
        uint16_t clstr = 0;
        uint32_t mi_s = 0, mi_e = 0;
        int bin = 0;
        while (seq_count < rec_size)
        {
            if (*ss == '>')
            {
                ss++;
                while (*ss != '|')
                {
                    ss++;
                }
                ss++;

                ref_s = 0;
                while (*ss != '|')
                {
                    ref_s = ref_s * 10 + (*ss - '0');
                    ss++;
                }
                ss++;

                ref_e = 0;
                while (*ss != '|')
                {
                    ref_e = ref_e * 10 + (*ss - '0');
                    ss++;
                }
                ss++;

                if (ref_s - ls_e > MAX_DIST)
                {

                    mi_ptr[bin] = mi_s;
                    bin++;
                    mi_ptr[bin] = mi_e;
                    bin++;

                    mi_s = mi_count;
                }

                ls_e = ref_e;

                clstr = 0;
                while (*ss != '\n')
                {
                    clstr = clstr * 10 + (*ss - '0');
                    ss++;
                }
                clstr = static_cast<uint16_t>(clstr + 1);
                ss++;
            }

            seq_s = ss;
            while (*ss != '\n' && ss < ee)
            {
                ss++;
            }
            ss++;
            seq_e = ss;
            int size_seq = seq_e - seq_s;
            uint64_t kmerKey = 0, baseCode = 0, min_h = 0xFFFFFFFFFFFFFFFF, tmp_h_key = 0;
            uint64_t kmerKey_r = 0, baseCode_r = 0, min_h_r = 0xFFFFFFFFFFFFFFFF, tmp_h_key_r = 0;
            int tmp = 0;

            for (int i = 0; i < kmer_l; ++i)
            {
                baseCode = charToHash[*(i + seq_s) - 'A'];
                baseCode_r = charToHash_r[*(seq_e - i - 2) - 'A'];

                kmerKey = ((kmerKey << 2) & ((1ULL << kmer_2l) - 1)) | baseCode;
                kmerKey_r = ((kmerKey_r << 2) & ((1ULL << kmer_2l) - 1)) | baseCode_r;
            }

            miRef[mi_count] = {hash64shift(kmerKey), clstr};
            miRef_r[mi_count] = {hash64shift(kmerKey_r), clstr};
            mi_count++;

            for (int i = kmer_l; i < size_seq; ++i)
            {
                baseCode = charToHash[*(i + seq_s) - 'A'];
                baseCode_r = charToHash_r[*(seq_e - i - 2) - 'A'];

                kmerKey = ((kmerKey << 2) & ((1ULL << kmer_2l) - 1)) | baseCode;
                kmerKey_r = ((kmerKey_r << 2) & ((1ULL << kmer_2l) - 1)) | baseCode_r;
                tmp = (i + 1) % window_s;
                if (tmp == 0)
                {
                    min_h = hash64shift(kmerKey);
                    min_h_r = hash64shift(kmerKey_r);
                }
                else if (tmp <= interval)
                {
                    tmp_h_key = hash64shift(kmerKey);
                    tmp_h_key_r = hash64shift(kmerKey_r);
                    if (tmp_h_key < min_h)
                    {
                        min_h = tmp_h_key;
                    }
                    if (tmp_h_key_r < min_h_r)
                    {
                        min_h_r = tmp_h_key_r;
                    }
                    if (tmp == interval || i == size_seq - 1)
                    {
                        miRef[mi_count] = {min_h, clstr};
                        miRef_r[mi_count] = {min_h_r, clstr};
                        mi_count++;
                    }
                }
            }
            mi_e = mi_count;
            seq_count++;
        }
        munmap(mbuf, len);
        close(fd);

        mi_ptr[bin] = mi_ptr[bin - 1];
        bin++;
        mi_ptr[bin] = mi_count;
        bin++;

        uint32_t *new_mi_ptr = new uint32_t[int(rec_size * 0.2)];
        int start_j = 0, new_b = 0, j = 0;
        if (mi_ptr[1] == 0)
        {
            start_j = 2;
            j = 2;
        }
        for (; j < bin / 2; j++)
        {
            if (mi_ptr[2 * j + 1] - mi_ptr[start_j] > 32767)
            {
                new_mi_ptr[new_b] = mi_ptr[start_j];
                new_b++;
                new_mi_ptr[new_b] = mi_ptr[2 * j - 1];
                new_b++;
                start_j = 2 * j;
            }
        }
        new_mi_ptr[new_b] = mi_ptr[start_j];
        new_b++;
        new_mi_ptr[new_b] = mi_ptr[bin - 1];
        new_b++;
        delete[] mi_ptr;

        uint64_t mo_number = MIN_B_SIZE > static_cast<uint64_t>(mi_count) ? MIN_B_SIZE : static_cast<uint64_t>(mi_count);

        for (int i = 0; i < new_b / 2; i++)
        {

            mi_s = new_mi_ptr[2 * i];
            mi_e = new_mi_ptr[2 * i + 1];

            uint64_t *miHashT = new uint64_t[mo_number];
            uint64_t *miHashT_r = new uint64_t[mo_number];
            for (size_t i = 0; i < mo_number; i++)
            {
                miHashT[i] = 0;
                miHashT_r[i] = 0;
            }
            uint64_t h_mi = 0;
            uint64_t cur = 0;
            uint64_16t mi_clstr = {0UL, 0};
            uint16_t cur_clstr = 0;
            for (uint32_t i = mi_s; i < mi_e; i++)
            {
                mi_clstr = miRef[i];

                cur = miHashT[h_mi];
                h_mi = mi_clstr.a % mo_number;

                cur_clstr = mi_clstr.b;
                if (cur == 0)
                {

                    miHashT[h_mi] = (((mi_clstr.a & mask5) << 11) | cur_clstr);
                }
                else if ((cur && mask11) == cur_clstr)
                {

                    continue;
                }
                else
                {
                    miHashT[h_mi] <<= 16;
                    miHashT[h_mi] |= (((mi_clstr.a & mask5) << 11) | cur_clstr);
                }

                mi_clstr = miRef_r[i];
                cur = miHashT_r[h_mi];
                h_mi = mi_clstr.a % mo_number;
                cur_clstr = mi_clstr.b;
                if (cur == 0)
                {

                    miHashT_r[h_mi] = (((mi_clstr.a & mask5) << 11) | cur_clstr);
                }
                else if ((cur && mask11) == cur_clstr)
                {

                    continue;
                }
                else
                {
                    miHashT_r[h_mi] <<= 16;
                    miHashT_r[h_mi] = (((mi_clstr.a & mask5) << 11) | cur_clstr);
                }
            }

            for (int q_id = 0; q_id < query_cnt; q_id++)
            {

                int query_mi_count = query_mi_cnts[q_id];
                uint64_t *miQuery = miQQuery[q_id];
                int *clstrs = new int[query_mi_count];
                int *clstrs_r = new int[query_mi_count];
                int shift_n = 0;
                for (int i = 0; i < query_mi_count; i++)
                {

                    uint64_t qmi = miQuery[i];

                    h_mi = qmi % mo_number;

                    cur = miHashT[h_mi];

                    // if (cur == 0 || ((cur >> 12) != (qmi & mask20)))
                    clstrs[i] = 0;
                    shift_n = 0;
                    while (cur != 0 && shift_n < 4)
                    {
                        if ((cur >> 11) == (qmi & mask5))
                        {
                            clstr = cur & mask11;
                            clstrs[i] = clstr;
                        }
                        else
                        {
                            cur >>= 16;
                        }
                        shift_n++;
                    }

                    cur = miHashT_r[h_mi];
                    clstrs_r[i] = 0;
                    shift_n = 0;
                    while (cur != 0 && shift_n < 4)
                    {
                        if ((cur >> 11) == (qmi & mask5))
                        {
                            clstr = cur & mask11;
                            clstrs_r[i] = clstr;
                        }
                        else
                        {
                            cur >>= 16;
                        }
                        shift_n++;
                    }
                }

                int up = 0, down = 0, up_r = 0, down_r = 0;
                int second_hits = 0;
                int hit_nums = 0;

                hit_nums = getHighFreq(clstrs, query_mi_count, up, down, second_hits);

                if (hit_nums > Thr_hits && second_hits > 1)
                {
                    if (!(up - down <= -10 || up - down >= 10))
                    {
                        res[q_id].flag = true;
                        res[q_id].align_res_array[f_id] = {up, down, hit_nums, second_hits, true};
                    }
                }
                delete[] clstrs;

                hit_nums = getHighFreq(clstrs_r, query_mi_count, up_r, down_r, second_hits);
                if (hit_nums > Thr_hits && second_hits > 1)
                {
                    if (!(up_r - down_r <= -10 || up_r - down_r >= 10))
                    {
                        res_r[q_id].flag = true;
                        res_r[q_id].align_res_array[f_id] = {up_r, down_r, hit_nums, second_hits, true};
                    }
                }
                delete[] clstrs_r;
            }

            delete[] miHashT;
            delete[] miHashT_r;
        }
        delete[] new_mi_ptr;
        delete[] miRef;
        delete[] miRef_r;
    }
}
