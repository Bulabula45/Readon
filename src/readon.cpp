#include "readon.hpp"
#include <iostream>
#include <thread>

using std::cout;
using std::endl;

size_t getPageSize()
{
    return sysconf(_SC_PAGE_SIZE);
}

size_t pageSize = getPageSize();

size_t alignToPageSize(size_t value, size_t pageSize)
{
    return (value - pageSize - 1) & ~(pageSize - 1);
}

off_t mapping(const int reads_cnt, int linesPerRecord, off_t reads_len, off_t offset, int reads_fd, std::ofstream &outputTsv, std::ofstream &outputReads)

{

    cout << "reads_len: " << reads_len << endl;
    cout << "offset: " << offset << endl;
    off_t offset_align = offset > 0 ? alignToPageSize(offset, pageSize) : offset;
    cout << "offset_align: " << offset_align << endl;
    cout << "reads_cnt: " << reads_cnt << endl;
    reads_len += offset - offset_align;
    char *reads_mbuf = (char *)mmap(NULL, reads_len, PROT_READ, MAP_PRIVATE, reads_fd, offset_align);
    if (reads_mbuf == MAP_FAILED)
    {
        perror("mmap");
        exit(1);
    }
    char *reads_ss = reads_mbuf + (offset - offset_align);

    char *reads_ee = reads_mbuf + reads_len;
    char *query_s = reads_ss, *header_s = reads_ss, *query_e = reads_ss;
    uint64_t **miQuery = new uint64_t *[reads_cnt];
    uint8_t **seqQuery = new uint8_t *[reads_cnt];
    int *sizeQuery = new int[reads_cnt];
    char **query_names = new char *[reads_cnt];
    int *query_mi_cnts = new int[reads_cnt];
    int q_cnt = 0;

    if (linesPerRecord == 2)
    {
        char hder = '>';
        while (q_cnt < reads_cnt && reads_ss < reads_ee)
        {
            if (*reads_ss == hder)
            {
                reads_ss++;
                header_s = reads_ss;
                while (*reads_ss != '\n')
                {
                    reads_ss++;
                }
                query_names[q_cnt] = new char[reads_ss - header_s + 1];
                memcpy(query_names[q_cnt], header_s, reads_ss - header_s);
                query_names[q_cnt][reads_ss - header_s] = '\0';
            }
            reads_ss++;
            query_s = reads_ss;
            while (*reads_ss != '\n' && reads_ss < reads_ee)
            {
                reads_ss++;
            }
            query_e = reads_ss;
            reads_ss++;

            int size_seq = query_e - query_s;
            int bytes_seq = (size_seq / 4) + (size_seq % 4 != 0);
            int query_mi_count = 0;
            miQuery[q_cnt] = new uint64_t[size_seq / window_s + 1];
            seqQuery[q_cnt] = new uint8_t[bytes_seq]();
            sizeQuery[q_cnt] = size_seq;
            uint8_t shift = 0;

            query_mi_count = 0;
            uint64_t kmerKey = 0, baseCode = 0, min_h = 0xFFFFFFFFFFFFFFFF, tmp_h_key = 0;
            int tmp = 0;
            for (int i = 0; i < kmer_l; ++i)
            {
                baseCode = charToHash[*(i + query_s) - 'A'];
                kmerKey = ((kmerKey << 2) & ((1ULL << kmer_2l) - 1)) | baseCode;

                shift = 6 - 2 * (i % 4);
                seqQuery[q_cnt][i / 4] |= baseCode << shift;
            }
            miQuery[q_cnt][query_mi_count] = hash64shift(kmerKey);
            query_mi_count++;

            for (int i = kmer_l; i < size_seq; ++i)
            {
                baseCode = charToHash[*(i + query_s) - 'A'];

                shift = 6 - 2 * (i % 4);
                seqQuery[q_cnt][i / 4] |= baseCode << shift;

                kmerKey = ((kmerKey << 2) & ((1ULL << kmer_2l) - 1)) | baseCode;
                tmp = (i + 1) % window_s;
                if (tmp == 0)
                {
                    min_h = hash64shift(kmerKey);
                }
                else if (tmp <= interval)
                {
                    tmp_h_key = hash64shift(kmerKey);
                    if (tmp_h_key < min_h)
                    {
                        min_h = tmp_h_key;
                    }
                    if (tmp == interval)
                    {
                        miQuery[q_cnt][query_mi_count] = min_h;
                        query_mi_count++;
                    }
                }
            }
            query_mi_cnts[q_cnt] = query_mi_count;
            q_cnt++;
        }
    }
    else
    {
        char hder = '@';
        while (q_cnt < reads_cnt && reads_ss < reads_ee)
        {
            if (*reads_ss == hder)
            {
                reads_ss++;
                header_s = reads_ss;
                while (*reads_ss != ' ')
                {
                    reads_ss++;
                }
                query_names[q_cnt] = new char[reads_ss - header_s + 1];
                memcpy(query_names[q_cnt], header_s, reads_ss - header_s);
                query_names[q_cnt][reads_ss - header_s] = '\0';
            }
            reads_ss++;

            while (*reads_ss != '\n')
            {
                reads_ss++;
            }
            reads_ss++;
            query_s = reads_ss;
            while (*reads_ss != '\n')
            {
                reads_ss++;
            }
            query_e = reads_ss;
            reads_ss++;

            int size_seq = query_e - query_s;
            int query_mi_count = 0;
            int bytes_seq = (size_seq / 4) + (size_seq % 4 != 0);
            miQuery[q_cnt] = new uint64_t[size_seq / window_s + 1];
            seqQuery[q_cnt] = new uint8_t[bytes_seq];
            sizeQuery[q_cnt] = size_seq;

            query_mi_count = 0;
            uint64_t kmerKey = 0, baseCode = 0, min_h = 0xFFFFFFFFFFFFFFFF, tmp_h_key = 0;
            int tmp = 0;
            uint8_t shift = 0;
            for (int i = 0; i < kmer_l; ++i)
            {
                baseCode = charToHash[*(i + query_s) - 'A'];
                kmerKey = ((kmerKey << 2) & ((1ULL << kmer_2l) - 1)) | baseCode;

                shift = 6 - 2 * (i % 4);
                seqQuery[q_cnt][i / 4] |= baseCode << shift;
            }
            miQuery[q_cnt][query_mi_count] = hash64shift(kmerKey);
            query_mi_count++;

            for (int i = kmer_l; i < size_seq; ++i)
            {
                baseCode = charToHash[*(i + query_s) - 'A'];

                shift = 6 - 2 * (i % 4);
                seqQuery[q_cnt][i / 4] |= baseCode << shift;

                kmerKey = ((kmerKey << 2) & ((1ULL << kmer_2l) - 1)) | baseCode;
                tmp = (i + 1) % window_s;
                if (tmp == 0)
                {
                    min_h = hash64shift(kmerKey);
                }
                else if (tmp <= interval)
                {
                    tmp_h_key = hash64shift(kmerKey);
                    if (tmp_h_key < min_h)
                    {
                        min_h = tmp_h_key;
                    }
                    if (tmp == interval)
                    {
                        miQuery[q_cnt][query_mi_count] = min_h;
                        query_mi_count++;
                    }
                }
            }
            query_mi_cnts[q_cnt] = query_mi_count;
            q_cnt++;

            while (*reads_ss != '\n')
            {
                reads_ss++;
            }
            reads_ss++;
            while (*reads_ss != '\n' && reads_ss < reads_ee)
            {
                reads_ss++;
            }
            reads_ss++;
        }
    }

    off_t curPosition = offset_align + (reads_ss - reads_mbuf);

    munmap(reads_mbuf, reads_len);
    cout << "load " << q_cnt << " reads done." << endl;

    res_record *res = new res_record[q_cnt];
    res_record *res_r = new res_record[q_cnt];
    int threads_num = 12;
    int file_per_thread = file_end / threads_num;
    std::vector<std::thread> works;
    works.reserve(threads_num);

    // cout << "line 237" << endl;

    for (int i = 0; i < threads_num; ++i)
    {
        int start = i * file_per_thread;
        int end = (i + 1) * file_per_thread;
        works.emplace_back(searchRefs, miQuery, query_mi_cnts, q_cnt, reFnames, recSize, start, end,
                           res, res_r);
    }
    for (auto &t : works)
    {
        t.join();
    }

    f_align_res *final_res = new f_align_res[q_cnt];
    f_align_res tmp_res;

    int hit_cnt = 0;
    for (int q_id = 0; q_id < q_cnt; q_id++)
    {
        int max_hits = 0;
        if (!res[q_id].flag && !res_r[q_id].flag)
        {
            continue;
        }
        if (res[q_id].flag)
        {
            for (int k = 0; k < 48; k++)
            {
                auto i = res[q_id].align_res_array[k];
                if (!i.valid)
                    continue;
                if (i.hits > max_hits)
                {
                    tmp_res.f_id = k;
                    tmp_res.hits = i.hits;
                    tmp_res.second_hits = i.second_hits;
                    tmp_res.q_id = q_id;
                    tmp_res.t1_id = i.t1_id;
                    tmp_res.t2_id = i.t2_id;
                    tmp_res.strand = true;
                    max_hits = i.hits;
                }
            }
        }

        if (res_r[q_id].flag)
        {

            for (int k = 0; k < 48; k++)
            {
                auto i = res_r[q_id].align_res_array[k];
                if (!i.valid)
                    continue;
                if (i.hits > max_hits)
                {

                    tmp_res.f_id = k;
                    tmp_res.hits = i.hits;
                    tmp_res.second_hits = i.second_hits;
                    tmp_res.q_id = q_id;
                    tmp_res.t1_id = i.t1_id;
                    tmp_res.t2_id = i.t2_id;
                    tmp_res.strand = false;
                    max_hits = i.hits;
                }
            }
        }

        if (tmp_res.t1_id != 0 && tmp_res.t2_id != 0)
        {
            final_res[hit_cnt++] = tmp_res;
        }
    }

    for (int i = 0; i < hit_cnt; i++)
    {
        tmp_res = final_res[i];
        if (tmp_res.t2_id == 0)
            continue;
        std::string name1 = clstrNames[tmp_res.f_id][tmp_res.t1_id - 1];
        std::string name2 = clstrNames[tmp_res.f_id][tmp_res.t2_id - 1];

        if (name1.find('-') != std::string::npos || name2.find('-') != std::string::npos)
            continue;

        if (tmp_res.strand)
            outputTsv << query_names[tmp_res.q_id] << "\t" << formatted_Ropt[tmp_res.f_id] << "\tF\t" << name1 << "\t" << name2 << "\t" << tmp_res.hits << "\t" << tmp_res.second_hits << endl;
        else
            outputTsv << query_names[tmp_res.q_id] << "\t" << formatted_Ropt[tmp_res.f_id] << "\tR\t" << name1 << "\t" << name2 << "\t" << tmp_res.hits << "\t" << tmp_res.second_hits << endl;

        if (true)
        {
            outputReads << '>' << query_names[tmp_res.q_id] << "," << name1 << "," << name2 << "\n";
            int len_seq = sizeQuery[tmp_res.q_id];
            uint8_t *seq = seqQuery[tmp_res.q_id];
            for (int i = 0; i < len_seq; i++)
            {
                int idx = i / 4;
                int shift = 6 - 2 * (i % 4);
                uint8_t _4mer = seq[idx];
                uint8_t baseCode = (_4mer >> shift) & 0b11;
                char baseChar = hashToChar[baseCode];
                outputReads << baseChar;
            }
            outputReads << endl;
        }
    }

    delete[] final_res;

    for (int i = 0; i < q_cnt; i++)
    {
        delete[] miQuery[i];
        delete[] query_names[i];
        delete[] seqQuery[i];
    }
    delete[] miQuery;
    delete[] query_names;
    delete[] query_mi_cnts;
    delete[] seqQuery;
    delete[] sizeQuery;
    delete[] res;
    delete[] res_r;

    return curPosition;
}
