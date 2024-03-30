#include "constants.hpp"
#include <fstream>
#include <cstring>


const int file_end = 48;
const char *reFnames[] = {
    "./hg38REF/chr21.reverse.cDNA.fa",
    "./hg38REF/chr1.forward.cDNA.fa",
    "./hg38REF/chr21.forward.cDNA.fa",
    "./hg38REF/chr1.reverse.cDNA.fa",
    "./hg38REF/chr13.reverse.cDNA.fa",
    "./hg38REF/chr2.forward.cDNA.fa",
    "./hg38REF/chr18.reverse.cDNA.fa",
    "./hg38REF/chr2.reverse.cDNA.fa",
    "./hg38REF/chr18.forward.cDNA.fa",
    "./hg38REF/chr3.forward.cDNA.fa",
    "./hg38REF/chr13.forward.cDNA.fa",
    "./hg38REF/chr3.reverse.cDNA.fa",
    "./hg38REF/chr22.reverse.cDNA.fa",
    "./hg38REF/chr17.reverse.cDNA.fa",
    "./hg38REF/chr22.forward.cDNA.fa",
    "./hg38REF/chr11.reverse.cDNA.fa",
    "./hg38REF/chr20.reverse.cDNA.fa",
    "./hg38REF/chr11.forward.cDNA.fa",
    "./hg38REF/chr20.forward.cDNA.fa",
    "./hg38REF/chr19.forward.cDNA.fa",
    "./hg38REF/chr9.reverse.cDNA.fa",
    "./hg38REF/chr19.reverse.cDNA.fa",
    "./hg38REF/chr14.reverse.cDNA.fa",
    "./hg38REF/chr12.forward.cDNA.fa",
    "./hg38REF/chr9.forward.cDNA.fa",
    "./hg38REF/chr17.forward.cDNA.fa",
    "./hg38REF/chr15.reverse.cDNA.fa",
    "./hg38REF/chr12.reverse.cDNA.fa",
    "./hg38REF/chr15.forward.cDNA.fa",
    "./hg38REF/chr5.forward.cDNA.fa",
    "./hg38REF/chr10.reverse.cDNA.fa",
    "./hg38REF/chr16.forward.cDNA.fa",
    "./hg38REF/chr10.forward.cDNA.fa",
    "./hg38REF/chr6.reverse.cDNA.fa",
    "./hg38REF/chr4.reverse.cDNA.fa",
    "./hg38REF/chr7.forward.cDNA.fa",
    "./hg38REF/chr8.forward.cDNA.fa",
    "./hg38REF/chr7.reverse.cDNA.fa",
    "./hg38REF/chr8.reverse.cDNA.fa",
    "./hg38REF/chr6.forward.cDNA.fa",
    "./hg38REF/chr14.forward.cDNA.fa",
    "./hg38REF/chr5.reverse.cDNA.fa",
    "./hg38REF/chr16.reverse.cDNA.fa",
    "./hg38REF/chr4.forward.cDNA.fa",
    "./hg38REF/chrX.forward.cDNA.fa",
    "./hg38REF/chrX.reverse.cDNA.fa",
    "./hg38REF/chrY.forward.cDNA.fa",
    "./hg38REF/chrY.reverse.cDNA.fa",
};

const int recSize[] = {
    1236, 10360, 1507, 9524, 1899, 8224, 1944, 8132, 1957, 7605, 2199, 7123, 2227,
    6946, 2492, 6740, 2513, 6709, 2763, 6651, 3652, 6412, 3657, 6204, 3794, 6114,
    3846, 5935, 3918, 5534, 3924, 5363, 4206, 5181, 4255, 5178, 4320, 5144, 4426,
    5082, 4515, 4751, 4533, 4603, 3710, 3595, 476, 390};

std::vector<std::vector<std::string>> clstrNames(file_end);

int clstrSizes[] = {288, 1954, 276, 1919, 549, 1491, 338, 1429, 338, 1246, 527, 1251, 433, 1026, 495, 1195, 505, 1198, 531, 1101,
                    808, 994, 782, 945, 794, 931, 656, 963, 707, 951, 811, 767, 850, 1030, 828, 1059, 713, 925, 795, 1054, 826,
                    922, 709, 883, 945, 939, 241, 231};

char **formatted_Ropt = new char *[file_end];


    
int kmer_l;
int kmer_2l;
int interval;
int window_s;



void initializeConstants()
{

    std::ifstream nameFile("./hg38REF/geneName.txt");

    for (int i = 0; i < file_end; ++i)
    {

        std::string geneName;

        int c_size = clstrSizes[i];

        clstrNames[i].resize(c_size);

        for (int j = 0; j < c_size && std::getline(nameFile, geneName); ++j)
        {
            clstrNames[i][j] = geneName;
        }
    }

    nameFile.close();

    for (int i = 0; i < file_end; ++i)
    {
        const char *geneName = reFnames[i];

        const char *chr_s = geneName + 10;

        const char *chr_e = strchr(chr_s, '.');

        size_t chr_l = chr_e - chr_s;

        char *chromosome = new char[chr_l + 1];

        strncpy(chromosome, chr_s, chr_l);
        
        chromosome[chr_l] = '\0';

        const char *ort = (strstr(geneName, "reverse") != nullptr) ? "-" : "+";
        
        size_t formattedStringLength = strlen(chromosome) + 1 + strlen(ort) + 1;
        
        formatted_Ropt[i] = new char[formattedStringLength];
        
        sprintf(formatted_Ropt[i], "%s\t%s", chromosome, ort);
        
        delete[] chromosome;
    }
}
