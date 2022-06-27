#include <fstream>
#include <iostream>
#include <sstream>

#include <map>
#include <string>
#include <vector>

size_t ham_dist(const std::string& s1, const std::string& s2) {
    if (s1.size() != s2.size()) {
        throw std::logic_error("Hamming distance assumes that string lengths are equal!");
    }

    size_t dist = 0;
    for (size_t i = 0; i < s1.size(); ++i) {
        if (s1[i] != s2[i]) {
            ++dist;
        }
    }

    return dist;
}

std::map<std::string, std::vector<size_t>> build_kmers_map(const std::string& genome, size_t kmer_len) {
    std::map<std::string, std::vector<size_t>> kmers;
    size_t pos = 0;
    for (size_t i = 0; i < genome.size() - kmer_len + 1; ++i) {
        std::string kmer = genome.substr(i, kmer_len);
        if (kmers.find(kmer) == kmers.end()) {
            kmers.insert({kmer, {i}});
        } else {
            kmers[kmer].push_back(i);
        }
    }

    return kmers;
}

std::map<std::string, int> simple_reads_mapping(const std::string& genome, std::vector<std::string> reads, size_t kmer_len) {
    std::map<std::string, int> mapped_reads;
    std::map<std::string, std::vector<size_t>> kmers_map = build_kmers_map(genome, kmer_len);
    size_t count = 0;
    for (const auto& read : reads) {
        if (count % 10000 == 0) {
            std::cout << count << " reads were processed" << std::endl;
        }

        std::string seed = read.substr(0, kmer_len);
        bool found = false;
        auto seed_it = kmers_map.find(seed);
        std::vector<size_t> poss = {};

        if (seed_it != kmers_map.end()) {
            poss = seed_it->second;
        }

        for (size_t pos : poss) {
            if (pos + read.size() < genome.size()) {
                size_t d = ham_dist(read, genome.substr(pos, read.size()));
                if (d == 0) {
                    mapped_reads.insert({read, static_cast<int>(pos)});
                    found = true;
                }
            }
        }

        if (!found) {
            mapped_reads.insert({ read, -1 }); // no perfect match for read in genome
        }
        ++count;
    }

    return mapped_reads;
}

int main() {
    std::ifstream genome_file;
    std::ifstream reads_file;

    //std::string genome_file_name = "test_example/genome.fasta";
    std::string genome_file_name = "../data/genome.fasta";
    genome_file.open(genome_file_name );
    if (!genome_file) {
        throw std::runtime_error("Unable to open file " + genome_file_name);
    }

    //std::string reads_file_name = "test_example/reads.txt";
    std::string reads_file_name = "../data/reads.txt";
    reads_file.open(reads_file_name );
    if (!reads_file) {
        throw std::runtime_error("Unable to open file " + reads_file_name);
    }

    std::string line;
    std::string genome = "";
    std::getline(genome_file, line); // skip header
    while(std::getline (genome_file, line)) {
      genome += line;
    }
    genome_file.close();

    std::vector<std::string> reads;
    while (std::getline(reads_file, line)) {
      reads.push_back(line);
    }
    reads_file.close();
    auto mapped_reads = simple_reads_mapping(genome, reads, 3);


    std::ofstream res("res.txt");
    for (const auto& item : mapped_reads) {
        res << item.first << " " << item.second << std::endl;
    }

    res.close();
}