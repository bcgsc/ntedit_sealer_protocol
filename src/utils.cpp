#include "utils.hpp"

#include "btllib/status.hpp"
#include "btllib/nthash.hpp"
#include "btllib/bloom_filter.hpp"
#include "btllib/counting_bloom_filter.hpp"

#include <algorithm>
#include <cstdlib>
#include <random>
#include <fstream>
#include <thread>
#include <chrono>
#include <memory>
#include <vector>
#include <utility>
#include <string>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

std::vector<size_t> get_random_indices(const size_t total_size, const size_t count) {
  btllib::check_error(count > total_size, "get_random_indices: count cannot be larger than total_size.");

  static std::random_device dev;
  static std::mt19937 rng(dev());

  std::vector<size_t> all_indices(total_size);
  std::iota(all_indices.begin(), all_indices.end(), 0);
  std::shuffle(all_indices.begin(), all_indices.end(), rng);
  decltype(all_indices) indices(all_indices.begin(), all_indices.begin() + count);

  return indices;
}

void load_index(Index& index, const std::string& filepath) {
  btllib::log_info(std::string("load_index: Loading index from ") + filepath + "... ");
  std::ifstream ifs(filepath);
  std::string token, id;
  unsigned long id_start, id_end, seq_end;
  unsigned long i = 0;
  while (ifs >> token) {
    switch (i % 4) {
      case 0: id = std::move(token); break;
      case 1: id_start = std::stoul(token); break;
      case 2: id_end = std::stoul(token); break;
      case 3: {
        seq_end = std::stoul(token);
        index.emplace(std::piecewise_construct, std::make_tuple(id), std::make_tuple(id_end + 1, seq_end - id_end - 1));
        break;
      }
      default: {
        btllib::log_error("load_index: Invalid switch branch.");
        std::exit(EXIT_FAILURE);
      }
    }
    i++;
  }
  btllib::log_info("load_index: Done!");
}

void load_reads_mapping(ReadsMapping& reads_mapping, const std::string& filepath, const unsigned mx_threshold) {
  btllib::log_info(std::string("load_reads_mapping: Loading contig reads from ") + filepath + "... ");
  std::ifstream ifs(filepath);
  std::string token, read_id, contig_id;
  unsigned long i = 0;
  while (ifs >> token) {
    switch (i % 3) {
      case 0: read_id = std::move(token); break;
      case 1: contig_id = std::move(token); break;
      case 2: {
        auto it = reads_mapping.find(contig_id);
        if (it == reads_mapping.end()) {
          const auto emplacement = reads_mapping.emplace(contig_id, std::vector<SeqId>());
          it = emplacement.first;
        }
        const auto minimizers = std::stoul(token);
        if (minimizers >= mx_threshold) {
          it->second.push_back(read_id);
        }
        break;
      }
      default: {
        btllib::log_error("load_reads_mapping: Invalid switch branch.");
        std::exit(EXIT_FAILURE);
      }
    }
    i++;
  }
  btllib::log_info("load_reads_mapping: Done!");
}

void fill_bfs(const char* seq,
              const size_t seq_len,
              const unsigned hash_num,
              const std::vector<unsigned>& ks,
              const unsigned kmer_threshold,
              std::vector<std::unique_ptr<btllib::KmerCountingBloomFilter8>>& cbfs,
              std::vector<std::unique_ptr<btllib::KmerBloomFilter>>& bfs)
{
  for (size_t i = 0; i < ks.size(); i++)
  {
    const auto& cbf = cbfs[i];
    const auto& bf = bfs[i];
    btllib::NtHash nthash(seq, seq_len, hash_num, ks[i]);
    while (nthash.roll()) {
      if (cbf->insert_thresh_contains(nthash.hashes(), kmer_threshold) >= kmer_threshold) {
        bf->insert(nthash.hashes());
      }
    }
  }
}

void start_watchdog() {
  (new std::thread([] () {
    while (getppid() != 1) {
      std::this_thread::sleep_for(std::chrono::seconds(1));
    }
    std::exit(EXIT_FAILURE);
  }))->detach();
}