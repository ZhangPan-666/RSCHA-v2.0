#pragma once
#include <mpi.h>
#include <algorithm>
#include <climits>
#include <cstddef>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

namespace MPIParallel {

inline void check(int error, const char *operation) {
  if (error == MPI_SUCCESS) return;
  char message[MPI_MAX_ERROR_STRING]{};
  int length = 0;
  MPI_Error_string(error, message, &length);
  throw std::runtime_error(std::string(operation) + ": " +
      std::string(message, static_cast<std::size_t>(length)));
}

class Environment {
public:
  Environment(int &argc, char **&argv) {
    int initialized = 0;
    check(MPI_Initialized(&initialized), "MPI_Initialized");
    if (!initialized) {
      check(MPI_Init(&argc, &argv), "MPI_Init");
      owns_mpi_ = true;
    }
    check(MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN), "MPI_Comm_set_errhandler");
    check(MPI_Comm_rank(MPI_COMM_WORLD, &rank_), "MPI_Comm_rank");
    check(MPI_Comm_size(MPI_COMM_WORLD, &size_), "MPI_Comm_size");
  }
  Environment(const Environment &) = delete;
  Environment &operator=(const Environment &) = delete;
  ~Environment() {
    int finalized = 0;
    MPI_Finalized(&finalized);
    if (owns_mpi_ && !finalized) MPI_Finalize();
  }
  int rank() const noexcept { return rank_; }
  int size() const noexcept { return size_; }
  bool is_root() const noexcept { return rank_ == 0; }
private:
  bool owns_mpi_ = false;
  int rank_ = 0, size_ = 1;
};

inline bool active() {
  int initialized = 0, finalized = 0;
  MPI_Initialized(&initialized);
  if (initialized) MPI_Finalized(&finalized);
  return initialized && !finalized;
}
inline int rank() {
  int value = 0;
  if (active()) check(MPI_Comm_rank(MPI_COMM_WORLD, &value), "MPI_Comm_rank");
  return value;
}
inline int size() {
  int value = 1;
  if (active()) check(MPI_Comm_size(MPI_COMM_WORLD, &value), "MPI_Comm_size");
  return value;
}
inline double time() { return active() ? MPI_Wtime() : 0.0; }
inline void barrier() {
  if (active() && size() > 1) check(MPI_Barrier(MPI_COMM_WORLD), "MPI_Barrier");
}
inline void phase_begin(const char *label) {
  if (rank() == 0) {
    std::printf("[MPI stage] BEGIN %s; time=%.6f\n", label, time());
    std::fflush(stdout);
  }
}
inline void phase_complete(const char *label, double started) {
  barrier();
  if (rank() == 0) {
    std::printf("[MPI stage] COMPLETE %s; all ranks arrived; elapsed=%.6f s\n",
        label, time()-started);
    std::fflush(stdout);
  }
}

inline std::vector<int> owned_indices(int total, int owner, int processes) {
  if (total < 0 || processes <= 0 || owner < 0 || owner >= processes)
    throw std::invalid_argument("Invalid MPI data partition.");
  const int base = total / processes, extra = total % processes;
  const int count = base + (owner < extra ? 1 : 0);
  const int begin = owner * base + std::min(owner, extra);
  std::vector<int> result(static_cast<std::size_t>(count));
  for (int i = 0; i < count; ++i) result[i] = begin+i;
  return result;
}

inline void broadcast(std::vector<double> &values) {
  if (!active() || size() == 1) return;
  for (std::size_t offset = 0; offset < values.size();) {
    const int count = static_cast<int>(std::min(values.size()-offset, static_cast<std::size_t>(INT_MAX)));
    check(MPI_Bcast(values.data()+offset, count, MPI_DOUBLE, 0, MPI_COMM_WORLD), "MPI_Bcast vector");
    offset += static_cast<std::size_t>(count);
  }
}

// Retain the public helper names; implement the Fortran gather/sum/broadcast architecture.
inline void allreduce_sum(double &value) {
  if (!active() || size() == 1) return;
  const int processes = size(), owner = rank();
  std::vector<int> counts(processes, 1), offsets(processes);
  for (int i = 0; i < processes; ++i) offsets[i] = i;
  std::vector<double> gathered(owner == 0 ? processes : 0);
  barrier();
  check(MPI_Gatherv(&value, 1, MPI_DOUBLE, owner == 0 ? gathered.data() : nullptr,
      counts.data(), offsets.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD), "MPI_Gatherv scalar");
  if (owner == 0) {
    value = 0.0;
    for (int i = 0; i < processes; ++i) value += gathered[i];
  }
  check(MPI_Bcast(&value, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD), "MPI_Bcast scalar");
}

inline void allreduce_sum(std::vector<double> &values) {
  if (!active() || size() == 1) return;
  const int processes = size(), owner = rank();
  const unsigned long long length = static_cast<unsigned long long>(values.size());
  std::vector<unsigned long long> lengths(owner == 0 ? processes : 0);
  check(MPI_Gather(&length, 1, MPI_UNSIGNED_LONG_LONG, owner == 0 ? lengths.data() : nullptr,
      1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD), "MPI_Gather vector lengths");
  int valid = 1;
  if (owner == 0)
    for (int i = 0; i < processes; ++i) if (lengths[i] != length) valid = 0;
  check(MPI_Bcast(&valid, 1, MPI_INT, 0, MPI_COMM_WORLD), "MPI_Bcast vector validity");
  if (!valid) throw std::invalid_argument("MPI sum vector lengths differ between ranks.");
  if (values.empty()) return;
#ifndef RSCHA_MPI_GATHER_BUFFER_VALUES
#define RSCHA_MPI_GATHER_BUFFER_VALUES (8u * 1024u * 1024u)
#endif
  // Limit root scratch storage to 64 MiB and keep all MPI counts/displacements in int range.
  const std::size_t block_limit = std::min(
      static_cast<std::size_t>(RSCHA_MPI_GATHER_BUFFER_VALUES) / processes,
      static_cast<std::size_t>(INT_MAX) / processes);
  if (block_limit == 0) throw std::runtime_error("MPI process count exceeds gather buffer capacity.");
  const std::size_t block_size = std::min(values.size(), block_limit);
  std::vector<double> gathered(owner == 0 ? block_size * processes : 0);
  std::vector<int> counts(processes), offsets(processes);
  barrier();
  for (std::size_t offset = 0; offset < values.size();) {
    const int count = static_cast<int>(std::min(values.size()-offset, block_size));
    for (int i = 0; i < processes; ++i) { counts[i] = count; offsets[i] = i*count; }
    check(MPI_Gatherv(values.data()+offset, count, MPI_DOUBLE, owner == 0 ? gathered.data() : nullptr,
        counts.data(), offsets.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD), "MPI_Gatherv sum block");
    if (owner == 0) {
      for (int j = 0; j < count; ++j) {
        double sum = 0.0;
        for (int i = 0; i < processes; ++i) sum += gathered[static_cast<std::size_t>(i)*count+j];
        values[offset+j] = sum;
      }
    }
    check(MPI_Bcast(values.data()+offset, count, MPI_DOUBLE, 0, MPI_COMM_WORLD), "MPI_Bcast sum block");
    offset += static_cast<std::size_t>(count);
  }
}

inline void allreduce_sum(std::vector<std::vector<double>> &matrix) {
  if (!active() || size() == 1 || matrix.empty()) return;
  const std::size_t columns = matrix.front().size();
  for (const auto &row : matrix)
    if (row.size() != columns) throw std::runtime_error("MPI matrix sum requires equal row lengths.");
  if (columns == 0) return;
  // Preserve the matrix helper interface; the inversion design matrix does not use it.
  const std::size_t rows_per_block = std::max<std::size_t>(1,
      std::min(matrix.size(), static_cast<std::size_t>(8u*1024u*1024u)/columns));
  std::vector<double> flat(rows_per_block*columns);
  for (std::size_t first = 0; first < matrix.size(); first += rows_per_block) {
    const std::size_t rows = std::min(rows_per_block, matrix.size()-first);
    flat.resize(rows*columns);
    for (std::size_t row = 0; row < rows; ++row)
      std::copy(matrix[first+row].begin(), matrix[first+row].end(), flat.begin()+row*columns);
    allreduce_sum(flat);
    for (std::size_t row = 0; row < rows; ++row)
      std::copy_n(flat.begin()+row*columns, columns, matrix[first+row].begin());
  }
}

// Gather only owned prediction entries, rather than summing a sparse full-data vector on every rank.
class IndexedGather {
public:
  IndexedGather(int total, const std::vector<int> &rows)
      : total_(total), local_(static_cast<int>(rows.size())), rows_(rows) {
    const int processes = size(), owner = rank();
    counts_.assign(processes, 0); offsets_.assign(processes, 0);
    int valid = 1;
    if (processes > 1) {
      check(MPI_Gather(&local_, 1, MPI_INT, counts_.data(), 1, MPI_INT, 0, MPI_COMM_WORLD), "MPI_Gather row counts");
    } else counts_[0] = local_;
    if (owner == 0) {
      long long count = 0;
      for (int i = 0; i < processes; ++i) {
        if (count > INT_MAX || counts_[i] < 0) valid = 0;
        offsets_[i] = count <= INT_MAX ? static_cast<int>(count) : 0;
        count += counts_[i];
      }
      if (count != total_ || total_ <= 0) valid = 0;
    }
    if (processes > 1) check(MPI_Bcast(&valid, 1, MPI_INT, 0, MPI_COMM_WORLD), "MPI_Bcast row count validity");
    if (!valid) throw std::invalid_argument("MPI data partition does not cover the global data count.");
    if (owner == 0) { order_.resize(total_); buffer_.resize(total_); }
    int dummy = 0;
    if (processes > 1)
      check(MPI_Gatherv(rows.empty() ? &dummy : rows.data(), local_, MPI_INT,
          owner == 0 ? order_.data() : nullptr, counts_.data(), offsets_.data(), MPI_INT,
          0, MPI_COMM_WORLD), "MPI_Gatherv global row indices");
    else order_ = rows;
    if (owner == 0) {
      std::vector<unsigned char> seen(total_, 0);
      for (int row : order_) {
        if (row < 0 || row >= total_ || seen[row]) { valid = 0; break; }
        seen[row] = 1;
      }
    }
    if (processes > 1) check(MPI_Bcast(&valid, 1, MPI_INT, 0, MPI_COMM_WORLD), "MPI_Bcast row index validity");
    if (!valid) throw std::invalid_argument("MPI data rows are missing, duplicated, or outside the global range.");
  }
  void collect(const std::vector<double> &local, std::vector<double> &global) {
    if (local.size() != static_cast<std::size_t>(local_) || global.size() != static_cast<std::size_t>(total_))
      throw std::invalid_argument("Inconsistent MPI prediction gather dimensions.");
    const int owner = rank();
    double dummy = 0.0;
    barrier();
    if (size() > 1)
      check(MPI_Gatherv(local.empty() ? &dummy : local.data(), local_, MPI_DOUBLE,
          owner == 0 ? buffer_.data() : nullptr, counts_.data(), offsets_.data(), MPI_DOUBLE,
          0, MPI_COMM_WORLD), "MPI_Gatherv predictions");
    else buffer_ = local;
    if (owner == 0)
      for (int i = 0; i < total_; ++i) global[order_[i]] = buffer_[i];
    broadcast(global);
  }
private:
  int total_, local_;
  std::vector<int> rows_, counts_, offsets_, order_;
  std::vector<double> buffer_;
};

inline void communication_check() {
  phase_begin("communication preflight");
  const double started = time(), local = double(rank())+1.0;
  const double expected = double(size())*(double(size())+1.0)/2.0;
  double value = local;
  allreduce_sum(value);
  if (value != expected) throw std::runtime_error("MPI scalar preflight failed.");
  for (int count : {72010, 144124}) {
    std::vector<double> values(count, local);
    allreduce_sum(values);
    for (double v : values) if (v != expected) throw std::runtime_error("MPI vector preflight failed.");
  }
  std::vector<int> rows{2*rank()+1, 2*rank()};
  IndexedGather gather(2*size(), rows);
  std::vector<double> local_values{double(rows[0]), double(rows[1])}, global(2*size());
  gather.collect(local_values, global);
  for (int i = 0; i < 2*size(); ++i)
    if (global[i] != i) throw std::runtime_error("MPI indexed gather preflight failed.");
  phase_complete("communication preflight PASSED", started);
}
}
