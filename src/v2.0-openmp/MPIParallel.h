#pragma once

#include <mpi.h>

#include <algorithm>
#include <climits>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace MPIParallel {

class Environment {
public:
  Environment(int &argc, char **&argv) {
    int initialized = 0;
    MPI_Initialized(&initialized);
    if (!initialized) {
      check(MPI_Init(&argc, &argv), "MPI_Init");
      owns_mpi_ = true;
    }
    check(MPI_Comm_rank(MPI_COMM_WORLD, &rank_), "MPI_Comm_rank");
    check(MPI_Comm_size(MPI_COMM_WORLD, &size_), "MPI_Comm_size");
  }

  Environment(const Environment &) = delete;
  Environment &operator=(const Environment &) = delete;

  ~Environment() {
    int finalized = 0;
    MPI_Finalized(&finalized);
    if (owns_mpi_ && !finalized)
      MPI_Finalize();
  }

  int rank() const noexcept { return rank_; }
  int size() const noexcept { return size_; }
  bool is_root() const noexcept { return rank_ == 0; }

private:
  static void check(int error, const char *operation) {
    if (error == MPI_SUCCESS)
      return;
    char message[MPI_MAX_ERROR_STRING]{};
    int length = 0;
    MPI_Error_string(error, message, &length);
    throw std::runtime_error(
        std::string(operation) + ": " +
        std::string(message, static_cast<std::size_t>(length)));
  }

  bool owns_mpi_ = false;
  int rank_ = 0;
  int size_ = 1;
};

inline bool active() {
  int initialized = 0;
  int finalized = 0;
  MPI_Initialized(&initialized);
  if (initialized)
    MPI_Finalized(&finalized);
  return initialized && !finalized;
}

inline int rank() {
  if (!active())
    return 0;
  int value = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &value);
  return value;
}

inline int size() {
  if (!active())
    return 1;
  int value = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &value);
  return value;
}

inline void allreduce_sum(double &value) {
  if (!active() || size() == 1)
    return;
  double result = 0.0;
  MPI_Allreduce(&value, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  value = result;
}

inline void allreduce_sum(std::vector<double> &values) {
  if (!active() || size() == 1 || values.empty())
    return;
  std::size_t offset = 0;
  while (offset < values.size()) {
    const int count = static_cast<int>(std::min<std::size_t>(
        values.size() - offset, static_cast<std::size_t>(INT_MAX)));
    MPI_Allreduce(MPI_IN_PLACE, values.data() + offset, count, MPI_DOUBLE,
                  MPI_SUM, MPI_COMM_WORLD);
    offset += static_cast<std::size_t>(count);
  }
}

inline void allreduce_sum(std::vector<std::vector<double>> &matrix) {
  if (!active() || size() == 1 || matrix.empty())
    return;
  const std::size_t columns = matrix.front().size();
  for (const auto &row : matrix) {
    if (row.size() != columns)
      throw std::runtime_error(
          "MPI matrix reduction requires equal row lengths");
  }
  if (columns == 0)
    return;
  const std::size_t target_values = 8u * 1024u * 1024u;
  const std::size_t rows_per_block = std::max<std::size_t>(
      1, std::min<std::size_t>(matrix.size(), target_values / columns));
  std::vector<double> flat(rows_per_block * columns, 0.0);
  for (std::size_t first = 0; first < matrix.size(); first += rows_per_block) {
    const std::size_t rows = std::min(rows_per_block, matrix.size() - first);
    flat.resize(rows * columns);
    for (std::size_t row = 0; row < rows; ++row)
      std::copy(matrix[first + row].begin(), matrix[first + row].end(),
                flat.begin() + static_cast<std::ptrdiff_t>(row * columns));
    allreduce_sum(flat);
    for (std::size_t row = 0; row < rows; ++row)
      std::copy_n(flat.begin() + static_cast<std::ptrdiff_t>(row * columns),
                  columns, matrix[first + row].begin());
  }
}

}
