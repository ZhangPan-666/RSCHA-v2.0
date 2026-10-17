#pragma once
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

// Presentation only: do not change numerical settings or MPI collectives here.
namespace RSCHAConsole {
inline void banner(const char* title, const char* backend) {
  std::cout << '\n' << std::string(72, '=') << '\n'
            << "  " << title << " | " << backend << '\n'
            << std::string(72, '=') << '\n';
}
inline void section(int number, const char* title) {
  std::ostringstream label;
  label << '[' << std::setw(2) << std::setfill('0') << number << "] " << title;
  std::cout << '\n' << label.str() << '\n' << std::string(72, '-') << '\n' << std::flush;
}
template<class T> inline void row(const char* label, const T& value, const char* unit = "") {
  std::ostringstream line;
  line << "  " << std::left << std::setw(32) << label << ": "
       << std::setprecision(8) << value;
  if (*unit) line << ' ' << unit;
  std::cout << line.str() << '\n';
}
inline void model(double lon, double lat, double theta, int ki, int ke,
                  int p, int m, double lower, double upper) {
  row("Cap centre longitude", lon, "deg"); row("Cap centre latitude", lat, "deg");
  row("Cap half-angle", theta, "deg"); row("Internal / external index", std::to_string(ki)+" / "+std::to_string(ke));
  row("Mehler / constant index", std::to_string(p)+" / "+std::to_string(m));
  row("Lower / upper radius", std::to_string(lower/1000)+" / "+std::to_string(upper/1000), "km");
}
inline void kernel(const char* label) { std::cout << "\n  Kernel: " << label << '\n' << std::flush; }
inline void done(const char* label) { std::cout << "  [OK] " << label << '\n' << std::flush; }
inline void complete(const char* title) {
  std::cout << '\n' << std::string(72, '=') << '\n'
            << "  [OK] " << title << '\n' << std::string(72, '=') << '\n' << std::flush;
}
}
