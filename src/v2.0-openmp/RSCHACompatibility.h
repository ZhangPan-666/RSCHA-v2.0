#pragma once

// Project-local replacements for library facilities missing on older clusters.
// No symbols are injected into namespace std; numerical comparisons stay unchanged.
#include <cerrno>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>
#include <sys/stat.h>
#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#include <direct.h>
#else
#include <unistd.h>
#endif

namespace RSCHACompat {
template<class T>
constexpr const T& clamp(const T& value, const T& lower, const T& upper) {
    return value < lower ? lower : upper < value ? upper : value;
}

class path {
public:
#ifdef _WIN32
    typedef wchar_t value_type;
    typedef std::wstring string_type;
    path(const wchar_t* value) : value_(value ? value : L"") {}
    path(const std::wstring& value) : value_(value) {}
    path(const char* value) : value_(widen(value ? std::string(value) : std::string())) {}
    path(const std::string& value) : value_(widen(value)) {}
#else
    typedef char value_type;
    typedef std::string string_type;
    path(const char* value) : value_(value ? value : "") {}
    path(const std::string& value) : value_(value) {}
#endif
    path() {}
    const value_type* c_str() const { return value_.c_str(); }
    const string_type& native() const { return value_; }
    bool empty() const { return value_.empty(); }
    std::string string() const {
#ifdef _WIN32
        if (value_.empty()) return std::string();
        const int count = WideCharToMultiByte(CP_ACP, 0, value_.c_str(), -1, nullptr, 0, nullptr, nullptr);
        if (!count) throw std::runtime_error("Cannot encode filesystem path");
        std::vector<char> buffer(static_cast<std::size_t>(count));
        if (!WideCharToMultiByte(CP_ACP, 0, value_.c_str(), -1, buffer.data(), count, nullptr, nullptr))
            throw std::runtime_error("Cannot encode filesystem path");
        return std::string(buffer.data());
#else
        return value_;
#endif
    }
    bool is_absolute() const {
        if (value_.empty()) return false;
#ifdef _WIN32
        return separator(value_[0]) || (value_.size() >= 3 && value_[1] == L':' && separator(value_[2]));
#else
        return value_[0] == '/';
#endif
    }
    path parent_path() const {
        if (value_.empty()) return path();
        std::size_t end = value_.size();
        while (end > 1 && separator(value_[end-1])) {
#ifdef _WIN32
            if (end == 3 && value_[1] == L':') break;
#endif
            --end;
        }
        std::size_t slash = end;
        while (slash > 0 && !separator(value_[slash-1])) --slash;
        if (!slash) return path();
        if (slash == 1) return from_native(value_.substr(0,1));
#ifdef _WIN32
        if (slash == 3 && value_[1] == L':') return from_native(value_.substr(0,3));
#endif
        return from_native(value_.substr(0,slash-1));
    }
    friend path operator/(const path& left, const path& right) {
        if (right.is_absolute() || left.empty()) return right;
        string_type joined = left.value_;
        if (!separator(joined.back())) joined.push_back(preferred_separator());
        joined += right.value_;
        return from_native(joined);
    }
private:
    string_type value_;
    static value_type preferred_separator() {
#ifdef _WIN32
        return L'\\';
#else
        return '/';
#endif
    }
    static bool separator(value_type ch) {
#ifdef _WIN32
        return ch == L'/' || ch == L'\\';
#else
        return ch == '/';
#endif
    }
    static path from_native(const string_type& value) { path result; result.value_=value; return result; }
#ifdef _WIN32
    static std::wstring widen(const std::string& value) {
        if (value.empty()) return std::wstring();
        const int count=MultiByteToWideChar(CP_ACP,0,value.c_str(),-1,nullptr,0);
        if (!count) throw std::runtime_error("Cannot decode filesystem path");
        std::vector<wchar_t> buffer(static_cast<std::size_t>(count));
        if (!MultiByteToWideChar(CP_ACP,0,value.c_str(),-1,buffer.data(),count))
            throw std::runtime_error("Cannot decode filesystem path");
        return std::wstring(buffer.data());
    }
#endif
};

inline path current_path() {
    for (std::size_t size=256; size <= 1048576; size*=2) {
        std::vector<path::value_type> buffer(size);
        errno=0;
#ifdef _WIN32
        if (_wgetcwd(buffer.data(),static_cast<int>(size))) return path(buffer.data());
#else
        if (::getcwd(buffer.data(),size)) return path(buffer.data());
#endif
        if (errno != ERANGE) throw std::runtime_error(std::string("Cannot get current directory: ")+std::strerror(errno));
    }
    throw std::runtime_error("Current directory exceeds supported path length");
}

inline bool path_status(const path& value, bool* directory=nullptr) {
    if (value.empty()) return false;
#ifdef _WIN32
    struct _stat64 info;
    const int status=_wstat64(value.c_str(),&info);
#else
    struct stat info;
    const int status=::stat(value.c_str(),&info);
#endif
    if (status == 0) {
        if (directory) {
#ifdef _WIN32
            *directory=(info.st_mode & _S_IFMT) == _S_IFDIR;
#else
            *directory=S_ISDIR(info.st_mode);
#endif
        }
        return true;
    }
    if (errno == ENOENT || errno == ENOTDIR) return false;
    throw std::runtime_error(std::string("Cannot inspect filesystem path: ")+std::strerror(errno));
}
inline bool exists(const path& value) { return path_status(value); }
inline bool create_directories(const path& value) {
    if (value.empty()) throw std::runtime_error("Cannot create empty directory path");
    bool directory=false;
    if (path_status(value,&directory)) {
        if (!directory) throw std::runtime_error("Output directory path already names a file");
        return false;
    }
    const path parent=value.parent_path();
    if (!parent.empty() && parent.native()!=value.native()) create_directories(parent);
#ifdef _WIN32
    const int status=_wmkdir(value.c_str());
#else
    const int status=::mkdir(value.c_str(),0777);
#endif
    if (status == 0) return true;
    if (errno == EEXIST && path_status(value,&directory) && directory) return false;
    throw std::runtime_error(std::string("Cannot create output directory: ")+std::strerror(errno));
}
} // namespace RSCHACompat
