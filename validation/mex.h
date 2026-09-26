// Minimal validation-only MEX array shim. This is not the MATLAB ABI.
#pragma once
#include <cstddef>
#include <vector>
#include <stdexcept>
struct mxArray
{
    std::size_t rows, cols;
    std::vector<double> values;
    mxArray(std::size_t m, std::size_t n) : rows(m), cols(n), values(m * n) {}
};
constexpr int mxREAL = 0;
inline double *mxGetPr(const mxArray *a)
{
    return const_cast<double *>(a->values.data());
}
inline std::size_t mxGetM(const mxArray *a)
{
    return a->rows;
}
inline std::size_t mxGetN(const mxArray *a)
{
    return a->cols;
}
inline mxArray *mxCreateDoubleMatrix(std::size_t m, std::size_t n, int)
{
    return new mxArray(m, n);
}

inline bool mxIsDouble(const mxArray *) { return true; }
inline bool mxIsComplex(const mxArray *) { return false; }
inline bool mxIsSparse(const mxArray *) { return false; }
inline std::size_t mxGetNumberOfElements(const mxArray *a) { return a->values.size(); }
[[noreturn]] inline void mexErrMsgIdAndTxt(const char *, const char *message)
{
    throw std::runtime_error(message);
}
