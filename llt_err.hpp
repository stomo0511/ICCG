#pragma once

#include <string>

#include "crs_io.hpp"
#include "ic0.hpp"

struct LLTErrorResult {
    double norm_A_F = 0.0;
    double norm_LLt_F = 0.0;
    double inner_A_LLt = 0.0;
    double norm_error_F = 0.0;
    double epsilon_F = 0.0;
};

LLTErrorResult EvaluateLLTError(const CRS& A_lower, const IC0& M);

double DenseCheckLLTEpsilon(const CRS& A_lower, const IC0& M);

CRS LowerCRSWithDiag(const CRS& A);

long long SymmetricNNZFromLower(const CRS& A_lower);

std::string MatrixNameFromPath(const std::string& path);
