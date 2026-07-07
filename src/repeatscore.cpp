#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cstring>
using namespace Rcpp;

// Repeat score for a vector of sequences, matching the pure-R repeatScore()
// exactly. For each sequence it counts, among all length-`wsize` windows, the
// fraction that belong to a k-mer occurring more than `dup_cutoff` times, i.e.
//   sum(cnt[cnt > dup_cutoff]) / (nchar - wsize + 1).
// The duplicate structure is invariant to how symbols are encoded, so each part
// gets its own compact symbol map (0..base-1); the window key is the base-`base`
// polynomial packed into a 64-bit integer (safe: base is the number of distinct
// bytes in the part, so base^wsize fits in uint64 for any realistic sequence).
// Parts shorter than `wsize` return NA_real_ so the caller can reproduce the
// exact (degenerate) R value for those rare tail parts.
// [[Rcpp::export]]
NumericVector repeatScoreVecCpp(CharacterVector s, int wsize = 11, int dup_cutoff = 2){
  R_xlen_t ns = s.size();
  NumericVector out(ns);
  std::vector<unsigned long long> keys, pw((size_t)wsize);
  for(R_xlen_t si = 0; si < ns; ++si){
    if(STRING_ELT(s, si) == NA_STRING){ out[si] = NA_REAL; continue; }
    const char* str = CHAR(STRING_ELT(s, si));
    R_xlen_t L = (R_xlen_t) std::strlen(str);
    R_xlen_t nw = L - wsize + 1;
    if(nw < 1){ out[si] = NA_REAL; continue; }

    // Per-part symbol map -> small codes 0..base-1
    int lut[256]; for(int i = 0; i < 256; i++) lut[i] = -1;
    int base = 0;
    for(R_xlen_t i = 0; i < L; i++){
      unsigned char c = (unsigned char) str[i];
      if(lut[c] < 0) lut[c] = base++;
    }
    pw[0] = 1ULL;
    for(int j = 1; j < wsize; j++) pw[j] = pw[j-1] * (unsigned long long) base;

    // Window keys
    keys.clear(); keys.reserve((size_t) nw);
    for(R_xlen_t i = 0; i < nw; i++){
      unsigned long long k = 0ULL;
      for(int j = 0; j < wsize; j++)
        k += (unsigned long long) lut[(unsigned char) str[i + j]] * pw[j];
      keys.push_back(k);
    }
    std::sort(keys.begin(), keys.end());

    // Sum run lengths of groups larger than dup_cutoff
    long long num = 0;
    R_xlen_t i = 0, m = (R_xlen_t) keys.size();
    while(i < m){
      R_xlen_t j = i + 1;
      while(j < m && keys[j] == keys[i]) ++j;
      long long cnt = (long long)(j - i);
      if(cnt > dup_cutoff) num += cnt;
      i = j;
    }
    out[si] = (double) num / (double) nw;
  }
  return out;
}
