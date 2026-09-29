#include <Rcpp.h>
#include <vector>
#include <cstring>
using namespace Rcpp;

// Dotplot hits, matching the pure-R mxComp(toupper(seq2mx(s1, wsize)),
// toupper(seq2mx(s2, wsize)), wsize, nmatch) exactly, including the
// column-major order of which(arr.ind = TRUE).
// score(i, j) = number of positions k in 0..wsize-1 with s1[i+k] == s2[j+k],
// counting only A/C/G/T (case-insensitive). Instead of the dense
// n1 x n2 matrix, one column of scores is kept and slid along the diagonal:
//   score(i, j) = score(i-1, j-1) + m(i+wsize-1, j+wsize-1) - m(i-1, j-1),
// which is O(n1 * n2) time and O(n1 + hits) memory.
// Returns list(row, col, values) with 1-based numeric row/col.
// [[Rcpp::export]]
List dotHitsCpp(CharacterVector s1, CharacterVector s2, int wsize, int nmatch){
  R_xlen_t L1 = s1.size(), L2 = s2.size();
  R_xlen_t n1 = L1 - wsize + 1, n2 = L2 - wsize + 1;
  if(wsize < 1 || n1 < 1 || n2 < 1) stop("Sequences must be at least wsize long");

  // Codes: A/C/G/T -> 0..3; anything else never matches (different per sequence)
  auto encode = [](CharacterVector s, int other){
    std::vector<int> x(s.size());
    for(R_xlen_t i = 0; i < s.size(); i++){
      const char* c = CHAR(STRING_ELT(s, i));
      switch(c[0]){
        case 'A': case 'a': x[i] = 0; break;
        case 'C': case 'c': x[i] = 1; break;
        case 'G': case 'g': x[i] = 2; break;
        case 'T': case 't': x[i] = 3; break;
        default: x[i] = other;
      }
      if(c[0] != '\0' && c[1] != '\0') x[i] = other;
    }
    return x;
  };
  std::vector<int> a = encode(s1, -1), b = encode(s2, -2);

  std::vector<int> score(n1);
  std::vector<double> out_row, out_col, out_val;

  for(R_xlen_t j = 0; j < n2; j++){
    if(j == 0){
      for(R_xlen_t i = 0; i < n1; i++){
        int sc = 0;
        for(int k = 0; k < wsize; k++) sc += (a[i + k] == b[k]);
        score[i] = sc;
      }
    } else {
      // Descending i: score[i-1] still holds column j-1
      for(R_xlen_t i = n1 - 1; i >= 1; i--){
        score[i] = score[i - 1] + (a[i + wsize - 1] == b[j + wsize - 1]) - (a[i - 1] == b[j - 1]);
      }
      int sc = 0;
      for(int k = 0; k < wsize; k++) sc += (a[k] == b[j + k]);
      score[0] = sc;
    }
    for(R_xlen_t i = 0; i < n1; i++){
      if(score[i] >= nmatch && score[i] != 0){
        out_row.push_back((double)(i + 1));
        out_col.push_back((double)(j + 1));
        out_val.push_back((double)score[i]);
      }
    }
    if((j & 1023) == 0) checkUserInterrupt();
  }

  return List::create(_["row"] = wrap(out_row),
                      _["col"] = wrap(out_col),
                      _["values"] = wrap(out_val));
}
