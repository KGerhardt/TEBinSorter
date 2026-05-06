// WFA_TEsorter
//
// Aligns multi-sequence query FASTA against multi-sequence target FASTA
// using the WaveFront Alignment library (WFA2-lib) and emits PAF format.
//
// Convention used here: WFA `pattern` is the QUERY, WFA `text` is the TARGET.
// WFA emits ops where 'I' consumes text (target) and 'D' consumes pattern
// (query). To produce SAM/PAF-compliant CIGAR and cs strings (where 'I'
// consumes query and 'D' consumes target), we invert I<->D when serialising.
//
// Build (one line; static-linked against WFA2-lib so the binary has no
// runtime libwfa2*.so dependency). Run from this src/ directory; assumes
// WFA2-lib sits at ./WFA2-lib with build/ already populated (libwfa2.a,
// libwfa2cpp.a). Adjust the -I and .a paths if your layout differs.
//   g++ -O3 -std=c++17 -fopenmp -Wall -I./WFA2-lib WFA_TEsorter.cpp -o WFA_TEsorter ./WFA2-lib/build/libwfa2cpp.a ./WFA2-lib/build/libwfa2.a -lpthread -lm

#include <algorithm>
#include <atomic>
#include <cctype>
#include <cerrno>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>
#include <mutex>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include <getopt.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "bindings/cpp/WFAligner.hpp"

using std::string;
using std::vector;
using std::ostream;
using std::cerr;
using std::cout;
using std::endl;

// ---------- FASTA -------------------------------------------------------------

struct FastaRec {
    string name;     // first whitespace-delimited token of the header
    string desc;     // remainder of the header (may be empty)
    string seq;      // upper-cased sequence
};

static void upper_inplace(string& s) {
    for (char& c : s) c = (char)std::toupper((unsigned char)c);
}

// Read a FASTA file into a vector of records. Sequence is upper-cased.
// Multi-line sequences and blank lines are handled. Errors abort the program.
static vector<FastaRec> read_fasta(const string& path) {
    std::ifstream in(path);
    if (!in) {
        cerr << "[ERROR] cannot open FASTA: " << path << " (" << std::strerror(errno) << ")\n";
        std::exit(1);
    }
    vector<FastaRec> recs;
    string line;
    FastaRec cur;
    bool have_cur = false;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (have_cur) { upper_inplace(cur.seq); recs.push_back(std::move(cur)); cur = FastaRec{}; }
            // parse name and description
            size_t p = 1;
            while (p < line.size() && !std::isspace((unsigned char)line[p])) ++p;
            cur.name = line.substr(1, p - 1);
            while (p < line.size() && std::isspace((unsigned char)line[p])) ++p;
            cur.desc = (p < line.size()) ? line.substr(p) : string();
            have_cur = true;
        } else {
            // strip trailing CR (windows line endings)
            while (!line.empty() && (line.back() == '\r' || line.back() == '\n')) line.pop_back();
            cur.seq.append(line);
        }
    }
    if (have_cur) { upper_inplace(cur.seq); recs.push_back(std::move(cur)); }
    if (recs.empty()) {
        cerr << "[ERROR] no sequences found in " << path << "\n";
        std::exit(1);
    }
    return recs;
}

// ---------- reverse complement ------------------------------------------------

static char comp_base(char c) {
    switch (c) {
        case 'A': return 'T'; case 'T': return 'A';
        case 'C': return 'G'; case 'G': return 'C';
        case 'U': return 'A';
        case 'R': return 'Y'; case 'Y': return 'R';
        case 'S': return 'S'; case 'W': return 'W';
        case 'K': return 'M'; case 'M': return 'K';
        case 'B': return 'V'; case 'V': return 'B';
        case 'D': return 'H'; case 'H': return 'D';
        case 'N': return 'N';
        case 'a': return 't'; case 't': return 'a';
        case 'c': return 'g'; case 'g': return 'c';
        case 'u': return 'a';
        case 'n': return 'n';
        default:  return 'N';
    }
}

static string revcomp(const string& s) {
    string out;
    out.resize(s.size());
    for (size_t i = 0; i < s.size(); ++i) out[s.size() - 1 - i] = comp_base(s[i]);
    return out;
}

// ---------- IUPAC nucleotide matching ----------------------------------------
//
// Each base is encoded as a 4-bit mask over {A=1, C=2, G=4, T=8}. Two bases
// "match" iff their masks share at least one bit, i.e. their possible-base
// sets intersect. Examples: A vs W={A,T} -> match; W vs S={C,G} -> mismatch;
// N matches everything; anything not in IUPAC has mask 0 and matches nothing.
// Used by the WFA lambda-match path (alignEnd2End / alignEndsFree /
// alignExtension lambda overloads) so the wavefront treats IUPAC-compatible
// positions as exact matches rather than mismatches.

static uint8_t g_base_mask[256];      // [byte] -> 4-bit base set
static uint8_t g_iupac_match[65536];  // [(qbyte<<8)|tbyte] -> 0/1 match flag

static void init_iupac_tables() {
    static bool inited = false;
    if (inited) return;
    inited = true;
    for (int i = 0; i < 256; ++i) g_base_mask[i] = 0;
    auto set_pair = [](char up, char lo, uint8_t m) {
        g_base_mask[(unsigned char)up] = m;
        g_base_mask[(unsigned char)lo] = m;
    };
    set_pair('A','a', 0x1);
    set_pair('C','c', 0x2);
    set_pair('G','g', 0x4);
    set_pair('T','t', 0x8);
    set_pair('U','u', 0x8);          // RNA U == T
    set_pair('R','r', 0x1|0x4);      // A|G
    set_pair('Y','y', 0x2|0x8);      // C|T
    set_pair('S','s', 0x2|0x4);      // C|G
    set_pair('W','w', 0x1|0x8);      // A|T
    set_pair('K','k', 0x4|0x8);      // G|T
    set_pair('M','m', 0x1|0x2);      // A|C
    set_pair('B','b', 0x2|0x4|0x8);      // C|G|T
    set_pair('D','d', 0x1|0x4|0x8);      // A|G|T
    set_pair('H','h', 0x1|0x2|0x8);      // A|C|T
    set_pair('V','v', 0x1|0x2|0x4);      // A|C|G
    set_pair('N','n', 0xF);              // any
    for (int a = 0; a < 256; ++a) {
        for (int b = 0; b < 256; ++b) {
            g_iupac_match[(a << 8) | b] = (g_base_mask[a] & g_base_mask[b]) ? 1 : 0;
        }
    }
}

static inline bool iupac_compatible(char a, char b) {
    return g_iupac_match[((unsigned char)a << 8) | (unsigned char)b] != 0;
}

// Context passed to the WFA lambda match callback.
struct LambdaCtx {
    const char* pattern;   // query bytes, [0, plen)
    const char* text;      // target bytes, [0, tlen)
};

// WFA calls this for each (v=pattern_pos, h=text_pos) probed during extension.
// Return non-zero iff the two bases are IUPAC-compatible; WFA then treats the
// position as a match (no penalty). The indices are guaranteed in-range by
// the caller (see wavefront_sequences_cmp), so no bounds check needed.
static int iupac_match_funct(int v, int h, void* args) {
    auto* ctx = static_cast<LambdaCtx*>(args);
    return g_iupac_match[((unsigned char)ctx->pattern[v] << 8)
                       | (unsigned char)ctx->text[h]];
}

// ---------- alignment op-string utilities -------------------------------------
//
// `ops` is the raw WFA operation string of M/X/I/D, one char per column.
// Convention here: WFA pattern == query, WFA text == target.
//   WFA 'M': match, advances both
//   WFA 'X': mismatch, advances both
//   WFA 'I': advances text (target)  -> SAM-style 'D' (deletion from query)
//   WFA 'D': advances pattern (query) -> SAM-style 'I' (insertion in query)

struct AlnStats {
    int matches    = 0;   // count of M (or = if show_eq) operations
    int mismatches = 0;   // X operations
    int ins_q      = 0;   // bases inserted in query (gap in target) = WFA 'D'
    int del_q      = 0;   // bases deleted from query (gap in query, base in target) = WFA 'I'
    int aln_len    = 0;   // total ops length (alignment block length)
    int q_consumed = 0;   // bases of query consumed (M+X+ins_q)
    int t_consumed = 0;   // bases of target consumed (M+X+del_q)
};

static AlnStats compute_stats(const string& ops) {
    AlnStats s;
    s.aln_len = (int)ops.size();
    for (char c : ops) {
        switch (c) {
            case 'M': ++s.matches; ++s.q_consumed; ++s.t_consumed; break;
            case 'X': ++s.mismatches; ++s.q_consumed; ++s.t_consumed; break;
            case 'I': ++s.del_q;     ++s.t_consumed; break;     // gap in query
            case 'D': ++s.ins_q;     ++s.q_consumed; break;     // gap in target
            default: break;
        }
    }
    return s;
}

// Run-length-encode WFA ops into a SAM CIGAR string with I<->D inversion so
// that the result describes query-vs-target in the standard way.
//   show_eq=false: emit M for both M and X (legacy SAM)
//   show_eq=true:  emit '=' for M and 'X' for X (extended SAM)
static string ops_to_cigar(const string& ops, bool show_eq) {
    string out;
    if (ops.empty()) return out;
    auto sam_op = [show_eq](char w) -> char {
        switch (w) {
            case 'M': return show_eq ? '=' : 'M';
            case 'X': return show_eq ? 'X' : 'M';
            case 'I': return 'D';   // WFA I -> SAM D
            case 'D': return 'I';   // WFA D -> SAM I
        }
        return '?';
    };
    char run_op = sam_op(ops[0]);
    int  run_len = 1;
    for (size_t i = 1; i < ops.size(); ++i) {
        char op = sam_op(ops[i]);
        if (op == run_op) { ++run_len; continue; }
        out += std::to_string(run_len); out += run_op;
        run_op = op; run_len = 1;
    }
    out += std::to_string(run_len); out += run_op;
    return out;
}

// Build minimap2-style cs tag.
//   short (default):  :LEN  *xy  +bases  -bases
//   long:             =BASES *xy +bases  -bases
// `query` and `target` are the (oriented) full sequences corresponding to the
// alignment. q0/t0 are the alignment start offsets in those sequences.
static string ops_to_cs(const string& ops,
                        const string& query, int q0,
                        const string& target, int t0,
                        bool long_form) {
    string out;
    int qi = q0, ti = t0;
    size_t i = 0;
    while (i < ops.size()) {
        char w = ops[i];
        if (w == 'M') {
            size_t j = i;
            while (j < ops.size() && ops[j] == 'M') ++j;
            int run = (int)(j - i);
            if (long_form) {
                out += '=';
                out.append(query, qi, run);
            } else {
                out += ':';
                out += std::to_string(run);
            }
            qi += run; ti += run; i = j;
        } else if (w == 'X') {
            // emit each mismatch separately as *xy where x=target base, y=query base
            // (lowercase per minimap2 short-cs convention)
            char tbase = (char)std::tolower((unsigned char)target[ti]);
            char qbase = (char)std::tolower((unsigned char)query[qi]);
            out += '*'; out += tbase; out += qbase;
            ++qi; ++ti; ++i;
        } else if (w == 'D') {
            // WFA 'D' -> SAM 'I' -> +bases (insertion to target = bases from query)
            size_t j = i;
            while (j < ops.size() && ops[j] == 'D') ++j;
            int run = (int)(j - i);
            out += '+';
            for (int k = 0; k < run; ++k) {
                out += (char)std::tolower((unsigned char)query[qi + k]);
            }
            qi += run; i = j;
        } else if (w == 'I') {
            // WFA 'I' -> SAM 'D' -> -bases (deletion from query = bases from target)
            size_t j = i;
            while (j < ops.size() && ops[j] == 'I') ++j;
            int run = (int)(j - i);
            out += '-';
            for (int k = 0; k < run; ++k) {
                out += (char)std::tolower((unsigned char)target[ti + k]);
            }
            ti += run; i = j;
        } else {
            ++i; // unknown op, skip
        }
    }
    return out;
}

// ---------- PAF record --------------------------------------------------------

struct PafRec {
    string qname; int qlen, qstart, qend;
    char strand;
    string tname; int tlen, tstart, tend;
    int matches;
    int aln_len;
    int mapq;
    int score;        // WFA penalty (typically <= 0)
    int nm;           // edit distance: X+I+D
    string cigar;     // optional; empty if not requested
    string cs;        // optional; empty if not requested
};

static void write_paf(ostream& os, const PafRec& p) {
    os << p.qname << '\t' << p.qlen << '\t' << p.qstart << '\t' << p.qend << '\t'
       << p.strand << '\t'
       << p.tname << '\t' << p.tlen << '\t' << p.tstart << '\t' << p.tend << '\t'
       << p.matches << '\t' << p.aln_len << '\t' << p.mapq
       << "\ttp:A:P"
       << "\tAS:i:" << p.score
       << "\tNM:i:" << p.nm;
    if (!p.cigar.empty()) os << "\tcg:Z:" << p.cigar;
    if (!p.cs.empty())    os << "\tcs:Z:" << p.cs;
    os << '\n';
}

// ---------- option parsing ----------------------------------------------------

enum class CsMode { None, Short, Long };
enum class AlnMode { End2End, EndsFree, Extension };
enum class Strand  { Both, Forward, Reverse };
enum class MemMode { High, Med, Low, Ultralow };

struct Opts {
    string query_path;
    string target_path;
    string output_path;     // empty => stdout
    string pairs_path;      // optional; empty => default all-vs-all behavior
    int threads     = 1;
    int match_pen   = 0;
    int mismatch    = 4;
    int gap_open    = 6;
    int gap_ext     = 2;
    bool emit_cigar = false;
    CsMode cs       = CsMode::None;
    bool show_eq    = false;
    AlnMode mode    = AlnMode::End2End;
    Strand strand   = Strand::Both;
    MemMode mem     = MemMode::High;
    bool has_min_score = false;
    int min_score   = 0;
    bool iupac      = false;
    bool verbose    = false;
};

static void usage_align(FILE* f = stderr) {
    std::fprintf(f,
"Usage: WFA_TEsorter -q QUERY.fa -t TARGET.fa [options] > out.paf\n"
"       WFA_TEsorter view --paf F.paf -q Q.fa -t T.fa --qname X --tname Y\n"
"\n"
"Required:\n"
"  -q, --query FILE       multi-seq query FASTA\n"
"  -t, --target FILE      multi-seq target FASTA\n"
"\n"
"Output:\n"
"  -o, --output FILE      PAF output (default: stdout)\n"
"      --cigar            emit CIGAR as cg:Z tag\n"
"      --cs [=short|long] emit minimap2 cs tag (default: short)\n"
"      --show-eq          use '='/'X' instead of 'M' in CIGAR\n"
"\n"
"Penalties (gap-affine, all positive integers):\n"
"      --match N          match score (default: 0)\n"
"      --mismatch N       mismatch penalty (default: 4)\n"
"      --gap-open N       gap-open penalty (default: 6)\n"
"      --gap-ext N        gap-extension penalty (default: 2)\n"
"\n"
"Alignment:\n"
"      --mode STR         end2end|ends-free|extension (default: end2end)\n"
"                         ends-free: query is anchored, target ends are free\n"
"      --strand STR       both|forward|reverse (default: both)\n"
"      --memory STR       high|med|low|ultralow (default: high)\n"
"      --min-score N      drop alignments with WFA score < N (default: keep all)\n"
"      --iupac            treat IUPAC-compatible bases as matches (e.g. A==W,\n"
"                         W==N, R==G). Off by default; turn on to align\n"
"                         consensus sequences carrying ambiguity codes. Uses\n"
"                         WFA's lambda-match path, which is slower than the\n"
"                         default literal-byte kernel.\n"
"      --pairs FILE       restrict to (qname<TAB>tname) pairs from FILE\n"
"                         (default: full all-vs-all over the FASTAs)\n"
"                         Self-pairs (qname == tname) are always skipped.\n"
"      --threads N        worker threads (default: 1)\n"
"\n"
"Misc:\n"
"  -v, --verbose          progress messages to stderr\n"
"  -h, --help             this message\n"
"\n"
"Subcommands:\n"
"  view                   pretty-print an alignment from a PAF (requires cg:Z)\n"
"                         see: WFA_TEsorter view -h\n"
);
}

static MemMode parse_mem(const string& s) {
    if (s == "high")     return MemMode::High;
    if (s == "med")      return MemMode::Med;
    if (s == "low")      return MemMode::Low;
    if (s == "ultralow") return MemMode::Ultralow;
    cerr << "[ERROR] unknown --memory value: " << s << "\n"; std::exit(2);
}
static AlnMode parse_mode(const string& s) {
    if (s == "end2end" || s == "end-to-end") return AlnMode::End2End;
    if (s == "ends-free")                    return AlnMode::EndsFree;
    if (s == "extension")                    return AlnMode::Extension;
    cerr << "[ERROR] unknown --mode value: " << s << "\n"; std::exit(2);
}
static Strand parse_strand(const string& s) {
    if (s == "both")    return Strand::Both;
    if (s == "forward" || s == "+") return Strand::Forward;
    if (s == "reverse" || s == "-") return Strand::Reverse;
    cerr << "[ERROR] unknown --strand value: " << s << "\n"; std::exit(2);
}

static Opts parse_align_opts(int argc, char** argv) {
    Opts o;
    static struct option long_opts[] = {
        {"query",      required_argument, nullptr, 'q'},
        {"target",     required_argument, nullptr, 't'},
        {"output",     required_argument, nullptr, 'o'},
        {"threads",    required_argument, nullptr,  1 },
        {"match",      required_argument, nullptr,  2 },
        {"mismatch",   required_argument, nullptr,  3 },
        {"gap-open",   required_argument, nullptr,  4 },
        {"gap-ext",    required_argument, nullptr,  5 },
        {"cigar",      no_argument,       nullptr,  6 },
        {"cs",         optional_argument, nullptr,  7 },
        {"show-eq",    no_argument,       nullptr,  8 },
        {"mode",       required_argument, nullptr,  9 },
        {"strand",     required_argument, nullptr, 10 },
        {"memory",     required_argument, nullptr, 11 },
        {"min-score",  required_argument, nullptr, 12 },
        {"pairs",      required_argument, nullptr, 13 },
        {"iupac",      no_argument,       nullptr, 14 },
        {"verbose",    no_argument,       nullptr, 'v'},
        {"help",       no_argument,       nullptr, 'h'},
        {nullptr, 0, nullptr, 0}
    };
    int c;
    while ((c = getopt_long(argc, argv, "q:t:o:vh", long_opts, nullptr)) != -1) {
        switch (c) {
            case 'q': o.query_path  = optarg; break;
            case 't': o.target_path = optarg; break;
            case 'o': o.output_path = optarg; break;
            case  1 : o.threads    = std::atoi(optarg); break;
            case  2 : o.match_pen  = std::atoi(optarg); break;
            case  3 : o.mismatch   = std::atoi(optarg); break;
            case  4 : o.gap_open   = std::atoi(optarg); break;
            case  5 : o.gap_ext    = std::atoi(optarg); break;
            case  6 : o.emit_cigar = true; break;
            case  7 : o.cs = (optarg && string(optarg) == "long") ? CsMode::Long : CsMode::Short; break;
            case  8 : o.show_eq = true; break;
            case  9 : o.mode    = parse_mode(optarg); break;
            case 10 : o.strand  = parse_strand(optarg); break;
            case 11 : o.mem     = parse_mem(optarg); break;
            case 12 : o.has_min_score = true; o.min_score = std::atoi(optarg); break;
            case 13 : o.pairs_path = optarg; break;
            case 14 : o.iupac = true; break;
            case 'v': o.verbose = true; break;
            case 'h': usage_align(stdout); std::exit(0);
            default : usage_align(stderr); std::exit(2);
        }
    }
    if (o.query_path.empty() || o.target_path.empty()) {
        cerr << "[ERROR] --query and --target are required\n";
        usage_align(stderr); std::exit(2);
    }
    if (o.threads < 1) o.threads = 1;
    if (o.mismatch < 0 || o.gap_open < 0 || o.gap_ext < 0) {
        cerr << "[ERROR] penalties must be non-negative\n"; std::exit(2);
    }
    return o;
}

// ---------- alignment runner --------------------------------------------------

static wfa::WFAligner::MemoryModel to_wfa_mem(MemMode m) {
    switch (m) {
        case MemMode::High:     return wfa::WFAligner::MemoryHigh;
        case MemMode::Med:      return wfa::WFAligner::MemoryMed;
        case MemMode::Low:      return wfa::WFAligner::MemoryLow;
        case MemMode::Ultralow: return wfa::WFAligner::MemoryUltralow;
    }
    return wfa::WFAligner::MemoryHigh;
}

struct AlnResult {
    bool   ok     = false;
    string ops;             // WFA M/X/I/D string for the *aligned* region only
    int    score  = 0;
    int    qstart = 0, qend = 0;
    int    tstart = 0, tend = 0;
};

// Run one alignment. For ends-free (query-anchored) mode the leading and
// trailing runs of WFA 'I' ops (free target-end skips) are stripped from `ops`
// and accounted for in tstart/tend. For extension mode qend/tend reflect the
// actual bases consumed by the partial alignment.
//
// When `iupac` is true, the lambda-match overloads are used so that bases
// connected by IUPAC ambiguity codes (e.g. W vs A) score as matches rather
// than mismatches. The lambda path bypasses the SIMD/64-bit-block extension
// kernel and is therefore slower than the default literal-byte path.
static AlnResult align_one(wfa::WFAlignerGapAffine& aligner, AlnMode mode,
                           const string& q, const string& t, bool iupac) {
    AlnResult r;
    int qlen = (int)q.size(), tlen = (int)t.size();
    wfa::WFAligner::AlignmentStatus st;
    if (iupac) {
        LambdaCtx ctx{q.data(), t.data()};
        if (mode == AlnMode::End2End) {
            st = aligner.alignEnd2End(iupac_match_funct, &ctx, qlen, tlen);
        } else if (mode == AlnMode::EndsFree) {
            st = aligner.alignEndsFree(iupac_match_funct, &ctx,
                                       qlen, 0, 0,
                                       tlen, tlen, tlen);
        } else {
            st = aligner.alignExtension(iupac_match_funct, &ctx, qlen, tlen);
        }
    } else if (mode == AlnMode::End2End) {
        st = aligner.alignEnd2End(q, t);
    } else if (mode == AlnMode::EndsFree) {
        // Query-anchored semi-global: query is anchored, target ends are free.
        st = aligner.alignEndsFree(q, 0, 0, t, tlen, tlen);
    } else {
        string qm = q, tm = t;
        st = aligner.alignExtension(qm, tm);
    }
    if ((int)st < 0) return r;
    string ops = aligner.getAlignment();
    if (ops.empty()) return r;

    int q0 = 0, q1 = qlen, t0 = 0, t1 = tlen;
    if (mode == AlnMode::EndsFree) {
        // Strip leading/trailing 'I' ops (WFA: text-advance) which represent
        // free skips of the target ends.
        size_t lo = 0, hi = ops.size();
        while (lo < hi && ops[lo] == 'I')      { ++t0; ++lo; }
        while (hi > lo && ops[hi - 1] == 'I')  { --t1; --hi; }
        ops = ops.substr(lo, hi - lo);
    } else if (mode == AlnMode::Extension) {
        AlnStats s = compute_stats(ops);
        q1 = q0 + s.q_consumed;
        t1 = t0 + s.t_consumed;
    }

    if (ops.empty()) return r;
    r.ok    = true;
    r.ops   = std::move(ops);
    r.score = aligner.getAlignmentScore();
    r.qstart = q0; r.qend = q1;
    r.tstart = t0; r.tend = t1;
    return r;
}

// Build a PafRec from an AlnResult. The aligned ops describe positions
// [qstart,qend) on `q_oriented` and [tstart,tend) on `t_oriented`.
static PafRec make_paf(const string& qname, const string& q_oriented, int qlen_full,
                       const string& tname, const string& t_oriented, int tlen_full,
                       char strand, const AlnResult& r,
                       const Opts& o) {
    AlnStats s = compute_stats(r.ops);
    PafRec p;
    p.qname = qname;  p.qlen = qlen_full;  p.qstart = r.qstart; p.qend = r.qend;
    p.strand = strand;
    p.tname = tname;  p.tlen = tlen_full;  p.tstart = r.tstart; p.tend = r.tend;
    p.matches = s.matches;
    p.aln_len = s.aln_len;
    p.mapq    = 60;
    p.score   = r.score;
    p.nm      = s.mismatches + s.ins_q + s.del_q;
    if (o.emit_cigar) p.cigar = ops_to_cigar(r.ops, o.show_eq);
    if (o.cs != CsMode::None) {
        p.cs = ops_to_cs(r.ops, q_oriented, r.qstart,
                                  t_oriented, r.tstart,
                                  o.cs == CsMode::Long);
    }
    return p;
}

static int run_align(int argc, char** argv) {
    Opts o = parse_align_opts(argc, argv);

    if (o.iupac) init_iupac_tables();

    if (o.verbose) cerr << "[INFO] reading FASTAs\n";
    auto queries = read_fasta(o.query_path);
    auto targets = read_fasta(o.target_path);
    if (o.verbose) {
        cerr << "[INFO] " << queries.size() << " query seqs, "
             << targets.size() << " target seqs"
             << (o.iupac ? " (IUPAC-aware matching)" : "") << "\n";
    }

    // Output sink
    std::ofstream fout;
    ostream* out = &cout;
    if (!o.output_path.empty()) {
        fout.open(o.output_path);
        if (!fout) { cerr << "[ERROR] cannot open output: " << o.output_path << "\n"; return 1; }
        out = &fout;
    }
    std::mutex out_mu;

    // Per-thread aligners (created lazily inside the parallel region so each
    // thread owns its own state).
    int nthreads = o.threads;
#ifdef _OPENMP
    omp_set_num_threads(nthreads);
#else
    nthreads = 1;
#endif

    // Optional: pre-load explicit (qname, tname) pairs from --pairs FILE.
    // When pairs_path is empty, fall through to the default all-vs-all loop.
    std::vector<std::pair<int,int>> pair_jobs;
    if (!o.pairs_path.empty()) {
        std::unordered_map<string,int> q_idx, t_idx;
        q_idx.reserve(queries.size() * 2);
        t_idx.reserve(targets.size() * 2);
        for (size_t i = 0; i < queries.size(); ++i) q_idx.emplace(queries[i].name, (int)i);
        for (size_t i = 0; i < targets.size(); ++i) t_idx.emplace(targets[i].name, (int)i);

        std::ifstream pin(o.pairs_path);
        if (!pin) {
            cerr << "[ERROR] cannot open --pairs file: " << o.pairs_path << "\n";
            return 1;
        }
        size_t miss_q = 0, miss_t = 0, parsed = 0, skipped_self = 0;
        string line;
        while (std::getline(pin, line)) {
            if (line.empty() || line[0] == '#') continue;
            // strip trailing CR
            while (!line.empty() && (line.back() == '\r' || line.back() == '\n')) line.pop_back();
            size_t tab = line.find('\t');
            if (tab == string::npos) continue;
            string qn = line.substr(0, tab);
            string rest = line.substr(tab + 1);
            size_t tab2 = rest.find('\t');
            string tn = (tab2 == string::npos) ? rest : rest.substr(0, tab2);
            ++parsed;
            if (qn == tn) { ++skipped_self; continue; }
            auto qit = q_idx.find(qn);
            auto tit = t_idx.find(tn);
            if (qit == q_idx.end()) { ++miss_q; continue; }
            if (tit == t_idx.end()) { ++miss_t; continue; }
            pair_jobs.emplace_back(qit->second, tit->second);
        }
        if (o.verbose) {
            cerr << "[INFO] --pairs parsed=" << parsed
                 << " kept=" << pair_jobs.size()
                 << " skipped-self=" << skipped_self
                 << " missing-qname=" << miss_q
                 << " missing-tname=" << miss_t << "\n";
        }
        if (pair_jobs.empty()) {
            cerr << "[ERROR] no usable pairs from " << o.pairs_path << "\n";
            return 1;
        }
    }

    // Pre-count self-pairs in the all-vs-all matrix so the progress total
    // reflects only the alignments we'll actually run.
    size_t self_pairs_avs = 0;
    if (pair_jobs.empty()) {
        std::unordered_map<string,int> tname_counts;
        tname_counts.reserve(targets.size() * 2);
        for (const auto& T : targets) ++tname_counts[T.name];
        for (const auto& Q : queries) {
            auto it = tname_counts.find(Q.name);
            if (it != tname_counts.end()) self_pairs_avs += (size_t)it->second;
        }
    }
    std::atomic<size_t> done{0};
    const size_t total = pair_jobs.empty()
        ? (queries.size() * targets.size() - self_pairs_avs)
        : pair_jobs.size();
    auto t_start = std::chrono::steady_clock::now();

    // Process one (qi, ti) job using the caller's per-thread aligner.
    auto process_job = [&](size_t qi, size_t ti, wfa::WFAlignerGapAffine& aligner) {
        const FastaRec& Q = queries[qi];
        const FastaRec& T = targets[ti];
        if (Q.seq.empty() || T.seq.empty()) { ++done; return; }
        // Skip self-alignments: a sequence aligned to itself yields no useful
        // information and is wasted compute (typical case: all-vs-all over a
        // single FASTA loaded as both --query and --target). These are
        // pre-excluded from `total`, so do not increment `done` here.
        if (Q.name == T.name) return;

        AlnResult res_f, res_r;
        string qrev;
        if (o.strand == Strand::Forward || o.strand == Strand::Both) {
            res_f = align_one(aligner, o.mode, Q.seq, T.seq, o.iupac);
        }
        if (o.strand == Strand::Reverse || o.strand == Strand::Both) {
            qrev = revcomp(Q.seq);
            res_r = align_one(aligner, o.mode, qrev, T.seq, o.iupac);
        }

        bool emit_f = false, emit_r = false;
        if (o.strand == Strand::Forward) emit_f = res_f.ok;
        else if (o.strand == Strand::Reverse) emit_r = res_r.ok;
        else {
            if (res_f.ok && res_r.ok) {
                (res_f.score >= res_r.score) ? (emit_f = true) : (emit_r = true);
            } else if (res_f.ok) emit_f = true;
            else if (res_r.ok) emit_r = true;
        }
        if (emit_f && o.has_min_score && res_f.score < o.min_score) emit_f = false;
        if (emit_r && o.has_min_score && res_r.score < o.min_score) emit_r = false;

        if (emit_f) {
            PafRec p = make_paf(Q.name, Q.seq, (int)Q.seq.size(),
                                T.name, T.seq, (int)T.seq.size(),
                                '+', res_f, o);
            std::lock_guard<std::mutex> lk(out_mu);
            write_paf(*out, p);
        }
        if (emit_r) {
            PafRec p = make_paf(Q.name, qrev, (int)qrev.size(),
                                T.name, T.seq, (int)T.seq.size(),
                                '-', res_r, o);
            std::lock_guard<std::mutex> lk(out_mu);
            write_paf(*out, p);
        }

        size_t now = ++done;
        if (o.verbose && (now % 1000 == 0 || now == total)) {
            auto el = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - t_start).count();
            std::lock_guard<std::mutex> lk(out_mu);
            cerr << "[INFO] " << now << "/" << total
                 << " pairs aligned (" << (now / std::max(el, 1e-9)) << "/s)\n";
        }
    };

    if (pair_jobs.empty()) {
        // Default behavior: full all-vs-all over the two FASTAs.
#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
            wfa::WFAlignerGapAffine aligner(
                o.match_pen, o.mismatch, o.gap_open, o.gap_ext,
                wfa::WFAligner::Alignment, to_wfa_mem(o.mem));
#ifdef _OPENMP
            #pragma omp for schedule(dynamic, 1) collapse(2)
#endif
            for (size_t qi = 0; qi < queries.size(); ++qi) {
                for (size_t ti = 0; ti < targets.size(); ++ti) {
                    process_job(qi, ti, aligner);
                }
            }
        }
    } else {
        // --pairs mode: align only the explicit (qi, ti) pairs.
#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
            wfa::WFAlignerGapAffine aligner(
                o.match_pen, o.mismatch, o.gap_open, o.gap_ext,
                wfa::WFAligner::Alignment, to_wfa_mem(o.mem));
#ifdef _OPENMP
            #pragma omp for schedule(dynamic, 64)
#endif
            for (size_t pi = 0; pi < pair_jobs.size(); ++pi) {
                process_job(pair_jobs[pi].first, pair_jobs[pi].second, aligner);
            }
        }
    }

    if (o.verbose) {
        auto el = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - t_start).count();
        cerr << "[INFO] done in " << el << " s\n";
    }
    return 0;
}

// ---------- view subcommand ---------------------------------------------------

struct ViewOpts {
    string paf_path;
    string query_path;
    string target_path;
    string qname;
    string tname;
    int    width = 80;     // wrap width for the pretty alignment
    bool   verbose = false;
};

static void usage_view(FILE* f = stderr) {
    std::fprintf(f,
"Usage: WFA_TEsorter view --paf F.paf -q Q.fa -t T.fa --qname X --tname Y\n"
"\n"
"Re-creates a human-readable alignment from a PAF record produced by\n"
"WFA_TEsorter (the record must contain a cg:Z CIGAR tag).\n"
"\n"
"Required:\n"
"      --paf FILE         PAF file written by WFA_TEsorter --cigar\n"
"  -q, --query FILE       query FASTA used to produce the PAF\n"
"  -t, --target FILE      target FASTA used to produce the PAF\n"
"      --qname STR        query sequence name to display\n"
"      --tname STR        target sequence name to display\n"
"\n"
"Optional:\n"
"      --width N          characters per row (default: 80)\n"
"  -h, --help             this message\n"
);
}

static ViewOpts parse_view_opts(int argc, char** argv) {
    ViewOpts v;
    static struct option lo[] = {
        {"paf",     required_argument, nullptr, 1},
        {"query",   required_argument, nullptr, 'q'},
        {"target",  required_argument, nullptr, 't'},
        {"qname",   required_argument, nullptr, 2},
        {"tname",   required_argument, nullptr, 3},
        {"width",   required_argument, nullptr, 4},
        {"help",    no_argument,       nullptr, 'h'},
        {nullptr, 0, nullptr, 0}
    };
    int c;
    while ((c = getopt_long(argc, argv, "q:t:h", lo, nullptr)) != -1) {
        switch (c) {
            case  1 : v.paf_path = optarg; break;
            case 'q': v.query_path = optarg; break;
            case 't': v.target_path = optarg; break;
            case  2 : v.qname = optarg; break;
            case  3 : v.tname = optarg; break;
            case  4 : v.width = std::max(20, std::atoi(optarg)); break;
            case 'h': usage_view(stdout); std::exit(0);
            default : usage_view(stderr); std::exit(2);
        }
    }
    if (v.paf_path.empty() || v.query_path.empty() || v.target_path.empty() ||
        v.qname.empty() || v.tname.empty()) {
        cerr << "[ERROR] --paf, --query, --target, --qname, --tname all required\n";
        usage_view(stderr); std::exit(2);
    }
    return v;
}

// Parse a SAM/PAF CIGAR string into pairs (length, op).
// Recognises M, I, D, =, X. Returns false on parse error.
static bool parse_cigar(const string& cg, vector<std::pair<int,char>>& out) {
    out.clear();
    int i = 0, n = (int)cg.size();
    while (i < n) {
        int v = 0;
        if (!std::isdigit((unsigned char)cg[i])) return false;
        while (i < n && std::isdigit((unsigned char)cg[i])) {
            v = v * 10 + (cg[i] - '0');
            ++i;
        }
        if (i >= n) return false;
        char op = cg[i++];
        if (op != 'M' && op != 'I' && op != 'D' && op != '=' && op != 'X') return false;
        out.emplace_back(v, op);
    }
    return true;
}

// Find a tag in a tab-separated PAF line; returns its value (after "xx:T:") or
// empty string if absent.
static string get_tag(const string& line, const string& key /* e.g. "cg:Z:" */) {
    // search tab-delimited fields
    size_t pos = 0;
    while (pos < line.size()) {
        size_t end = line.find('\t', pos);
        size_t f_end = (end == string::npos) ? line.size() : end;
        if (f_end - pos >= key.size() &&
            line.compare(pos, key.size(), key) == 0) {
            return line.substr(pos + key.size(), f_end - pos - key.size());
        }
        if (end == string::npos) break;
        pos = end + 1;
    }
    return string();
}

// Build a vector of FastaRec keyed by name, returning a pointer (or nullptr).
static const FastaRec* find_rec(const vector<FastaRec>& v, const string& name) {
    for (const auto& r : v) if (r.name == name) return &r;
    return nullptr;
}

// Print the pretty alignment given the (oriented) query/target slices and CIGAR.
static void print_pretty(ostream& os,
                         const string& qname, const string& q_aln,
                         int q_off,
                         const string& tname, const string& t_aln,
                         int t_off,
                         char strand,
                         const vector<std::pair<int,char>>& cig,
                         int width) {
    // Build aligned strings (with '-' gaps).
    string qline, mline, tline;
    int qi = 0, ti = 0;
    for (auto [len, op] : cig) {
        for (int k = 0; k < len; ++k) {
            if (op == 'M' || op == '=' || op == 'X') {
                char qc = q_aln[qi++], tc = t_aln[ti++];
                qline += qc; tline += tc;
                // For legacy 'M' (match-or-mismatch), use IUPAC compatibility
                // so consensus bases like W vs A are drawn as matches when the
                // alignment was produced under --iupac. Pure ACGT pairs reduce
                // to literal equality.
                bool match = (op == '=') ||
                             (op == 'M' && iupac_compatible(qc, tc));
                mline += match ? '|' : ' ';
            } else if (op == 'I') {     // query has extra base (gap in target)
                qline += q_aln[qi++]; tline += '-'; mline += ' ';
            } else if (op == 'D') {     // target has extra base (gap in query)
                qline += '-'; tline += t_aln[ti++]; mline += ' ';
            }
        }
    }

    // Coordinate width is sized to fit the largest position printed.
    int q_max = q_off + qi;
    int t_max = t_off + ti;
    int pos_w = (int)std::max(std::to_string(q_max).size(),
                              std::to_string(t_max).size());

    auto pad_left = [](int v, int w) {
        std::string s = std::to_string(v);
        if ((int)s.size() < w) s.insert(0, w - s.size(), ' ');
        return s;
    };

    // Header
    os << "QUERY  : " << qname << "  (" << strand << " strand"
       << ", " << qi << " bp aligned starting at " << q_off << ")\n"
       << "TARGET : " << tname << "  ("
       << ti << " bp aligned starting at " << t_off << ")\n\n";

    // Each printed row has:  "QUERY  " (7) + <pos_w> + "  " (2) + alignment
    const int label_w = 7 + pos_w + 2;
    const std::string mid_pad(label_w, ' ');

    int n = (int)qline.size();
    int q_pos = q_off;
    int t_pos = t_off;
    for (int p = 0; p < n; p += width) {
        int len = std::min(width, n - p);
        int q_step = 0, t_step = 0;
        for (int k = 0; k < len; ++k) {
            if (qline[p + k] != '-') ++q_step;
            if (tline[p + k] != '-') ++t_step;
        }
        os << "QUERY  " << pad_left(q_pos + 1, pos_w) << "  "
           << qline.substr(p, len) << "  " << (q_pos + q_step) << "\n";
        os << mid_pad << mline.substr(p, len) << "\n";
        os << "TARGET " << pad_left(t_pos + 1, pos_w) << "  "
           << tline.substr(p, len) << "  " << (t_pos + t_step) << "\n\n";
        q_pos += q_step;
        t_pos += t_step;
    }
}

static int run_view(int argc, char** argv) {
    ViewOpts v = parse_view_opts(argc, argv);

    init_iupac_tables();

    auto qs = read_fasta(v.query_path);
    auto ts = read_fasta(v.target_path);
    const FastaRec* Q = find_rec(qs, v.qname);
    const FastaRec* T = find_rec(ts, v.tname);
    if (!Q) { cerr << "[ERROR] query name not found: " << v.qname << "\n"; return 1; }
    if (!T) { cerr << "[ERROR] target name not found: " << v.tname << "\n"; return 1; }

    std::ifstream pin(v.paf_path);
    if (!pin) { cerr << "[ERROR] cannot open PAF: " << v.paf_path << "\n"; return 1; }

    string line;
    int hits = 0;
    while (std::getline(pin, line)) {
        if (line.empty() || line[0] == '#') continue;
        // split first 12 fields
        vector<string> f;
        f.reserve(16);
        size_t pos = 0;
        while (pos <= line.size() && f.size() < 12) {
            size_t end = line.find('\t', pos);
            f.push_back(line.substr(pos, (end == string::npos) ? line.size() - pos : end - pos));
            if (end == string::npos) break;
            pos = end + 1;
        }
        if (f.size() < 12) continue;
        if (f[0] != v.qname || f[5] != v.tname) continue;

        char strand = f[4].empty() ? '+' : f[4][0];
        int  qstart = std::atoi(f[2].c_str());
        int  qend   = std::atoi(f[3].c_str());
        int  tstart = std::atoi(f[7].c_str());
        int  tend   = std::atoi(f[8].c_str());
        (void)qend; (void)tend;

        string cg = get_tag(line, "cg:Z:");
        if (cg.empty()) {
            cerr << "[WARN] PAF row for " << v.qname << " vs " << v.tname
                 << " has no cg:Z tag; skipping\n";
            continue;
        }
        vector<std::pair<int,char>> cig;
        if (!parse_cigar(cg, cig)) {
            cerr << "[WARN] could not parse CIGAR; skipping\n";
            continue;
        }
        // Build oriented query for this strand
        string q_oriented = (strand == '-') ? revcomp(Q->seq) : Q->seq;
        const string& t_oriented = T->seq;

        // Slice the aligned region from oriented sequences.
        string q_aln = q_oriented.substr(qstart, q_oriented.size() - qstart);
        string t_aln = t_oriented.substr(tstart, t_oriented.size() - tstart);

        ++hits;
        cout << "# alignment " << hits << " : " << v.qname << " ("
             << strand << ") vs " << v.tname << "\n";
        print_pretty(cout, v.qname, q_aln, qstart, v.tname, t_aln, tstart,
                     strand, cig, v.width);
    }
    if (hits == 0) {
        cerr << "[WARN] no PAF records found for " << v.qname << " vs " << v.tname << "\n";
        return 1;
    }
    return 0;
}

// ---------- entry -------------------------------------------------------------

int main(int argc, char** argv) {
    if (argc >= 2 && string(argv[1]) == "view") {
        // Shift past subcommand for getopt
        return run_view(argc - 1, argv + 1);
    }
    if (argc >= 2 && string(argv[1]) == "align") {
        return run_align(argc - 1, argv + 1);
    }
    return run_align(argc, argv);
}
