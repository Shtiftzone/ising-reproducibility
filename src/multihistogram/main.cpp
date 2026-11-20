// Circular moving block bootstrap for the 2D Ising multihistogram.
// Supports square and triangular lattices.
// Analyzes |EC_sym|/L^2, |M|/L^2, and |D|/L^2.
// Uses D = M + 2*EC_sym (square) or D = M + EC_sym (triangular).
// Compatible output columns with the original program.
// Usage: program temperatures_file L dimension results_file T_min T_max n_points output_file input_data_dir lattice_type// Optional environment variables: ISING_BLOCK_FACTOR (default 10),
// ISING_BOOTSTRAP_SEED (default 123456), ISING_BOOTSTRAP_REPS (default 200).
#include "MultiHistRW.h"
#include "Autocorr.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace std;

static double average(const vector<double>& v) {
    if (v.empty()) throw runtime_error("Cannot average an empty series");
    double sum = 0.0;
    for (double x : v) sum += x;
    return sum / static_cast<double>(v.size());
}

static double positive_env_double(const char* name, double fallback) {
    const char* s = getenv(name);
    if (!s) return fallback;
    char* end = nullptr;
    double v = strtod(s, &end);
    if (!end || *end != '\0' || !isfinite(v) || v <= 0.0)
        throw runtime_error(string("Invalid value for ") + name);
    return v;
}

static int positive_env_int(const char* name, int fallback) {
    const char* s = getenv(name);
    if (!s) return fallback;
    char* end = nullptr;
    long v = strtol(s, &end, 10);
    if (!end || *end != '\0' || v <= 0 || v > 10000000)
        throw runtime_error(string("Invalid value for ") + name);
    return static_cast<int>(v);
}

static unsigned long long rng_seed() {
    const char* s = getenv("ISING_BOOTSTRAP_SEED");
    if (!s) return 123456ULL;
    char* end = nullptr;
    unsigned long long v = strtoull(s, &end, 10);
    if (!end || *end != '\0' || *s == '-')
        throw runtime_error("Invalid ISING_BOOTSTRAP_SEED");
    return v;
}

static vector<double> squared(const vector<double>& x) {
    vector<double> result;
    result.reserve(x.size());
    for (double v : x) result.push_back(v * v);
    return result;
}

static double integrated_tau(const vector<double>& series, const string& name,
                             double T, int L) {
    if (series.size() < 4)
        throw runtime_error("Too few measurements for autocorrelation");
    const auto mm = minmax_element(series.begin(), series.end());
    if (*mm.first == *mm.second)
        throw runtime_error("Constant series in AutoCorr: " + name);
    double tau = 0.0, err = 0.0;
    AutoCorr(series, tau, err);
    if (!isfinite(tau) || !isfinite(err) || tau <= 0.0)
        throw runtime_error("AutoCorr failed for " + name + " at L=" +
                            to_string(L) + " T=" + to_string(T));
    return tau;
}

// Generates N indices from overlapping circular blocks of length block_len.
// Blocks are sampled with replacement. The final block is truncated to N.
static vector<size_t> circular_indices(size_t N, size_t block_len,
                                      mt19937_64& rng) {
    vector<size_t> indices;
    indices.reserve(N);
    uniform_int_distribution<size_t> start_dist(0, N - 1);
    while (indices.size() < N) {
        const size_t start = start_dist(rng);
        for (size_t k = 0; k < block_len && indices.size() < N; ++k)
            indices.push_back((start + k) % N);
    }
    return indices;
}

static vector<double> take(const vector<double>& data,
                           const vector<size_t>& indices) {
    vector<double> v;
    v.reserve(indices.size());
    for (size_t j : indices) v.push_back(data[j]);
    return v;
}

struct Series {
    double T = 0.0;
    vector<double> E, ec, mag, d;
    size_t block_len = 0;
    double tau_max = 0.0;
    string slowest;
};

int main(int argc, char** argv) {
    try {
        if (argc != 11) {
            cerr << "Usage: " << argv[0]
                 << " temperatures_file L dimension results_file T_min T_max n_points"
                    " output_file input_data_dir lattice_type\n";
            return 1;
        }
        const int L = stoi(argv[2]);
        const int dim = stoi(argv[3]);
        // Match the original repository interface: T_min followed by T_max.
        const double Tmin = stod(argv[5]);
        const double Tmax = stod(argv[6]);
        const int npoints = stoi(argv[7]);
        const double block_factor = positive_env_double("ISING_BLOCK_FACTOR", 10.0);
        const int nrep = positive_env_int("ISING_BOOTSTRAP_REPS", 200);
        const auto seed = rng_seed();
        if (L <= 0 || dim != 2 || npoints < 2 ||
            !(Tmin > 0.0 && Tmax > Tmin))
            throw runtime_error("Require L>0, dimension=2, n_points>=2, 0<T_min<T_max");
        const double vol = pow(static_cast<double>(L), dim);
        const string base = argv[9];
        const string lattice_type = argv[10];
        
        if (lattice_type != "square" && lattice_type != "triangular")
            throw runtime_error("lattice_type must be square or triangular");

        const double d_factor = (lattice_type == "square") ? 2.0 : 1.0;
        ifstream input(argv[1]);
        if (!input) throw runtime_error(string("Cannot open temperatures file: ") + argv[1]);
        ofstream outf(argv[4]), outf0(argv[8]);
        if (!outf || !outf0) throw runtime_error("Cannot open output file(s)");
        vector<Series> series;
        string line;
        while (getline(input, line)) {
            if (line.find_first_not_of(" \t\r") == string::npos ||
                line[line.find_first_not_of(" \t\r")] == '#') continue;
            istringstream ss(line);
            double T; int run;
            if (!(ss >> T >> run) || T <= 0.0)
                throw runtime_error("Invalid temperatures file row: " + line);
            ostringstream filename;
            filename << base << "/size_" << L << "/conf-"
                     << fixed << setprecision(4) << T << "-" << run << ".dat";
            ifstream in(filename.str());
            if (!in) throw runtime_error("Cannot open " + filename.str());
            Series s;
            s.T = T;
            double M, E, chip, chin, f0p, f1p, f2p, f0n, f1n, f2n;
            while (in >> M >> E >> chip >> chin >> f0p >> f1p >> f2p
                      >> f0n >> f1n >> f2n) {
                      // ec_sym = (chi_plus - chi_minus) / 2.
                      // d = M + 2*ec_sym (square) or M + ec_sym (triangular).
                      // Store absolute observables per lattice site.
                      const double q = (chip - chin) / 2.0;
                      const double d = M + d_factor * q;
                s.E.push_back(E);
                s.ec.push_back(abs(q) / vol);
                s.mag.push_back(abs(M) / vol);
                s.d.push_back(abs(d) / vol);
            }
            if (!in.eof()) throw runtime_error("Malformed measurements in " + filename.str());
            if (s.E.size() < 4) throw runtime_error("Too few measurements in " + filename.str());
            series.push_back(std::move(s));
        }
        if (series.empty()) throw runtime_error("No input series");
        vector<double> beta(series.size(), 0.0);
        vector<double> logZ(series.size(), 0.0);
        for (size_t i = 0; i < series.size(); ++i) {
            Series& s = series[i];
            beta[i] = 1.0 / s.T;
            // Estimate autocorrelation times for the energy, three absolute
            // observables, and their squares. Use the slowest time scale
            // to choose a single block length for all measurements at T.
            const vector<pair<string, vector<double>>> vars = {
                {"energy", s.E}, {"abs_m", s.mag}, {"abs_ec", s.ec},
                {"abs_d", s.d}, {"m_squared", squared(s.mag)},
                {"ec_squared", squared(s.ec)}, {"d_squared", squared(s.d)}
            };
            s.tau_max = 0.0;
            for (const auto& kv : vars) {
                double tau = integrated_tau(kv.second, kv.first, s.T, L);
                if (tau > s.tau_max) { s.tau_max = tau; s.slowest = kv.first; }
            }
            // Block length = ceil(block_factor * max(0.5, tau_max)).
            // Report the effective number of blocks below as a diagnostic;
            // do not tune the block length to enforce a block-count threshold.
            const double raw_len = ceil(block_factor * max(0.5, s.tau_max));
            if (raw_len > static_cast<double>(numeric_limits<size_t>::max()))
                throw runtime_error("Block size overflow");
            s.block_len = static_cast<size_t>(raw_len);
            if (s.block_len > s.E.size())
                throw runtime_error("Block length exceeds series length");
            cerr << setprecision(10) << "L=" << L << " T=" << s.T
                 << " N=" << s.E.size() << " tau_max=" << s.tau_max
                 << " slowest=" << s.slowest << " block=" << s.block_len
                 << " blocks_per_replicate="
                 << (s.E.size() + s.block_len - 1) / s.block_len << "\n";
        }
        cerr << "Circular block bootstrap: seed=" << seed
             << " reps=" << nrep << " factor=" << block_factor << "\n";
        mt19937_64 rng(seed);
        // Keep both output files and column order identical to the old main.
        // Each replicate uses independently sampled blocks at each temperature.
        for (int rep = 0; rep < nrep; ++rep) {
            vector<vector<double>> ebs(series.size()), obs(series.size()),
                                   mbs(series.size()), dbs(series.size());
            for (size_t i = 0; i < series.size(); ++i) {
                const Series& s = series[i];
                // Reuse the same indices for energy and all observables to
                // retain correlations within each saved configuration.
                const auto indices = circular_indices(s.E.size(), s.block_len, rng);
                ebs[i] = take(s.E, indices);
                obs[i] = take(s.ec, indices);
                mbs[i] = take(s.mag, indices);
                dbs[i] = take(s.d, indices);
                // Bootstrap estimates at the original simulation temperatures.
                const double ecm = average(obs[i]);
                const double mm = average(mbs[i]);
                const double dm = average(dbs[i]);
                const double ecvar = average(squared(obs[i])) - ecm*ecm;
                const double mvar = average(squared(mbs[i])) - mm*mm;
                const double dvar = average(squared(dbs[i])) - dm*dm;
                outf0 << setprecision(10) << rep << ' ' << s.T << ' '
                      << ecm << ' ' << ecvar*vol << ' '
                      << mm << ' ' << mvar*vol << ' '
                      << dm << ' ' << dvar*vol << '\n';
            }
            // Multihistogram reweighting using this same bootstrap replicate.
            // The existing multihistogram functions expect arrays of vectors.
            MultiHistRw(static_cast<int>(series.size()), beta.data(), ebs.data(), logZ.data());
            // Beta decreases as temperature increases.
            const double bmin = 1.0 / Tmax, bmax = 1.0 / Tmin;
            for (int j = 0; j < npoints; ++j) {
                const double b = bmin + (bmax-bmin)*j/static_cast<double>(npoints-1);
                const int nr = static_cast<int>(series.size());
                const double freeNRG = LnZ(b, nr, beta.data(), ebs.data(), logZ.data());
                const double ecmean = exp(LnO_n(obs.data(), 1.0, b, nr, beta.data(), ebs.data(), logZ.data()) - freeNRG);
                const double ecvar = exp(LnO_n(obs.data(), 2.0, b, nr, beta.data(), ebs.data(), logZ.data()) - freeNRG) - ecmean*ecmean;
                const double mmean = exp(LnO_n(mbs.data(), 1.0, b, nr, beta.data(), ebs.data(), logZ.data()) - freeNRG);
                const double mvar = exp(LnO_n(mbs.data(), 2.0, b, nr, beta.data(), ebs.data(), logZ.data()) - freeNRG) - mmean*mmean;
                const double dmean = exp(LnO_n(dbs.data(), 1.0, b, nr, beta.data(), ebs.data(), logZ.data()) - freeNRG);
                const double dvar = exp(LnO_n(dbs.data(), 2.0, b, nr, beta.data(), ebs.data(), logZ.data()) - freeNRG) - dmean*dmean;
                outf << setprecision(10) << rep << ' ' << 1.0/b << ' ' << freeNRG << ' '
                     << ecmean << ' ' << ecvar*vol << ' '
                     << mmean << ' ' << mvar*vol << ' '
                     << dmean << ' ' << dvar*vol << '\n';
            }
        }
        return 0;
    } catch (const exception& e) {
        cerr << "ERROR: " << e.what() << '\n';
        return 1;
    }
}
