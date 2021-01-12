#include <filesystem>
#include <fstream>
#include <functional>
#include <random>

#include "utils.hpp"

#include "distances/dtw/dtw.hpp"
#include "distances/dtw/cdtw.hpp"
#include "distances/dtw/wdtw.hpp"
#include "distances/elementwise/elementwise.hpp"
#include "distances/erp/erp.hpp"
#include "distances/lcss/lcss.hpp"
#include "distances/msm/msm.hpp"
#include "distances/twe/twe.hpp"

#include "distances/dtw/lowerbounds/envelopes.hpp"
#include "distances/dtw/lowerbounds/lb_keogh.hpp"
#include "distances/dtw/lowerbounds/lb_webb.hpp"

#include "tseries/tseries.hpp"
#include "tseries/readers/tsreader/tsreader.hpp"

#define STRINGIFY_DETAIL(x) #x
#define STRINGIFY(x) STRINGIFY_DETAIL(x)

namespace fs = std::filesystem;
using namespace std;
using PRNG = mt19937_64;

/** Print the usage/help on cout */
void print_usage(const string &execname, ostream &out) {
    out << "NN1 classification for the UCR archive - Monash University, Melbourne, Australia, 2021" << endl;
    out << "Command:" << endl;
    out << "  " << execname << " <dataset> <distance> <base|eap|eap_la> [-out outfile]" << endl;
    out << "Dataset:" << endl;
    out << "  -ucr <path to dir>    Path to a directory containing datasets in ts format." << endl;
    out << "                        A dataset is itself a directory containing a _TRAIN.ts and a _TEST.ts files." << endl;
    out << "  -dataset <string>     Limit experiment to a given dataset. Must match the dataset's folder name." << endl;
    out << "Distance: -dist <distance name and args>:" << endl;
    out << "  dtw [LB]              DTW distance, optional lower bound" << endl;
    out << "  cdtw <wr> [LB]        CDTW distance with a window ratio 0<=wr<=1, optional lower bound" << endl;
    out << "    LB:" << endl;
    out << "      lb-none           Do not use any lower bound (default)" << endl;
    out << "      lb-keogh          LB-Keogh between the query and the envelopes of the candidate" << endl;
    out << "      lb-keogh2         LB-Keogh both way" << endl;
    out << "      lb-webb           LB-WEBB" << endl;
    out << "  wdtw <g>              WDTW distance with a weight factor 0<=g" << endl;
    out << "  sqed                  Squared Euclidean distance" << endl;
    out << "  erp <gv> <wr>         ERP distance with a gap value 0<=gv and window ratio 0<=wr<=1" << endl;
    out << "  lcss <e> <wr>         LCSS distance with an epsilon 0<=e and window ratio 0<=wr<=1" << endl;
    out << "  lcss <e> int <w>      LCSS distance with an epsilon 0<=e and integer window size 0<=w" << endl;
    out << "  msm <cost>            MSM distance with a Split and Merge cost 0<=cost" << endl;
    out << "  twe <nu> <lambda>     TWE distance with a stiffness 0<=nu and a constant delete penalty 0<=lambda" << endl;
    out << "Note: windows are computed using the window ratio and the longest series in the dataset." << endl;
    out << "Run mode:" << endl;
    out << "  base                  Base implementation, on a double buffer where applicable (default)" << endl;
    out << "  pru                   Pruning only, using EAP and cutoff based on the diagonal (no early abandoning!). Not for LCSS, SQED. " << endl;
    out << "  pru_la                PRU with last alignment cutoff adjustment. Not for LCSS, SQED" << endl;
    out << "  eap                   Early abandoned and pruned. SQED and LCSS are only early abandoned, not pruned" << endl;
    out << "  eap_la                EAP with last alignment cutoff adjustment. Not for LCSS." << endl;
    out << "Create an output file: '-out <outfile>'" << endl;
    out << "Examples:" << endl;
    out << "  " << execname << " -dist cdtw 0.2 lb-keogh2 -ucr /path/to/Univariate_ts -dataset Crop" << endl;
    out << endl;
}

/** Print an error followed by usage, then exit */
void print_error_exit(const string& invoc_string, const string& msg, int exit_code){
    std::cerr << "Error: " << msg << endl;
    print_usage(invoc_string, std::cerr);
    exit(exit_code);
}

/** Check if a path represent a readable file */
inline bool is_read_open(const fs::path& p) {
    return fs::exists(p) && ifstream(p);
}

enum RUNMODE{
    BASE,
    PRU,
    PRU_LA,
    EAP,
    EAP_LA,
};

string to_string(RUNMODE rm){
    switch(rm){
        case BASE: return "base";
        case PRU: return "pru";
        case PRU_LA: return "pru_la";
        case EAP: return "eap";
        case EAP_LA: return "eap_la";
    }
    return "";
}

enum DTW_LB{
    NONE,
    KEOGH,
    KEOGH2,
    WEBB,
};

string to_string(DTW_LB lb){
    switch(lb) {
        case NONE: return "lb-none";
        case KEOGH: return "lb-keogh";
        case KEOGH2: return "lb-keogh2";
        case WEBB: return "lb-webb";
    }
    return "";
}

typedef const TSData& TD;
#define TRAIN_TEST TD train, size_t idtrain, TD test, size_t idtest
typedef std::function<double(TRAIN_TEST, double ub)> distfun_t;
#define USE(s) s.data(), s.length()
#define TRAIN USE(train[idtrain])
#define TEST USE(test[idtest])

struct config {
    optional<fs::path> path_ucr{};
    optional<string> dataset{};
    optional<fs::path> path_train{};
    optional<fs::path> path_test{};
    optional<string> distance_name{};

    optional<fs::path> outpath{};

    DTW_LB dtw_lb{NONE};
    double cdtw_wratio{1}; // 1 by default: allow to use without change for dtw
    double wdtw_g{-1};
    double erp_gv{-1};
    double erp_wratio{-1};
    double lcss_epsilon{-1};
    double lcss_wratio{-1};
    bool lcss_int{false};
    double msm_cost{-1};
    double twe_nu{-1};
    double twe_lambda{-1};
    RUNMODE runmode{BASE};
};

/// Type of an envelope
typedef vector<double> env_t;

/** Helper used to compute the upper and lower envelopes of a dataset */
void compute_envelopes(TD dataset, vector<tuple<env_t, env_t>>& envelopes, size_t w){
    for(const auto& s:dataset.series){
        const auto l = s.length();
        vector<double> up(l);
        vector<double> lo(l);
        get_envelopes(s.data(), l, up.data(), lo.data(), w);
        envelopes.emplace_back(up, lo);
    }
}

/** Helper used to compute envelopes for lb Webb*/
void compute_envelopes_Webb(TD dataset, vector<tuple<env_t, env_t, env_t, env_t>>& envelopes, size_t w){
    for(const auto& s:dataset.series){
        const auto l = s.length();
        env_t up(l);
        env_t lo(l);
        env_t lo_up(l);
        env_t up_lo(l);
        get_envelopes_Webb(s.data(), l, up.data(), lo.data(), lo_up.data(), up_lo.data(), w);
        envelopes.emplace_back(up, lo, lo_up, up_lo);
    }
}

/** Only for DTW and CDTW, all series with same length.
 * Embed a distance under a lower bound
 */
variant<string, distfun_t> do_lb(distfun_t&& innerdist, const config& conf, size_t maxl, size_t minl, size_t w){
    // Pre-check
    if(conf.dtw_lb != NONE && minl != maxl){
        return {"Lower bound require same-length series"};
    }
    // Only for DTW and CDTW
    if(conf.distance_name=="dtw" || conf.distance_name=="cdtw"){
        switch(conf.dtw_lb){
            case NONE:{ return {innerdist}; }

            // --- --- ---
            case KEOGH:{
                // Captured state
                vector<tuple<vector<double>, vector<double>>> env_train;

                // Definition
                distfun_t dfun = [innerdist, env_train, w](TRAIN_TEST, double ub) mutable {
                    // Compute envelopes once
                    if(env_train.empty()){ compute_envelopes(train, env_train, w); }
                    // Lb Keogh
                    const auto& [up, lo] = env_train[idtrain];
                    double v = lb_Keogh(TEST, up.data(), lo.data(), ub);
                    if(v<ub){ v = innerdist(train, idtrain, test, idtest, ub); }
                    return v;
                };

                return {dfun};
            }

            // --- --- ---
            case KEOGH2:{
                // Captured state
                vector<tuple<vector<double>, vector<double>>> env_train;
                map<size_t, tuple<vector<double>, vector<double>>> env_test;

                // Definition
                distfun_t dfun = [innerdist, env_train, env_test, w](TRAIN_TEST, double ub) mutable {
                    // Compute envelopes once
                    if(env_train.empty()){ compute_envelopes(train, env_train, w); }
                    // Lb Keogh 1
                    const auto& [up, lo] = env_train[idtrain];
                    double v = lb_Keogh(TEST, up.data(), lo.data(), ub);
                    if(v<ub){
                        // LB Keogh 2
                        const auto lq = test[idtest].length();
                        if(!contains(env_test, idtest)){
                            vector<double> upq(lq);
                            vector<double> loq(lq);
                            get_envelopes(TEST, upq.data(), loq.data(), w);
                            env_test.emplace(idtest, std::tuple(std::move(upq), std::move(loq)));
                        }
                        const auto& [upq, loq] = env_test[idtest];
                        v = lb_Keogh(TRAIN, upq.data(), loq.data(), ub);
                        // Distance
                        if(v<ub){
                            v = innerdist(train, idtrain, test, idtest, ub);
                        }
                    }
                    return v;
                };

                return {dfun};
            }

            // --- --- ---
            case WEBB:{
                // Captured state
                vector<tuple<env_t, env_t, env_t, env_t>> env_train;
                map<size_t, tuple<env_t, env_t, env_t, env_t>> env_test;

                // Definition
                distfun_t dfun = [innerdist, env_train, env_test, w](TRAIN_TEST, double ub) mutable {
                    // Compute envelopes once for the train
                    if(env_train.empty()){ compute_envelopes_Webb(train, env_train, w); }
                    // Retrieve candidate
                    const auto& [up, lo, lo_up, up_lo] = env_train[idtrain];
                    // Compute envelopes of the query
                    const auto lq = test[idtest].length();
                    if(!contains(env_test, idtest)){
                        env_t q_up(lq);
                        env_t q_lo(lq);
                        env_t q_lo_up(lq);
                        env_t q_up_lo(lq);
                        get_envelopes_Webb(TEST, q_up.data(), q_lo.data(), q_lo_up.data(), q_up_lo.data(), w);
                        env_test.emplace(idtest, std::tuple(std::move(q_up), std::move(q_lo), std::move(q_lo_up), std::move(q_up_lo)));
                    }
                    const auto& [q_up, q_lo, q_lo_up, q_up_lo] = env_test[idtest];
                    //
                    double v = lb_Webb(
                            TEST, q_up.data(), q_lo.data(), q_lo_up.data(), q_up_lo.data(),
                            TRAIN, up.data(), lo.data(), lo_up.data(), up_lo.data(),
                            w, ub);
                    // Distance
                    if(v<ub){ v = innerdist(train, idtrain, test, idtest, ub); }
                    return v;
                };

                return dfun;
            }
        }
    } else {
        return {"Lower bound is incompatible with " + conf.distance_name.value()};
    }
    return {"Should not happen at " __FILE__ ": " STRINGIFY(__LINE__)};
}

/** Generate a distance function based on the configuration and the maximum length of the series.
 * The max length is used to compute the actual window sized, which is specified by a ratio */
variant<string, distfun_t> get_distfun(const config& conf, size_t maxl, size_t minl){
    distfun_t dfun;
    if(conf.distance_name == "dtw"){
        switch(conf.runmode){
            case BASE:{
                dfun = [](TRAIN_TEST, [[maybe_unused]]double ub){ return dtw_base(TRAIN, TEST); };
                break;
            }
            case PRU:{
                dfun = [](TRAIN_TEST, [[maybe_unused]]double ub){ return dtw<false>(TRAIN, TEST); };
                break;
            }
            case PRU_LA:{
                dfun = [](TRAIN_TEST, [[maybe_unused]]double ub){ return dtw<true>(TRAIN, TEST); };
                break;
            }
            case EAP:{
                dfun = [](TRAIN_TEST, double ub){ return dtw<false>(TRAIN, TEST, ub); };
                break;
            }
            case EAP_LA:{
                dfun = [](TRAIN_TEST, double ub){ return dtw<true>(TRAIN, TEST, ub); };
                break;
            }
            default:{ return {conf.distance_name.value() + " distance not compatible with " +  to_string(conf.runmode) }; }
        }
        size_t w = conf.cdtw_wratio * maxl;
        return do_lb(std::move(dfun), conf, maxl, minl, w);
    }
    else if (conf.distance_name == "cdtw"){
        size_t w = conf.cdtw_wratio * maxl;
        switch(conf.runmode){
            case BASE:{
                dfun = [w](TRAIN_TEST, [[maybe_unused]]double ub){ return cdtw_base(TRAIN, TEST, w); };
                break;
            }
            case PRU:{
                dfun = [w](TRAIN_TEST, [[maybe_unused]]double ub){ return cdtw<false>(TRAIN, TEST, w); };
                break;
            }
            case PRU_LA:{
                dfun = [w](TRAIN_TEST, [[maybe_unused]]double ub){ return cdtw<true>(TRAIN, TEST, w); };
                break;
            }
            case EAP:{
                dfun = [w](TRAIN_TEST, double ub){ return cdtw<false>(TRAIN, TEST, w, ub); };
                break;
            }
            case EAP_LA:{
                dfun = [w](TRAIN_TEST, double ub){ return cdtw<true>(TRAIN, TEST, w, ub); };
                break;
            }
            default:{ return {conf.distance_name.value() + " distance not compatible with " +  to_string(conf.runmode) }; }
        }
        return do_lb(std::move(dfun), conf, maxl, minl, w);
    }
    else if (conf.distance_name == "wdtw"){
        auto weights = generate_weights(conf.wdtw_g, maxl);
        switch(conf.runmode){
            case BASE:{
                dfun = [weights](TRAIN_TEST, [[maybe_unused]]double ub){ return wdtw_base(TRAIN, TEST, weights.data()); };
                break;
            }
            case PRU:{
                dfun = [weights](TRAIN_TEST, [[maybe_unused]]double ub){ return wdtw<false>(TRAIN, TEST, weights.data()); };
                break;
            }
            case PRU_LA:{
                dfun = [weights](TRAIN_TEST, [[maybe_unused]]double ub){ return wdtw<true>(TRAIN, TEST, weights.data()); };
                break;
            }
            case EAP:{
                dfun = [weights](TRAIN_TEST, double ub){ return wdtw<false>(TRAIN, TEST, weights.data(), ub); };
                break;
            }
            case EAP_LA:{
                dfun = [weights](TRAIN_TEST, double ub){ return wdtw<true>(TRAIN, TEST, weights.data(), ub); };
                break;
            }
            default:{ return {conf.distance_name.value() + " distance not compatible with " +  to_string(conf.runmode) }; }
        }
    }
    else if (conf.distance_name == "sqed"){
        switch(conf.runmode){
            case BASE:{ dfun = [](TRAIN_TEST, [[maybe_unused]]double ub){ return elementwise(TRAIN, TEST); }; break; }
            case EAP:{ dfun = [](TRAIN_TEST, double ub){ return elementwise<false>(TRAIN, TEST, ub); }; break; }
            case EAP_LA:{ dfun = [](TRAIN_TEST, double ub){ return elementwise<true>(TRAIN, TEST, ub); }; break; }
            default:{ return {conf.distance_name.value() + " distance not compatible with " +  to_string(conf.runmode) }; }
        }
    }
    else if (conf.distance_name == "erp"){
        auto gv = conf.erp_gv;
        auto w = conf.erp_gv * maxl;
        switch(conf.runmode){
            case BASE:{
                dfun = [gv, w](TRAIN_TEST, [[maybe_unused]]double ub){ return erp_base(TRAIN, TEST, gv, w); };
                break;
            }
            case PRU:{
                dfun = [gv, w](TRAIN_TEST, [[maybe_unused]]double ub){ return erp<false>(TRAIN, TEST, gv, w); };
                break;
            }
            case PRU_LA:{
                dfun = [gv, w](TRAIN_TEST, [[maybe_unused]]double ub){ return erp<true>(TRAIN, TEST, gv, w); };
                break;
            }
            case EAP:{
                dfun = [gv, w](TRAIN_TEST, double ub){ return erp<false>(TRAIN, TEST, gv, w, ub); };
                break;
            }
            case EAP_LA:{
                dfun = [gv, w](TRAIN_TEST, double ub){ return erp<true>(TRAIN, TEST, gv, w, ub); };
                break;
            }
            default:{ return {conf.distance_name.value() + " distance not compatible with " +  to_string(conf.runmode) }; }
        }
    }
    else if (conf.distance_name == "lcss"){
        double e = conf.lcss_epsilon;
        double w = conf.lcss_wratio * maxl;
        if(conf.lcss_int){ w = conf.lcss_wratio; } // Hack for LCSS EE parameter
        switch(conf.runmode) {
            case BASE:{
                dfun = [e, w](TRAIN_TEST, [[maybe_unused]]double ub){ return lcss(TRAIN, TEST, e, w); };
                break;
            }
            case EAP:{
                dfun = [e, w](TRAIN_TEST, double ub){ return lcss(TRAIN, TEST, e, w, ub); };
                break;
            }
            default:{ return {conf.distance_name.value() + " distance not compatible with " +  to_string(conf.runmode) }; }
        }
    }
    else if (conf.distance_name == "msm"){
        double c = conf.msm_cost;
        switch(conf.runmode){
            case BASE:{
                dfun = [c](TRAIN_TEST, [[maybe_unused]]double ub){ return msm_base(TRAIN, TEST, c); };
                break;
            }
            case PRU:{
                dfun = [c](TRAIN_TEST, [[maybe_unused]]double ub){ return msm<false>(TRAIN, TEST, c); };
                break;
            }
            case PRU_LA:{
                dfun = [c](TRAIN_TEST, [[maybe_unused]]double ub){ return msm<true>(TRAIN, TEST, c); };
                break;
            }
            case EAP:{
                dfun = [c](TRAIN_TEST, double ub){ return msm<false>(TRAIN, TEST, c, ub); };
                break;
            }
            case EAP_LA:{
                dfun = [c](TRAIN_TEST, double ub){ return msm<true>(TRAIN, TEST, c, ub); };
                break;
            }
            default:{ return {conf.distance_name.value() + " distance not compatible with " +  to_string(conf.runmode) }; }
        }
    }
    else if (conf.distance_name == "twe"){
        double n = conf.twe_nu;
        double l = conf.twe_lambda;
        switch(conf.runmode){
            case BASE:{
                dfun = [n, l](TRAIN_TEST, [[maybe_unused]]double ub){ return twe_base(TRAIN, TEST, n, l); };
                break;
            }
            case PRU:{
                dfun = [n, l](TRAIN_TEST, [[maybe_unused]]double ub){ return twe<false>(TRAIN, TEST, n, l); };
                break;
            }
            case PRU_LA:{
                dfun = [n, l](TRAIN_TEST, [[maybe_unused]]double ub){ return twe<true>(TRAIN, TEST, n, l); };
                break;
            }
            case EAP:{
                dfun = [n, l](TRAIN_TEST, double ub){ return twe<false>(TRAIN, TEST, n, l, ub); };
                break;
            }
            case EAP_LA:{
                dfun = [n, l](TRAIN_TEST, double ub){ return twe<true>(TRAIN, TEST, n, l, ub); };
                break;
            }
            default:{ return {conf.distance_name.value() + " distance not compatible with " +  to_string(conf.runmode) }; }
        }
    }

    return dfun;
}

/** NN1, in case of ties, first found win
 * Return a tuple (nb correct, accuracy, duration)
 * where accuracy = nb correct/test size*/
variant<string, tuple<size_t, double, duration_t>> NN1(const config& conf, const TSData& train, const TSData& test){
    // --- --- --- Get the distance function
    auto maxl = max(train.longest_length, test.longest_length);
    auto minl = min(train.longest_length, test.longest_length);
    variant<string, distfun_t> v_dfun = get_distfun(conf, maxl, minl);
    distfun_t dfun;
    if(v_dfun.index() == 0){ return {get<0>(v_dfun)}; }
    dfun = get<1>(v_dfun);

    // --- --- --- NN1 loop
    double nb_correct{0};
    duration_t duration{0};
    for(size_t q_idx=0; q_idx < test.series.size(); q_idx++){
        double bsf = POSITIVE_INFINITY;
        const TSeries* bcandidates = nullptr;
        for(size_t c_idx=0; c_idx < train.series.size(); c_idx++){
            // Time the call
            auto start = now();
            double res = dfun(train, c_idx, test, q_idx, bsf);
            auto stop = now();
            duration += (stop-start);
            // update BSF
            if(res<bsf){
                bsf = res;
                bcandidates = &train[c_idx];
            }
        }
        if(bcandidates!= nullptr && bcandidates->label().value() == test[q_idx].label().value()){
            nb_correct++;
        }
    }

    return {tuple<size_t, double, duration_t>{nb_correct, nb_correct/(test.series.size()), duration}};
}

/** Config to JSON */
string conf_to_JSON(const config& conf){
    stringstream ss;
    ss << R"({"name":")" << conf.distance_name.value() <<'"';
    ss << R"(, "runmode":)";
    ss << '"' + to_string(conf.runmode) + '"';

    if(conf.distance_name == "dtw"){
        ss << R"(, "lb":)" << '"'+to_string(conf.dtw_lb)+'"';
    }
    else if (conf.distance_name == "cdtw"){
        ss << R"(, "wratio":)" << conf.cdtw_wratio;
        ss << R"(, "lb":)" << '"'+to_string(conf.dtw_lb)+'"';
    }
    else if (conf.distance_name == "wdtw"){
        ss << R"(, "g":)" << conf.wdtw_g;
    }
    else if (conf.distance_name == "erp"){
        ss << R"(, "gv":)" << conf.erp_gv;
        ss << R"(, "wratio":)" << conf.erp_wratio;
    }
    else if (conf.distance_name == "lcss"){
        ss << R"(, "epsilon":)" << conf.lcss_epsilon;
        ss << R"(, "wratio":)" << conf.lcss_wratio;
    }
    else if (conf.distance_name == "msm"){
        ss << R"(, "cost":)" << conf.msm_cost;
    }
    else if (conf.distance_name == "twe"){
        ss << R"(, "nu":)" << conf.twe_nu;
        ss << R"(, "lambda":)" << conf.twe_lambda;
    }
    ss << "}";
    return ss.str();
}

/** Dataset to JSON */
string ds_to_JSON(const TSData& tsd){
    stringstream ss;
    ss << R"({"size":)" << tsd.series.size();
    ss << R"(, "minl":)" << tsd.shortest_length;
    ss << R"(, "maxl":)" << tsd.longest_length;
    ss << R"(, "missing":)" << (tsd.missing.value()?"true":"false");
    ss << R"(, "nbclasses":)" << (tsd.labels.size());
    ss << "}";
    return ss.str();
}

int main(int argc, char** argv){
    config conf;

    // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
    // Command line parsing
    // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
    {
        // No argument: print the usage
        if (argc == 1) { print_usage(argv[0], cout); exit(0); }

        // Else, read the args. Declare local var.
        optional<string> error;
        int i{1};
        string arg;

        // Lambda: arg reading
        auto next_arg = [argv, &i]() { string res{argv[i]}; ++i; return res; };

        // Argument parsing loop
        while (i < argc) {
            arg = next_arg();

            // We work with the full dataset directory
            if (arg == "-ucr") {
                if (i < argc) {
                    fs::path path(next_arg());
                    if (fs::exists(path) && fs::is_directory(path)) {
                        conf.path_ucr = fs::canonical(path);
                    } else { error = {"Cannot find directory " + string(path)}; }
                } else { error = {"Directory path expected after '-ucr'"}; }
            }

                // We can limit to a given dataset
            else if(arg == "-dataset"){
                if (i < argc) {
                    conf.dataset = next_arg();
                } else { error = {"Dataset name expected after '-dataset'"}; }
            }

                // Specify output file
            else if(arg == "-out"){
                if(i<argc){
                    arg = next_arg();
                    fs::path path(arg);
                    try {
                        if(path.has_parent_path()) { fs::create_directories(path.parent_path()); }
                        ofstream of(path);
                        of.close();
                        conf.outpath = {path};
                    }  catch (std::exception& e) {
                        error = {e.what()};
                    }
                } else { error = {"Directory path expected after '-out'"}; }
            }

                // Specify the distance to run
            else if(arg == "-dist"){
                if(i<argc){
                    auto distname = next_arg();
                    std::transform(distname.begin(), distname.end(), distname.begin(),
                                   [](unsigned char c){ return std::tolower(c); });
                    conf.distance_name = {distname};
                    // --- --- ---
                    if(distname == "dtw") {
                        /* Maybe a lower bound */
                        if(i<argc){
                            arg=next_arg();
                            if(arg=="lb-none"){conf.dtw_lb=NONE;}
                            else if( arg == "lb-keogh"){conf.dtw_lb=KEOGH;}
                            else if( arg == "lb-keogh2"){conf.dtw_lb=KEOGH2;}
                            else if( arg == "lb-webb"){conf.dtw_lb=WEBB;}
                            else{--i; /* revert */}
                        }
                    } // --- --- ---
                    else if(distname == "cdtw") {
                        if(i<argc){
                            arg = next_arg();
                            auto wr = as_double(arg);
                            if(wr && 0<= wr.value() && wr.value() <= 1){ conf.cdtw_wratio = wr.value(); }
                            else { error = {"cdtw expects a windows value 0<=w<=1"}; }
                        } else { error = {"cdtw expects a windows value 0<=w<=1"}; }
                        /* Maybe a lower bound */
                        if(i<argc){
                            arg=next_arg();
                            if(arg=="lb-none"){conf.dtw_lb=NONE;}
                            else if( arg == "lb-keogh"){conf.dtw_lb=KEOGH;}
                            else if( arg == "lb-keogh2"){conf.dtw_lb=KEOGH2;}
                            else if( arg == "lb-webb"){conf.dtw_lb=WEBB;}
                            else{--i; /* revert */}
                        }
                    } // --- --- ---
                    else if (distname == "wdtw") {
                        if(i<argc){
                            arg = next_arg();
                            auto g = as_double(arg);
                            if(g && 0<= g.value()) { conf.wdtw_g = g.value(); }
                            else {error={"wdtw expects a weight factor 0<=g"}; }
                        } else {error={"wdtw expects a weight factor 0<=g"}; }
                    } // --- --- --
                    else if (distname == "sqed"){
                        /* no parameter */
                    } // --- --- ---
                    else if(distname == "erp"){
                        if(i<argc){
                            arg = next_arg();
                            auto gv = as_double(arg);
                            arg = next_arg();
                            auto wr = as_double(arg);
                            if(gv && 0<=gv.value() && wr && 0<=wr.value() && wr.value()<=1){
                                conf.erp_gv = gv.value();
                                conf.erp_wratio = wr.value();
                            }
                            else {error={"erp expects a gap value 0<=gv followed by a window ratio 0<=w<=1"};}
                        } else {error={"erp expects a gap value 0<=gv followed by a window ratio 0<=w<=1"};}
                    } // --- --- ---
                    else if(distname == "lcss"){
                        if(i+1<argc){
                            arg = next_arg();
                            auto epsilon = as_double(arg);
                            typeof(epsilon) wr;
                            arg = next_arg();
                            if(arg == "int"){ // Hack to support integer window from EE
                                arg=next_arg();
                                conf.lcss_int = true;
                            }
                            wr = as_double(arg);
                            if(epsilon && 0<=epsilon.value() && wr && ((!conf.lcss_int && 0<=wr.value() && wr.value()<=1)||(0<=wr.value()))){
                                conf.lcss_epsilon = epsilon.value();
                                conf.lcss_wratio = wr.value();
                            }
                            else {error= {"lcss expects an epsilon value followed by a window ratio 0<=w<=1"}; }
                        } else {error= {"lcss expects an epsilon value followed by a window ratio 0<=w<=1"}; }
                    } // --- --- ---
                    else if(distname == "msm"){
                        if(i<argc){
                            arg = next_arg();
                            auto cost = as_double(arg);
                            if(cost && 0<=cost.value()){
                                conf.msm_cost = cost.value();
                            } else { error={"msm expects a cost 0<=c"}; }
                        } else { error={"msm expects a cost 0<=c"}; }
                    } // --- --- ---
                    else if(distname == "twe"){
                        if(i+1<argc) {
                            arg = next_arg();
                            auto nu = as_double(arg);
                            arg = next_arg();
                            auto lambda = as_double(arg);
                            if(nu && 0<= nu.value() && lambda && 0<=lambda.value()){
                                conf.twe_nu = nu.value();
                                conf.twe_lambda = lambda.value();
                            }
                        } else { error={"twe expects a stiffness parameter 0<=nu followed by a cost parameter 0<=lambda"}; }
                    } else { error={"twe expects a stiffness parameter 0<=nu followed by a cost parameter 0<=lambda"}; }
                } else {error = {"Distance name expected after '-dist'"}; }
            }

                // --- --- --- Checking how to run it
            else if(arg=="base"){ conf.runmode = BASE; }
            else if(arg=="pru"){ conf.runmode = PRU; }
            else if(arg=="pru_la"){ conf.runmode = PRU_LA; }
            else if(arg=="eap"){ conf.runmode = EAP; }
            else if(arg=="eap_la"){ conf.runmode = EAP_LA; }
                // --- --- --- Unkwnon args
            else { error = {"Unkwnon arg: " + arg}; }

            // --- --- --- Error checking
            if (error.has_value()) { print_error_exit(argv[0], error.value(), 1); }
        }

        // Argument post check
        if(!conf.distance_name){
            print_error_exit(argv[0], string("No distance specified"), 1);
        }
        if(!conf.path_ucr){
            print_error_exit(argv[0], string("No UCR path specified"), 1);
        }
        if(conf.dataset.has_value()){
            auto p = conf.path_ucr.value() / conf.dataset.value();
            if(!fs::is_directory(p)){ print_error_exit(argv[0], string(p)+" is not a directory", 1); }
            auto ptrain = p / (conf.dataset.value() + "_TRAIN.ts");
            auto ptest = p / (conf.dataset.value() + "_TEST.ts");
            if(!is_read_open(ptrain)){
                print_error_exit(argv[0], "Cannot read " + string(ptrain), 1);
            }
            if(!is_read_open(ptest)){
                print_error_exit(argv[0], "Cannot read " + string(ptest), 1);
            }
            conf.path_train = {ptrain};
            conf.path_test = {ptest};
        }
    }

    // Parameter recall
    {
        cout << "Parameters = ";
        cout << " -ucr " << conf.path_ucr.value();
        if(conf.dataset) { cout << " -dataset " << conf.dataset.value(); }
        cout << " -dist " << conf.distance_name.value() << " ";
        if(conf.distance_name == "dtw"){ cout << to_string(conf.dtw_lb); }
        else if(conf.distance_name == "cdtw"){ cout << conf.cdtw_wratio << " " << to_string(conf.dtw_lb); }
        else if(conf.distance_name == "wdtw"){ cout << conf.wdtw_g; }
        else if(conf.distance_name == "sqed"){ /* */ }
        else if(conf.distance_name == "erp"){ cout << conf.erp_gv << " " << conf.erp_wratio; }
        else if(conf.distance_name == "lcss"){ cout << conf.lcss_epsilon << " " << conf.lcss_wratio; }
        else if(conf.distance_name == "msm"){ cout << conf.msm_cost; }
        else if(conf.distance_name == "twe"){ cout << conf.twe_nu << " " << conf.twe_lambda; }
        cout << " " << to_string(conf.runmode) << endl;
    }

    // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
    // Collect dataset (paths to train and test files) to launch
    // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
    vector<tuple<string, fs::path, fs::path>> train_test;
    if(conf.dataset){
        train_test.emplace_back(conf.dataset.value(), conf.path_train.value(), conf.path_test.value());
    } else {
        cout << "Looking for UCR dataset in ts format in " << conf.path_ucr.value() << ":" << endl;
        for(const auto& entry: fs::directory_iterator(conf.path_ucr.value())){
            if(entry.is_directory()){
                const auto dataset_name = entry.path().filename().string();
                const auto train_p = entry.path() / (dataset_name+"_TRAIN.ts");
                const auto test_p = entry.path() / (dataset_name+"_TEST.ts");
                if(fs::exists(train_p) && fs::is_regular_file(train_p) && fs::exists(test_p) && fs::is_regular_file(test_p)){
                    cout << "  Adding " << entry << endl;
                    train_test.emplace_back(dataset_name, train_p, test_p);
                } else {
                    cout << "  Skipping " << entry << " (could not find "+dataset_name+"_TRAIN.ts/_TEST.ts file)" << endl;
                }
            } else {
                cout << "  Skipping " << entry << " (not a directory)" << endl;
            }
        }
    }

    // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
    // Launch per dataset
    // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
    vector<string> vstr;
    for(const auto& entry :train_test){
        const auto& [name, trainp, testp] = entry;
        cout << "Loading data for '" << name << "'... ";
        flush(cout);
        auto itrain = ifstream(trainp);
        auto itest = ifstream(testp);
        auto start = now();
        auto res_train = TSReader::read(itrain);
        auto res_test = TSReader::read(itest);
        auto stop = now();
        if(res_train.index() == 0){ print_error_exit(argv[0], get<0>(res_train), 2); }
        if(res_test.index() == 0){ print_error_exit(argv[0], get<0>(res_test), 2); }
        cout << "Done in "; printDuration(cout, stop-start); cout << endl;
        const auto& train = get<1>(res_train);
        const auto& test = get<1>(res_test);
        cout << "TRAIN: " << train.series.size() << " x ";
        if(train.has_equallength()){ cout << train.shortest_length; }
        else { cout << "[" << train.shortest_length << ".." << train.longest_length << "]"; }
        cout << endl;
        cout << "TEST: " << test.series.size() << " x ";
        if(test.has_equallength()){ cout << train.shortest_length; }
        else { cout << "[" << test.shortest_length << ".." << train.longest_length << "]"; }
        cout << endl;

        // --- --- --- NN1
        auto res = NN1(conf, train, test);
        switch (res.index()) {
            case 0: {
                print_error_exit(argv[0], get<0>(res), 3);
            }
            case 1: {
                // Hack
                if(conf.lcss_int){
                    double m = std::max<double>(train.longest_length, test.longest_length);
                    conf.lcss_wratio = conf.lcss_wratio/m;
                }
                auto[nbcorrect, acc, duration] = get<1>(res);
                stringstream ss;
                ss << "{" << endl;
                ss << R"(  "type":"NN1",)" << endl;
                ss << R"(  "nb_correct":")" << nbcorrect << "\"," << endl;
                ss << R"(  "accuracy":")" << acc << "\"," << endl;
                ss << R"(  "dataset":")" << name << "\"," << endl;
                ss << R"(  "distance":)" << conf_to_JSON(conf) << ',' << endl;
                ss << R"(  "train":)" << ds_to_JSON(train) << "," << endl;
                ss << R"(  "test":)" << ds_to_JSON(test) << "," << endl;
                ss << R"(  "timing_ns":)" << duration.count() << "," << endl;
                ss << R"(  "timing":")"; printDuration(ss, duration); ss << "\"" << endl;
                ss << "}" << endl;
                string str = ss.str();
                std::cout << str;
                vstr.push_back(str);
            }
        }

        if(conf.outpath) {
            ofstream of(conf.outpath.value());
            if(vstr.size()>1) {
                of << "[" << std::endl;
                auto i = vstr.size();
                for (const auto &str:vstr) {
                    of << str;
                    if (i > 1) {
                        --i;
                        of << ',' << std::endl;
                    }
                }
                of << "]";
            } else {
                of<<vstr[0];
            }
            of.close();
        }

    }
}
