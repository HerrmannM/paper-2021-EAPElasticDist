#ifndef EXP_TSREADER_HPP
#define EXP_TSREADER_HPP

#include <functional>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <sstream>
#include <variant>
#include <vector>

#include "../../tseries.hpp"


/** Structure for the TS file format */
class TSData : private Uncopyable {
public:

    // --- --- --- --- --- --- --- --- --- --- -- --- --- --- --- -- --- --- --- --- -- --- --- --- --- -- --- --- ---
    // Fields:
    // --- --- --- --- --- --- --- --- --- --- -- --- --- --- --- -- --- --- --- --- -- --- --- --- --- -- --- --- ---

    // --- --- --- Header
    std::optional<std::string> problem_name{};
    std::optional<bool> timestamps{};
    std::optional<bool> missing{};
    std::optional<bool> univariate{};
    std::optional<bool> equallength{};
    std::optional<size_t> serieslength{};
    std::set<std::string> labels{};

    // --- --- --- Data
    std::vector<TSeries> series;

    // --- --- --- Extra
    size_t nb_dimensions{0};
    size_t shortest_length{std::numeric_limits<size_t>::max()};
    size_t longest_length{0};
    std::vector<size_t> series_with_missing_values{};

    // --- --- --- --- --- --- --- --- --- --- -- --- --- --- --- -- --- --- --- --- -- --- --- --- --- -- --- --- ---
    // Constructors/Destructors/Copy (uncopyable!)/Movement
    // --- --- --- --- --- --- --- --- --- --- -- --- --- --- --- -- --- --- --- --- -- --- --- --- --- -- --- --- ---

    TSData() = default;

    ~TSData() = default;

    TSData(TSData &&other) = default;

    TSData &operator=(TSData &&other) = default;

    // --- --- --- --- --- --- --- --- --- --- -- --- --- --- --- -- --- --- --- --- -- --- --- --- --- -- --- --- ---
    // Helpers
    // --- --- --- --- --- --- --- --- --- --- -- --- --- --- --- -- --- --- --- --- -- --- --- --- --- -- --- --- ---

    [[nodiscard]] inline const TSeries& operator[](std::vector<double>::size_type index) const { return series[index]; }

    [[nodiscard]] inline bool use_timestemps() const { return timestamps.value_or(false); }

    [[nodiscard]] inline bool is_univariate() const { return univariate.value_or(false); }

    [[nodiscard]] inline bool has_missings() const { return missing.value_or(false); }

    [[nodiscard]] inline bool has_equallength() const { return equallength.value_or(false); }

    [[nodiscard]] inline bool has_labels() const { return !labels.empty(); }
};


/** Allow to read an input stream into a TSData structure. */
class TSReader {
public:

    /** Read an input stream into a TSData.
     *  Returns a variant containing:
     *  * either an error-string if an error occured
     *  * or the data on success
     */
    static std::variant<std::string, TSData> read(std::istream &input);

private:
    // --- --- --- Private constructor
    explicit TSReader(std::istream &input) : input(input), state(nullptr) {}

    // --- --- --- Alias

    // The presence of a string signify an error
    using Result = std::optional<std::string>;
    // A state is a member function returning a result
    using State = Result (TSReader::*)();

    // --- --- --- Internal state
    std::istream &input;
    std::string buffer;
    State state;
    TSData data;
    // Length of the first series
    size_t length1st{};

    // --- --- --- Header's directives parsing tools
    // Directive witch code: used to simulate a switch on string through a map
    enum class DirectiveCode {
        dir_problem_name,
        dir_timestamp,
        dir_missing,
        dir_univariate,
        dir_equal_length,
        dir_series_length,
        dir_class_label,
        dir_data,
    };

    // Constants in lower case.
    inline static std::string str_problem_name = "problemname";
    inline static std::string str_timestamp = "timestamps";
    inline static std::string str_missing = "missing";
    inline static std::string str_univariate = "univariate";
    inline static std::string str_equallength = "equallength";
    inline static std::string str_serieslength = "serieslength";
    inline static std::string str_classlabel = "classlabel";
    inline static std::string str_data = "data";

    // Map (directive string |-> directive switch code)
    inline static std::map<std::string, DirectiveCode> directive_map = {
            {str_problem_name, DirectiveCode::dir_problem_name},
            {str_timestamp,    DirectiveCode::dir_timestamp},
            {str_missing,      DirectiveCode::dir_missing},
            {str_univariate,   DirectiveCode::dir_univariate},
            {str_equallength,  DirectiveCode::dir_equal_length},
            {str_serieslength, DirectiveCode::dir_series_length},
            {str_classlabel,   DirectiveCode::dir_class_label},
            {str_data,         DirectiveCode::dir_data}
    };

    // --- --- --- Methods

    // Main reading loop
    std::variant<std::string, TSData> read();

    // Loop over the header, skipping whites and line comments.
    // Change the state to read_directive when reading a '@'
    Result read_header();

    // Read the word w following a '@' (coming from read_header)
    // If w is a directive, read its argument and change the state back to read_header.
    // Else, it is an error.
    Result read_directive();

    // Read classes: space separated integer.
    Result read_classes();

    // Read the data
    Result read_data();

    // Working method for read_data
    Result read_data_(std::istream &in);
};

#endif
