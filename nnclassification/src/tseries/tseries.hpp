#ifndef EXP_TSERIES_HPP
#define EXP_TSERIES_HPP

#include "../utils.hpp"

class TSeries : private Uncopyable {
private:
    // --- --- --- --- --- --- --- --- --- --- -- --- --- --- --- -- --- --- --- --- -- --- --- --- --- -- --- --- ---
    // Only when owning data
    // --- --- --- --- --- --- --- --- --- --- -- --- --- --- --- -- --- --- --- --- -- --- --- --- --- -- --- --- ---
    std::vector<double> data_v_{};

    // --- --- --- --- --- --- --- --- --- --- -- --- --- --- --- -- --- --- --- --- -- --- --- --- --- -- --- --- ---
    // Common internal interface
    // --- --- --- --- --- --- --- --- --- --- -- --- --- --- --- -- --- --- --- --- -- --- --- --- --- -- --- --- ---
    const double *data_{nullptr};
    size_t length_{0};
    size_t nb_dimensions_{0};
    bool has_missing_{false};
    std::optional<std::string> label_{};

public:

    // --- --- --- --- --- --- --- --- --- --- -- --- --- --- --- -- --- --- --- --- -- --- --- --- --- -- --- --- ---
    // Constructors/Destructors/Copy (uncopyable!)/Movement
    // --- --- --- --- --- --- --- --- --- --- -- --- --- --- --- -- --- --- --- --- -- --- --- --- --- -- --- --- ---

    /** Default constructor, creating an empty series */
    TSeries() : data_(data_v_.data()) {}


    /*
    template<typename FloatType>
    static void derivative(const FloatType *series, size_t length, FloatType *out) {
        if (length > 2) {
            for (size_t i{1}; i < length - 1; ++i) {
                out[i] = ((series[i] - series[i - 1]) + ((series[i + 1] - series[i - 1]) / 2.0)) / 2.0;
            }
            out[0] = out[1];
            out[length - 1] = out[length - 2];
        } else {
            std::copy(series, series + length, out);
        }
    }*/

    /** Constructor with a moving vector: the new instance owns the data and will clear them when destroyed */
    TSeries(std::vector<double> &&data, size_t nb_dimensions, bool has_missing, std::optional<std::string> label)
            : data_v_(std::move(data)), nb_dimensions_(nb_dimensions), has_missing_(has_missing),
              label_(std::move(label)) {
        assert((data_v_.size() % nb_dimensions_) == 0);
        length_ = data_v_.size() / nb_dimensions_;
        data_v_.shrink_to_fit();
        data_ = data_v_.data();
        // std::vector<double> nd(data_v_.size());
        // derivative<double>(data_v_.data(), data_v_.size(), nd.data());
        // data_ = nd.data();
        // data_v_ = std::move(nd);
    }

    /** Constructor with a raw pointer: the new instance does not own the data and will not free the memory when destroyed.
     * To be use with FFI, e.g. Python */
    TSeries(const double *data_ptr, size_t length, size_t nb_dimensions, bool has_missing,
            std::optional<std::string> label)
            : data_(data_ptr), length_(length), nb_dimensions_(nb_dimensions), has_missing_(has_missing),
              label_(std::move(label)) {
        assert(length == 0 || data_ptr != nullptr);
        assert((length_ % nb_dimensions_) == 0);
    }

    /** Constructor from an other time series: the new instance does not own the data.
     * Be sure that the backing instance lives longer than the new one. */
    TSeries(TSeries &backing) :
            Uncopyable(),
            data_(backing.data_), length_(backing.length_), nb_dimensions_(backing.nb_dimensions_),
            has_missing_(backing.has_missing_), label_(backing.label_) {}

    /** Construct a new series owning its data, using the information of another time series. */
    TSeries(std::vector<double> &&data, const TSeries &info_source) :
            data_v_(std::move(data)),
            length_(info_source.length_), nb_dimensions_(info_source.nb_dimensions_),
            has_missing_(info_source.has_missing_), label_(info_source.label_) {
        data_v_.shrink_to_fit();
        data_ = data_v_.data();
    }


    // Destructor
    ~TSeries() = default;


    // Movement
    TSeries(TSeries &&other) noexcept {
        // Use the move-assignment operator
        *this = std::move(other);
    }

    // Move-assignment
    TSeries &operator=(TSeries &&other) noexcept {
        if (this != &other) {
            if (other.is_owning()) {
                data_v_ = std::move(other.data_v_);
                data_ = data_v_.data();
            } else {
                data_ = other.data_;
            }
            nb_dimensions_ = other.nb_dimensions_;
            length_ = other.length_;
            has_missing_ = other.has_missing_;
            label_ = std::move(other.label_);
        }
        return *this;
    }



    // --- --- --- --- --- --- --- --- --- --- -- --- --- --- --- -- --- --- --- --- -- --- --- --- --- -- --- --- ---
    // Methods
    // --- --- --- --- --- --- --- --- --- --- -- --- --- --- --- -- --- --- --- --- -- --- --- --- --- -- --- --- ---

    // Basic accessors:
    [[nodiscard]] inline size_t length() const { return length_; }

    [[nodiscard]] inline size_t nb_dimensions() const { return nb_dimensions_; }

    [[nodiscard]] inline bool has_missing() const { return has_missing_; }

    [[nodiscard]] inline const std::optional<std::string> &label() const { return label_; }

    [[nodiscard]] inline const double *data() const { return data_; }

    /** Access to the start of a dimension (first dimension is 0) */
    [[nodiscard]] inline const double *operator[](size_t dim) const { return data_ + (nb_dimensions_ * dim); }

    /** Access a value using a pair of coordinate (Dimension,index) */
    [[nodiscard]] inline double operator()(size_t dim, size_t idx) const {
        return *(data_ + (nb_dimensions_ * dim + idx));
    }

    /** return true if the series is owning its data */
    [[nodiscard]] inline bool is_owning() { return data_v_.data() == data_; }

};

inline bool operator==(const TSeries &lhs, const TSeries &rhs) {
    bool res1 = lhs.nb_dimensions() == rhs.nb_dimensions()
                && lhs.length() == rhs.length()
                && lhs.has_missing() == rhs.has_missing()
                && lhs.label() == rhs.label();
    if (res1) {
        const auto *ld = lhs.data();
        const auto *rd = rhs.data();
        if (ld == rd) {
            return true; // Same pointer, so all good
        } else { // Else, compare item one by one
            bool same = true;
            size_t index = 0;
            while (same && index < lhs.length()) {
                same = ld[index] == rd[index];
                ++index;
            }
            return same;
        }
    } else {
        return false;
    }
}


inline bool operator!=(const TSeries &lhs, const TSeries &rhs) { return !operator==(lhs, rhs); }


using TSTransform = TSeries(*)(const TSeries &);



#endif