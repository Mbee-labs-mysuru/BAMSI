#pragma once

#include "bamsi/status.hpp"

#include <stdexcept>
#include <utility>

namespace bamsi {

template <typename T>
class Result {
public:
    Result(const T& value)
        : status_(Status::Ok()), value_(value), has_value_(true) {}

    Result(T&& value)
        : status_(Status::Ok()), value_(std::move(value)), has_value_(true) {}

    Result(const Status& status)
        : status_(status), value_(), has_value_(false) {}

    bool ok() const { return status_.ok(); }
    const Status& status() const { return status_; }
    bool has_value() const { return has_value_; }

    const T& value() const {
        if (!has_value_) {
            throw std::logic_error("Result has no value");
        }
        return value_;
    }

    T& value() {
        if (!has_value_) {
            throw std::logic_error("Result has no value");
        }
        return value_;
    }

private:
    Status status_;
    T value_{};
    bool has_value_{false};
};

}  // namespace bamsi
