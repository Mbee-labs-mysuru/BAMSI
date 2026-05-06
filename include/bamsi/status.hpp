#pragma once

#include <string>
#include <utility>

namespace bamsi {

enum class StatusCode {
    kOk = 0,
    kInvalidArgument,
    kNotFound,
    kInternalError,
    kNotImplemented
};

class Status {
public:
    Status() : code_(StatusCode::kOk), message_("OK") {}
    Status(StatusCode code, std::string message)
        : code_(code), message_(std::move(message)) {}

    bool ok() const { return code_ == StatusCode::kOk; }
    StatusCode code() const { return code_; }
    const std::string& message() const { return message_; }

    static Status Ok() { return Status(); }

private:
    StatusCode code_;
    std::string message_;
};

}  // namespace bamsi
