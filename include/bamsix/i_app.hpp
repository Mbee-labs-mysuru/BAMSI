#pragma once

namespace bamsix {

class IApp {
public:
    virtual ~IApp() = default;
    virtual int run() = 0;
};

}  // namespace bamsix
