#pragma once

namespace bamsi {

class IApp {
public:
    virtual ~IApp() = default;
    virtual int run() = 0;
};

}  // namespace bamsi
