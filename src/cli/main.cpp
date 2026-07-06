#include "bamsix/app.hpp"

int main(int argc, char** argv) {
    return bamsix::RunApp(argc, argv).ok() ? 0 : 1;
}
