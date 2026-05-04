#include "bamsi/app.hpp"

int main(int argc, char** argv) {
    return bamsi::RunApp(argc, argv).ok() ? 0 : 1;
}
