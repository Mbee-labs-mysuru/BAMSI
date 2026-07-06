#include <cstdlib>

int main() {
    int rc = std::system("./bamsix --help > /tmp/bamsix-help.txt 2>&1");
    return rc == 0 ? 0 : 1;
}
