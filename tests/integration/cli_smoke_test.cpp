#include <cstdlib>

int main() {
    int rc = std::system("./bamsi --help > /tmp/bamsi-help.txt 2>&1");
    return rc == 0 ? 0 : 1;
}
