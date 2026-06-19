#include <iostream>

extern "C" {
#include <igraph/igraph.h>
}

int main() {
    std::cout << IGRAPH_VERSION << std::endl;
    std::cout << IGRAPH_VERSION_MAJOR << "."
              << IGRAPH_VERSION_MINOR << "."
              << IGRAPH_VERSION_PATCH << std::endl;
    return 0;
}
