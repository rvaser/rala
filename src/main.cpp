#include <stdio.h>
#include <stdlib.h>

#include "read.hpp"

using namespace RALAY;

int main(int argc, char** argv) {

    printf("Oh look is dat boi! OH SHITT WADDUP?!\n");

    auto read = createRead(0, "wtf", "hju");

    printf("%d %s %s\n", read->id(), read->sequence().c_str(), read->quality().c_str());

    return 0;
}
