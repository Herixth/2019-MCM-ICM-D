#include "Graph_def.h"

using namespace std;

const char* inFilename  = "data.txt";

int main() {
    srand(time(NULL) + rand());

    ESC_graph G(inFilename);
    
    G.init();

    G.simulate();

    system("pause");
    return 0;
}