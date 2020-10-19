#include <cstdlib>
#include <string>
#include <iostream>

#include "frontend/Frontend.h"
#include "utils/Utils.h"

int main(int argc, char **argv) {
    std::string cmd(argv[0]);
    if (argc < 2) {
        std::cerr << "Usage: " << cmd << " [mode] (mode dependent parameters). " << std::endl;
        std::cerr << "Type " << cmd << " --help to see available modes" << std::endl;
        exit(EXIT_FAILURE);
    }

    // We now shift the arguments, and pretend, that the new first (command) is "cmd mode"
    // Different modes can then parse the arguments separately and will think that the whole argv[0] is "cmd mode"
    std::string mode(argv[1]);
    std::string cmdAndMode = cmd + " " + mode;
    argv++;
    argc--;
    argv[0] = cmdAndMode.data();

    Frontend frontend(std::cout);
    if (mode == "-h" || mode == "--help") {
        frontend.printGeneralHelp(cmd);
    } else if (mode == "ed") {
        frontend.ed(argc, argv);
    } else if (mode == "analyze") {
        frontend.analyze(argc, argv);
    } else if (mode == "chebyshev") {
        frontend.chebyshev(argc, argv);
    } else if (mode == "quench") {
        frontend.quench(argc, argv);
    } else if (mode == "random_states") {
        frontend.randomStates(argc, argv);
    } else {
        die("Unknown mode " + mode);
    }

    return EXIT_SUCCESS;
}




