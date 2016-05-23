#include "interface/input_parameters.h"
#include "dcel/mesh.h"
#include "interface/input_manager.h"
#include "computation.h"
#include "debug.h"
#include <interface/file_writer.h>

#include "docopt/docopt.h"

static const char USAGE[] =
  R"(RIVET: Rank Invariant Visualization and Exploration Tool

     The RIVET console  application computes an augmented arrangement for
     2D persistent homology, which can be visualized with the RIVET GUI app.

    Usage:
      rivet_console (-h | --help)
      rivet_console --version
      rivet_console <input_file> <output_file> [-H <dimension>] [-V <verbosity>] [-x <xbins>] [-y <ybins>]

    Options:
      -h --help                                Show this screen
      --version                                Show the version
      -H <dimension> --homology=<dimension>    Dimension of homology to compute [default: 0]
      -x <xbins> --xbins=<xbins>               Number of bins in the x direction [default: 0]
      -y <ybins> --ybins=<ybins>               Number of bins in the y direction [default: 0]
      -V <verbosity> --verbosity=<verbosity>   Verbosity level: 0 (no console output) to 10 (lots of output) [default: 2]
)";

unsigned int get_uint_or_die(std::map<std::string, docopt::value> &args, const std::string &key) {
  try {
    return static_cast<unsigned int>(args[key].asLong());
  } catch (std::exception &e) {
    std::cerr << "Argument " << key << " must be an integer" << std::endl;
    exit(1);
  }
}

int main(int argc, char *argv[])
{
    debug() << "CONSOLE RIVET" << std::endl;

    InputParameters params;   //parameter values stored here


    std::map<std::string, docopt::value> args = docopt::docopt(USAGE, {argv + 1, argv + argc}, true, "RIVET Console 0.4");

    // for (auto const &arg : args) {
    //   std::cout << arg.first << ":" << arg.second << std::endl;
    // }

    params.fileName = args["<input_file>"].asString();
    params.outputFile = args["<output_file>"].asString();
    params.dim = get_uint_or_die(args, "--homology");
    params.x_bins = get_uint_or_die(args, "--xbins");
    params.y_bins = get_uint_or_die(args, "--ybins");
    params.verbosity = get_uint_or_die(args, "--verbosity");


    InputManager inputManager(params);
    Progress progress;
    Computation computation(params, progress);
    std::shared_ptr<InputData> input = inputManager.start(progress);
    std::shared_ptr<ComputationResult> result = computation.compute(*input);
    auto arrangement = result->arrangement;
     //TESTING: print arrangement info and verify consistency
     arrangement->print_stats();
     arrangement->test_consistency();

     if(params.verbosity >= 2) { debug() << "COMPUTATION FINISHED."; }

     //if an output file has been specified, then save the arrangement
     if(!params.outputFile.empty())
     {
         std::ofstream file(params.outputFile);
         if(file.is_open())
         {
             debug() << "Writing file:" << params.outputFile;

             FileWriter fw(params, *(arrangement), input->x_exact, input->y_exact, result->xi_support);
             fw.write_augmented_arrangement(file);
         }
         else
         {
             debug() << "Error: Unable to write file:" << params.outputFile;
         }
         ///TODO: error handling?
     }
    debug() << "CONSOLE RIVET: Goodbye";
    return 0;
}

//TODO: this was copied from driver.cpp, may need to merge the two versions.
// void console_augmented_arrangement_ready(Mesh* arrangement)
// {

//     //update status
//     if(input_params.verbosity >= 2) { qDebug() << "COMPUTATION FINISHED."; }

//     //if an output file has been specified, then save the arrangement
//     if(!input_params.outputFile.isEmpty())
//     {
//         QFile file(input_params.outputFile);
//         if(file.open(QIODevice::ReadWrite | QIODevice::Truncate))
//         {
//             qDebug() << "Writing file:" << input_params.outputFile;

//             FileWriter fw(input_params, arrangement, x_exact, y_exact, xi_support);
//             fw.write_augmented_arrangement(file);
//         }
//         else
//         {
//             qDebug() << "Error: Unable to write file:" << input_params.outputFile;
//         }
//         ///TODO: error handling?
//     }

//     //must call quit when complete
//     quit();

// }//end augmented_arrangement_ready()
